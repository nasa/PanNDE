/*! \headerfile MPIHaloExchanger.hpp "modules/NetMPI/include/internal/MPIHaloExchanger.hpp"
* "MPIHaloExchanger.hpp" implements efficient communication of ghost region data between
* mesh partitions in parallel simulations using MPI derived datatypes and persistent
* communication requests.
*
*/
/******************************************************
# Notices:
# Copyright 2021 United States Government as represented by the 
# Administrator of the National Aeronautics and Space Administration. 
# No copyright is claimed in the United States under Title 17, U.S. Code. 
# All Other Rights Reserved.
# 
# Googletest is a product of Google Inc. and is subject to the following:
# 
# Copyright 2008, Google Inc. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met:
#  * Redistributions of source code must retain the above copyright notice, 
#    this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright notice, 
#    this list of conditions and the following disclaimer in the documentation 
#    and/or other materials provided with the distribution.
#  * Neither the name of Google Inc. nor the names of its contributors may be used 
#    to endorse or promote products derived from this software without specific 
#    prior written permission.
#
# Disclaimers
# No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, 
# EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT 
# THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT 
# SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO 
# THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY 
# GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE 
# PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT 
# AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE 
# ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS." 
#
# Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES 
# GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S 
# USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES 
# ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S 
# USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, 
# ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  
# RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS 
# AGREEMENT.
******************************************************/

#pragma once

#include <vector>
#include <cstdint>
#include <set>
#include <unordered_map>
#include <string>

#include <mpi.h>

#include "MPIType.hpp"

#include "Mesh.hpp"
#include "Array.hpp"
#include "Field.hpp"

namespace NetMPI {
  /*! \class MPIHaloExchanger MPIHaloExchanger.hpp "modules/NetMPI/include/internal/MPIHaloExchanger.hpp"
  *
  * Implements ghost region (halo) exchange between partitioned meshes in parallel simulations.
  * 
  * In parallel mesh-based simulations, each process owns a portion (partition) of the overall mesh.
  * However, to accurately compute values near partition boundaries, processes need data from 
  * adjacent partitions. These shared boundary regions are called "ghost regions" or "halos."
  * 
  * MPIHaloExchanger analyzes mesh connectivity to automatically determine which nodes need to be
  * exchanged between processes. It then creates optimized MPI derived datatypes to represent
  * these non-contiguous data regions and establishes persistent communication patterns for
  * efficient updates during the simulation.
  * 
  * This class provides:
  * 1. Automatic detection of partition boundaries requiring communication
  * 2. Creation of optimized communication patterns using MPI datatypes
  * 3. Setup of persistent communication requests for field data
  * 4. Asynchronous initiation and completion of halo exchanges
  *
  */
  class MPIHaloExchanger {
    public:
      /*!
      * Constructs a halo exchanger with a specified MPI communicator.
      * 
      * \param comm MPI_Comm The MPI communicator to use (default: MPI_COMM_WORLD)
      */
      MPIHaloExchanger(MPI_Comm comm=MPI_COMM_WORLD){
        setupMPI(comm);
      };

      /*!
      * Determines the halo exchange pattern based on mesh connectivity.
      * 
      * This method analyzes the mesh to identify which nodes need to be sent to
      * and received from other processes. It then creates MPI derived datatypes
      * to represent these non-contiguous data regions for efficient communication.
      * 
      * \param mesh std::shared_ptr<PanNDE::Mesh> The partitioned mesh to analyze
      */
      void determineExchangePattern(std::shared_ptr<PanNDE::Mesh> mesh){
        //build block starts and lengths. applies to all state node vectors
        //  place into the mpi datatype map
        (this->mesh)=mesh;
        
        //logic: loop over all cells, if any non-local ownership nodes attached,
        //  indicate those nodes for recv from their rank, and all local-owned nodes
        //  for send to that rank

        std::unordered_map<int32_t/*send to rank*/,
                           std::set<int32_t>/*nodes to send*/> sendSet;
        std::unordered_map<int32_t/*recv from rank*/,
                           std::set<int32_t>/*nodes to recv*/> recvSet;

        parseHaloNodesToCommSets(sendSet,recvSet);

        for(auto it=sendSet.begin();it!=sendSet.end();it++){
          sendhalos.emplace(it->first,new MPI_Datatype());
          printf("[%i] send %lu nodes to rank %i\n",
                  rank,it->second.size(),it->first);
          makeMPIHaloExchangeDataTypes(it->second,sendhalos.at(it->first));
        };
        for(auto it=recvSet.begin();it!=recvSet.end();it++){
          recvhalos.emplace(it->first,new MPI_Datatype());
          printf("[%i] recv %lu nodes from rank %i\n",
                  rank,it->second.size(),it->first);
          makeMPIHaloExchangeDataTypes(it->second,recvhalos.at(it->first));
        };
      };

      /*!
      * Registers a field for halo exchange using the previously determined pattern.
      * 
      * \param keyname std::string Name identifier for the field
      * \param field std::shared_ptr<PanNDE::Field> The field to register
      */
      void setupDataLinks(std::string keyname,std::shared_ptr<PanNDE::Field> field){
        exchanging_fields.emplace(keyname,field);
        if(PanNDE::Field::NODE==field->type()){
          setupNodeDataLinks(keyname,field);
        };
        //if cell comms become a thing, do it
      };

      /*!
      * Sets up halo exchange for all fields in a field bundle.
      * 
      * This is the preferred method for configuring halo exchanges as it
      * handles all fields in the bundle at once.
      * 
      * \param fields std::shared_ptr<PanNDE::FieldBundle> Bundle containing fields to register
      */
      void setupDataLinks(std::shared_ptr<PanNDE::FieldBundle> fields){
        determineExchangePattern(fields->mesh());
        for(int k=0;k<fields->fieldCount();k++){
          auto name=fields->fieldName(k);
          setupDataLinks(name,fields->field(name));
        };
      };
      
      /*!
      * Initiates asynchronous halo exchange for a specific field.
      * 
      * This starts the non-blocking communication operations to update
      * ghost regions for the named field but does not wait for completion.
      * 
      * \param keyname std::string Name of the field to update
      */
      void startHaloExchange(std::string keyname){
        auto send_targets=send_halo_requests.equal_range(keyname);
        for(auto it=send_targets.first;it!=send_targets.second;it++){
          MPI_Start(&(it->second.request));
        };
        auto recv_targets=recv_halo_requests.equal_range(keyname);
        for(auto it=recv_targets.first;it!=recv_targets.second;it++){
          MPI_Start(&(it->second.request));
        };
      };

      /*!
      * Waits for completion of a previously initiated halo exchange.
      * 
      * This blocks until the asynchronous halo exchange for the named field
      * has completed, ensuring that ghost regions are fully updated.
      * 
      * \param keyname std::string Name of the field being updated
      */
      void waitUntilDone(std::string keyname){
        auto send_targets=send_halo_requests.equal_range(keyname);
        for(auto it=send_targets.first;it!=send_targets.second;it++){
          MPI_Wait(&(it->second.request),MPI_STATUS_IGNORE);
        };
        auto recv_targets=recv_halo_requests.equal_range(keyname);
        for(auto it=recv_targets.first;it!=recv_targets.second;it++){
          MPI_Wait(&(it->second.request),MPI_STATUS_IGNORE);
        };
      };

    private:
      /*!
      * Structure representing a persistent communication request.
      */
      struct CommRequest{
        //! Default constructor
        CommRequest(){};
        
        /*!
        * Constructor with remote rank and tag.
        * \param remote int The remote process rank
        * \param tag int The communication tag
        */
        CommRequest(int remote,int tag){remote_rank=remote;comm_tag=tag;};
        
        //! Rank of the remote process
        int remote_rank;
        //! Communication tag
        int comm_tag;
        //! MPI request handle
        MPI_Request request;
      };

      /*!
      * Sets up the MPI environment and request tracking structures.
      * \param comm MPI_Comm The MPI communicator to use
      */
      void setupMPI(MPI_Comm comm){
        MPI_Comm_dup(comm,&halo_comm);
        MPI_Comm_rank(halo_comm,&rank);
        MPI_Comm_size(halo_comm,&comm_size);
        setupCounters();
      };

      /*!
      * Analyzes mesh cells to determine which nodes need to be exchanged.
      * \param sendSet std::unordered_map<int32_t,std::set<int32_t>>& Map of nodes to send by destination rank
      * \param recvSet std::unordered_map<int32_t,std::set<int32_t>>& Map of nodes to receive by source rank
      */
      void parseHaloNodesToCommSets(std::unordered_map<int32_t,std::set<int32_t>>& sendSet,
                                    std::unordered_map<int32_t,std::set<int32_t>>& recvSet){
        int32_t box[8];
        for(int kc=0;kc<mesh->cellCount();kc++){
          mesh->cell(kc,box);
          parseHaloNodesOnCellToCommSets(box,sendSet,recvSet);
        };
      };

      /*!
      * Analyzes a single cell to determine which nodes need to be exchanged.
      * \param box int32_t[8] Array of node indices for the cell
      * \param sendSet std::unordered_map<int32_t,std::set<int32_t>>& Map of nodes to send
      * \param recvSet std::unordered_map<int32_t,std::set<int32_t>>& Map of nodes to receive
      */
      void parseHaloNodesOnCellToCommSets(int32_t box[8],std::unordered_map<int32_t,std::set<int32_t>>& sendSet,
                                          std::unordered_map<int32_t,std::set<int32_t>>& recvSet){
        int32_t ownership[8];
        getCellNodeOwnership(ownership,box);
        if(isCellHalo(ownership)){
          for(int kb=0;kb<8;kb++){
            placeNodeForComms(kb,ownership,box,sendSet,recvSet);
          };
        };
      };

      /*!
      * Gets the owner partition for each node in a cell.
      * \param ownership int32_t[8] Output array to receive owner partition IDs
      * \param box int32_t[8] Array of node indices for the cell
      */
      void getCellNodeOwnership(int32_t ownership[8],int32_t box[8]){
        for(int kb=0;kb<8;kb++){ownership[kb]=mesh->nodeHomePartition(box[kb]);};
      };

      /*!
      * Determines if a cell spans partition boundaries.
      * \param ownership int32_t[8] Array of owner partition IDs for cell nodes
      * \return bool True if the cell has nodes owned by different partitions
      */
      bool isCellHalo(int32_t ownership[8]){
        for(int kb=0;kb<8;kb++){
          if(ownership[kb]!=mesh->partitionId()){return true;};
        };
        return false;
      };

      /*!
      * Records a node for sending or receiving based on ownership pattern.
      * \param kba int Index of the node within the cell
      * \param ownership int32_t[8] Array of owner partition IDs for cell nodes
      * \param box int32_t[8] Array of node indices for the cell
      * \param sendSet std::unordered_map<int32_t,std::set<int32_t>>& Map of nodes to send
      * \param recvSet std::unordered_map<int32_t,std::set<int32_t>>& Map of nodes to receive
      */
      void placeNodeForComms(int kba,int32_t ownership[8],int32_t box[8],
                             std::unordered_map<int32_t,std::set<int32_t>>& sendSet,
                             std::unordered_map<int32_t,std::set<int32_t>>& recvSet){
        if(ownership[kba]!=mesh->partitionId()){
          auto emplr=recvSet.emplace(ownership[kba],std::set<int32_t>());
          emplr.first->second.emplace(box[kba]);
        }else{
          for(int kbc=0;kbc<8;kbc++){
            if(ownership[kbc]!=mesh->partitionId()){
              auto empls=sendSet.emplace(ownership[kbc],std::set<int32_t>());
              empls.first->second.emplace(box[kba]);
            };
          };
        };
      };

      /*!
      * Creates an MPI derived datatype for non-contiguous halo regions.
      * \param nodeset std::set<int32_t> Set of node indices to include in the datatype
      * \param haloType MPI_Datatype* Output parameter to receive the created datatype
      */
      void makeMPIHaloExchangeDataTypes(std::set<int32_t> nodeset,MPI_Datatype* haloType){
        std::vector<int32_t> block_offsets;
        std::vector<int32_t> block_sizes;
        makeOffsetsAndSizes(nodeset,block_offsets,block_sizes);
        MPI_Type_indexed(block_offsets.size(),block_sizes.data(),block_offsets.data(),MPI_DOUBLE,haloType);
        MPI_Type_commit(haloType);
      };

      /*!
      * Prepares arrays of offsets and sizes for creating MPI derived datatypes.
      * \param nodeset std::set<int32_t> Set of node indices
      * \param offsets std::vector<int32_t>& Output vector of offsets
      * \param sizes std::vector<int32_t>& Output vector of block sizes
      */
      void makeOffsetsAndSizes(std::set<int32_t> nodeset,std::vector<int32_t>& offsets,std::vector<int32_t>& sizes){
        auto set_it=nodeset.begin();
        offsets.resize(0);
        sizes.resize(0);
        while(set_it!=nodeset.end()){
          offsets.push_back(*set_it);
          sizes.push_back(1);
          //while(*(std::next(set_it))==(1+(*set_it))){
            set_it++;
          //  sizes.back()++;
          //};
        };
      };

      /*!
      * Initializes communication counters for all processes.
      */
      void setupCounters(){
        send_comm_counter.resize(0);send_comm_counter.resize(comm_size,0);
        recv_comm_counter.resize(0);recv_comm_counter.resize(comm_size,0);
      };

      /*!
      * Sets up persistent communication requests for node field data.
      * \param keyname std::string Name of the field
      * \param field std::shared_ptr<PanNDE::Field> Field to register for halo exchange
      */
      void setupNodeDataLinks(std::string keyname,std::shared_ptr<PanNDE::Field> field){
        int remoteRank;
        for(auto it=sendhalos.begin();it!=sendhalos.end();it++){
          remoteRank=it->first;
          auto dtype=it->second;
          auto send_req=createSendHaloTag(keyname,remoteRank);
          MPI_Send_init(field->data(),1,*dtype,send_req->remote_rank,send_req->comm_tag,
                        halo_comm,&(send_req->request));
        };
        for(auto it=recvhalos.begin();it!=recvhalos.end();it++){
          remoteRank=it->first;
          auto dtype=it->second;
          auto recv_req=createRecvHaloTag(keyname,remoteRank);
          MPI_Recv_init(field->data(),1,*dtype,recv_req->remote_rank,recv_req->comm_tag,
                        halo_comm,&(recv_req->request));
        };
      };

      /*!
      * Creates a persistent send request for a field to a specific process.
      * \param keyname std::string Field name
      * \param remote_rank int Destination process rank
      * \return CommRequest* Pointer to the created request object
      */
      CommRequest* createSendHaloTag(std::string keyname,int remote_rank){
        int tag=send_comm_counter.at(remote_rank)*comm_size+remote_rank;
        send_comm_counter.at(remote_rank)++;
        auto empl=send_halo_requests.emplace(keyname,CommRequest(remote_rank,tag));
        return &(empl->second);
      };

      /*!
      * Creates a persistent receive request for a field from a specific process.
      * \param keyname std::string Field name
      * \param remote_rank int Source process rank
      * \return CommRequest* Pointer to the created request object
      */
      CommRequest* createRecvHaloTag(std::string keyname,int remote_rank){
        int tag=recv_comm_counter.at(remote_rank)*comm_size+rank;
        recv_comm_counter.at(remote_rank)++;
        auto empl=recv_halo_requests.emplace(keyname,CommRequest(remote_rank,tag));
        return &(empl->second);
      };

      //! The mesh for which halo exchanges are being performed
      std::shared_ptr<PanNDE::Mesh> mesh=nullptr;

      //! Counters for generating unique receive tags for each remote process
      std::vector<int> recv_comm_counter;
      //! Counters for generating unique send tags for each remote process
      std::vector<int> send_comm_counter;

      //! Map of MPI datatypes for sending halo data to each process
      std::unordered_map<int32_t,MPI_Datatype*> sendhalos;
      //! Map of MPI datatypes for receiving halo data from each process
      std::unordered_map<int32_t,MPI_Datatype*> recvhalos;
      //! Map of fields registered for halo exchange
      std::map<std::string,std::shared_ptr<PanNDE::Field>> exchanging_fields;

      //! Multimap of persistent send requests for each field
      std::multimap<std::string/*keyname*/,CommRequest> send_halo_requests;
      //! Multimap of persistent receive requests for each field
      std::multimap<std::string/*keyname*/,CommRequest> recv_halo_requests;

      //! Process rank within the communicator
      int rank;
      //! Total number of processes in the communicator
      int comm_size;
      //! The MPI communicator used for halo exchanges
      MPI_Comm halo_comm;
  };
};