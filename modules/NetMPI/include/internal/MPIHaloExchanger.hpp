//
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

#include <mpi.h>

#include "MPIType.hpp"

#include "Mesh.hpp"
#include "Array.hpp"
#include "Field.hpp"

namespace NetMPI {
  class MPIHaloExchanger {
    public:
      MPIHaloExchanger(MPI_Comm comm=MPI_COMM_WORLD){
        setupMPI(comm);
      };
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
          printf("[%i] send %i nodes to rank %i\n",
                  rank,it->second.size(),it->first);
          makeMPIHaloExchangeDataTypes(it->second,sendhalos.at(it->first));
        };
        for(auto it=recvSet.begin();it!=recvSet.end();it++){
          recvhalos.emplace(it->first,new MPI_Datatype());
          printf("[%i] recv %i nodes from rank %i\n",
                  rank,it->second.size(),it->first);
          makeMPIHaloExchangeDataTypes(it->second,recvhalos.at(it->first));
        };
      };
      void setupDataLinks(std::string keyname,std::shared_ptr<PanNDE::Field> field){
        exchanging_fields.emplace(keyname,field);
        if(PanNDE::Field::NODE==field->type()){
          setupNodeDataLinks(keyname,field);
        };
        //if cell comms become a thing, do it
      };
      //prefered halo setup
      void setupDataLinks(std::shared_ptr<PanNDE::FieldBundle> fields){
        determineExchangePattern(fields->mesh());
        for(int k=0;k<fields->fieldCount();k++){
          auto name=fields->fieldName(k);
          setupDataLinks(name,fields->field(name));
        };
      };
      
      //directly testable
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
      struct CommRequest{
        CommRequest(){};
        CommRequest(int remote,int tag){remote_rank=remote;comm_tag=tag;};
        int remote_rank;
        int comm_tag;
        MPI_Request request;
      };

      void setupMPI(MPI_Comm comm){
        MPI_Comm_dup(comm,&halo_comm);
        MPI_Comm_rank(halo_comm,&rank);
        MPI_Comm_size(halo_comm,&comm_size);
        setupCounters();
      };
      void parseHaloNodesToCommSets(std::unordered_map<int32_t,std::set<int32_t>>& sendSet,
                                    std::unordered_map<int32_t,std::set<int32_t>>& recvSet){
        int32_t box[8];
        for(int kc=0;kc<mesh->cellCount();kc++){
          mesh->cell(kc,box);
          parseHaloNodesOnCellToCommSets(box,sendSet,recvSet);
        };
      };
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
      void getCellNodeOwnership(int32_t ownership[8],int32_t box[8]){
        for(int kb=0;kb<8;kb++){ownership[kb]=mesh->nodeHomePartition(box[kb]);};
      };
      bool isCellHalo(int32_t ownership[8]){
        for(int kb=0;kb<8;kb++){
          if(ownership[kb]!=mesh->partitionId()){return true;};
        };
        return false;
      };
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
      void makeMPIHaloExchangeDataTypes(std::set<int32_t> nodeset,MPI_Datatype* haloType){
        std::vector<int32_t> block_offsets;
        std::vector<int32_t> block_sizes;
        makeOffsetsAndSizes(nodeset,block_offsets,block_sizes);
        MPI_Type_indexed(block_offsets.size(),block_sizes.data(),block_offsets.data(),MPI_DOUBLE,haloType);
        MPI_Type_commit(haloType);
      };
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
      void setupCounters(){
        send_comm_counter.resize(0);send_comm_counter.resize(comm_size,0);
        recv_comm_counter.resize(0);recv_comm_counter.resize(comm_size,0);
      };
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

      CommRequest* createSendHaloTag(std::string keyname,int remote_rank){
        int tag=send_comm_counter.at(remote_rank)*comm_size+remote_rank;
        send_comm_counter.at(remote_rank)++;
        auto empl=send_halo_requests.emplace(keyname,CommRequest(remote_rank,tag));
        return &(empl->second);
      };
      CommRequest* createRecvHaloTag(std::string keyname,int remote_rank){
        int tag=recv_comm_counter.at(remote_rank)*comm_size+rank;
        recv_comm_counter.at(remote_rank)++;
        auto empl=recv_halo_requests.emplace(keyname,CommRequest(remote_rank,tag));
        return &(empl->second);
      };

      std::shared_ptr<PanNDE::Mesh> mesh=nullptr;

      std::vector<int> recv_comm_counter;
      std::vector<int> send_comm_counter;

      std::unordered_map<int32_t,MPI_Datatype*> sendhalos;
      std::unordered_map<int32_t,MPI_Datatype*> recvhalos;
      std::map<std::string,std::shared_ptr<PanNDE::Field>> exchanging_fields;

      std::multimap<std::string/*keyname*/,CommRequest> send_halo_requests;
      std::multimap<std::string/*keyname*/,CommRequest> recv_halo_requests;

      int rank,comm_size;
      MPI_Comm halo_comm;
  };
};