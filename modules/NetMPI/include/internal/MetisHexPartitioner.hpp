/*! \headerfile MetisHexPartitioner.hpp "modules/NetMPI/internal/MetisHexPartitioner.hpp"
* "MetisHexPartitioner.hpp" contains a class implementation for partitioning 
* hexahedral meshes using the METIS library for domain decomposition.
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

#include <cstdint>
#include <vector>
#include <set>
#include <memory>

#include "metis.h"

#include "Mesh.hpp"

namespace NetMPI {
  /*! \class MetisHexPartitioner MetisHexPartitioner.hpp "modules/NetMPI/internal/MetisHexPartitioner.hpp"
  *
  * Implements mesh partitioning for hexahedral meshes using the METIS graph partitioning library.
  * This class is responsible for decomposing a global mesh into distributed partitions for
  * parallel processing, tracking both node and cell ownership.
  *
  */
  class MetisHexPartitioner {
    public:
      /*!
      * Partitions a global mesh into the specified number of partitions using METIS.
      * \param global_mesh std::shared_ptr<PanNDE::Mesh>& The global mesh to be partitioned
      * \param Npartitions int The number of partitions to create
      */
      void partition(std::shared_ptr<PanNDE::Mesh>& global_mesh,int Npartitions){
        setup(global_mesh,Npartitions);
        doPartition();
      };

      /*!
      * Retrieves the set of node indices that belong to the specified partition.
      * \param partitionId int The ID of the partition to query
      * \return std::set<int64_t>* Pointer to the set of node indices in the partition
      */
      std::set<int64_t>* getPartitionNodes(int partitionId){
        return &(nodesOnPartition.at(partitionId));
      };

      /*!
      * Retrieves the set of cell indices that belong to the specified partition.
      * \param partitionId int The ID of the partition to query
      * \return std::set<int64_t>* Pointer to the set of cell indices in the partition
      */
      std::set<int64_t>* getPartitionCells(int partitionId){
        return &(cellsOnPartition.at(partitionId));
      };

      /*!
      * Gets the ownership information for all nodes in the domain.
      * \return std::vector<int32_t>* Pointer to vector containing partition ID for each node
      */
      std::vector<int32_t>* getFullDomainNodeOwnership(){
        return &nodeOwnership;
      };

      /*!
      * Gets the ownership information for all cells in the domain.
      * \return std::vector<int32_t>* Pointer to vector containing partition ID for each cell
      */
      std::vector<int32_t>* getFullDomainCellOwnership(){
        return &cellOwnership;
      };

    private:
      //! Reference to the global mesh being partitioned
      std::shared_ptr<PanNDE::Mesh> gmesh;

      //! Vector storing partition ownership for each cell (index is cell ID, value is partition ID)
      std::vector<int32_t> cellOwnership;
      
      //! Vector storing partition ownership for each node (index is node ID, value is partition ID)
      std::vector<int32_t> nodeOwnership;

      //! Collection of sets containing node indices for each partition
      std::vector<std::set<int64_t>> nodesOnPartition;
      
      //! Collection of sets containing cell indices for each partition
      std::vector<std::set<int64_t>> cellsOnPartition;

      //! Number of partitions to create
      int Nparts=0;

      /*!
      * Sets up internal data structures for partitioning.
      * \param mesh std::shared_ptr<PanNDE::Mesh> The global mesh to partition
      * \param Npartitions int The number of partitions to create
      */
      void setup(std::shared_ptr<PanNDE::Mesh> mesh,int Npartitions){
        gmesh=mesh;
        Nparts=Npartitions;
        nodesOnPartition.resize(Nparts);
        cellsOnPartition.resize(Nparts);
        //partitionNodeHaloOutput.resize(Nparts);
        //partitionCellsOutput.resize(Nparts);
      };

      /*!
      * Executes the partitioning process.
      */
      void doPartition(){
        executeMetis();
        parseMetis();
      };

      /*!
      * Processes partition results from METIS into internal data structures.
      */
      void parseMetis(){
        sortNodes();
        getCells();
        getNodes();
      };

      /*!
      * Organizes nodes into partition-specific sets based on ownership information.
      */
      void sortNodes(){
        for(int64_t kn=0;kn<nodeOwnership.size();kn++){
          nodesOnPartition.at(nodeOwnership.at(kn)).emplace(kn);
        };
      };

      /*!
      * Populates cell sets for all partitions.
      */
      void getCells(){
        for(int64_t kp=0;kp<Nparts;kp++){
          getCellsOnPartition(kp);
        };
      };

      /*!
      * Determines which cells belong to a specific partition based on node ownership.
      * \param partitionId int64_t The partition ID to process
      */
      void getCellsOnPartition(int64_t partitionId){
        int32_t box[8];
        for(auto it=nodesOnPartition.at(partitionId).begin();
                 it!=nodesOnPartition.at(partitionId).end();it++){
          gmesh->connectedCells(*it,box);
          putCellsIntoSet(partitionId,box);
        };
      };

      /*!
      * Adds cell indices to the specified partition's cell set.
      * \param partitionId int64_t The partition ID to add cells to
      * \param box int32_t[8] Array of cell indices to add
      */
      void putCellsIntoSet(int64_t partitionId,int32_t box[8]){
        for(int kb=0;kb<8;kb++){
          if(-1!=box[kb]){cellsOnPartition.at(partitionId).emplace(box[kb]);};
        };
      };

      /*!
      * Populates node sets for all partitions.
      */
      void getNodes(){
        for(int64_t kp=0;kp<Nparts;kp++){
          getNodesOnPartition(kp);
        };
      };

      /*!
      * Determines which nodes belong to a specific partition based on cell ownership.
      * \param partitionId int64_t The partition ID to process
      */
      void getNodesOnPartition(int64_t partitionId){
        int32_t box[8];
        for(auto it=cellsOnPartition.at(partitionId).begin();
                 it!=cellsOnPartition.at(partitionId).end();it++){
          gmesh->cell(*it,box);
          putNodesIntoSet(partitionId,box);
        };
      };

      /*!
      * Adds node indices to the specified partition's node set.
      * \param partitionId int64_t The partition ID to add nodes to
      * \param box int32_t[8] Array of node indices to add
      */
      void putNodesIntoSet(int64_t partitionId,int32_t box[8]){
        for(int kb=0;kb<8;kb++){
          nodesOnPartition.at(partitionId).emplace(box[kb]);
        };
      };

      /*!
      * Calls the METIS library to perform the actual partitioning.
      * Sets up the required data structures for METIS, calls the partitioning
      * routine, and stores the results in the ownership vectors.
      */
      void executeMetis(){
        idx_t ne=gmesh->cellCount();
        idx_t nn=gmesh->nodeCount();
        std::vector<idx_t> eptr(ne+1);
        std::vector<idx_t> eind(8*ne);
        idx_t ncommon=1;
        idx_t nparts=Nparts;
        idx_t objval;
        std::vector<idx_t> epart(ne);
        std::vector<idx_t> npart(nn);
        nodeOwnership.resize(nn);
        cellOwnership.resize(ne);

        int32_t box[8];
        for(int k=0;k<ne;k++){
          eptr.at(k)=8*k;
          gmesh->cell(k,box);
          for(int kb=0;kb<8;kb++){
            eind.at(8*k+kb)=box[kb];
          };
        };
        eptr.at(ne)=8*ne;

        //partition
        METIS_PartMeshDual(&ne,&nn,eptr.data(),eind.data(),
                           NULL,NULL,&ncommon,&nparts,
                           NULL,NULL,&objval,
                           epart.data(),npart.data());
        
        for(int k=0;k<nn;k++){
          nodeOwnership.at(k)=(int32_t)(npart.at(k));
        };
        for(int k=0;k<ne;k++){
          cellOwnership.at(k)=(int32_t)(epart.at(k));
        };
      };
  };
};