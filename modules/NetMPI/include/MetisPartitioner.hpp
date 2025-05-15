/*! \headerfile MetisPartitioner.hpp "modules/NetMPI/include/MetisPartitioner.hpp"
* "MetisPartitioner.hpp" implements domain decomposition and parallel distribution
* of meshes and associated fields using the METIS graph partitioning library.
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

#include <memory>
#include <map>

#include "Mesh.hpp"
#include "Field.hpp"
#include "Array.hpp"

#include "Communicator.hpp"
#include "MPICommunicator.hpp"

#include "Partitioner.hpp"

#include "internal/MetisHexPartitioner.hpp"

namespace NetMPI {
  /*! \class MetisPartitioner MetisPartitioner.hpp "modules/NetMPI/include/MetisPartitioner.hpp"
  *
  * Implements domain decomposition and distribution for parallel simulations.
  * 
  * MetisPartitioner uses the METIS graph partitioning library to divide a single
  * global mesh into multiple partitions optimized for parallel computing. It handles
  * both the mesh partitioning process and the distribution of mesh partitions and
  * associated field data to appropriate MPI processes.
  * 
  * The partitioner is designed to:
  * 1. Minimize communication by optimizing partition boundaries
  * 2. Balance computational load across processes
  * 3. Preserve mesh connectivity relationships across partition boundaries
  * 4. Manage proper data distribution for mesh-associated fields
  * 
  * This class implements the Controller::Partitioner interface and serves as
  * the bridge between PanNDE's abstract domain decomposition requirements and
  * the specific METIS implementation details.
  *
  */
  class MetisPartitioner : public Controller::Partitioner {
    public:
      /*!
      * Constructs a mesh partitioner with necessary factory objects.
      * 
      * The factories are used to create the various data arrays needed during
      * the partitioning and distribution process.
      * 
      * \param i32_maker std::shared_ptr<PanNDE::ArrayFactory<int32_t>> Factory for 32-bit integer arrays
      * \param i64_maker std::shared_ptr<PanNDE::ArrayFactory<int64_t>> Factory for 64-bit integer arrays
      * \param dbl_maker std::shared_ptr<PanNDE::ArrayFactory<double>> Factory for double-precision arrays
      * \param communicator std::shared_ptr<PanNDE::Communicator> Communicator for inter-process data exchange (default: creates a new MPICommunicator)
      */
      MetisPartitioner(std::shared_ptr<PanNDE::ArrayFactory<int32_t>> i32_maker,
                       std::shared_ptr<PanNDE::ArrayFactory<int64_t>> i64_maker,
                       std::shared_ptr<PanNDE::ArrayFactory<double>> dbl_maker,
                       std::shared_ptr<PanNDE::Communicator>communicator=nullptr){
        configure(communicator);
        i32mfg=i32_maker;
        i64mfg=i64_maker;
        dblmfg=dbl_maker;
      };

      /*!
      * Creates a shared pointer to a MetisPartitioner.
      * 
      * \param i32_maker std::shared_ptr<PanNDE::ArrayFactory<int32_t>> Factory for 32-bit integer arrays
      * \param i64_maker std::shared_ptr<PanNDE::ArrayFactory<int64_t>> Factory for 64-bit integer arrays
      * \param dbl_maker std::shared_ptr<PanNDE::ArrayFactory<double>> Factory for double-precision arrays
      * \param communicator std::shared_ptr<PanNDE::Communicator> Communicator for inter-process data exchange (default: creates a new MPICommunicator)
      * \return std::shared_ptr<NetMPI::MetisPartitioner> Shared pointer to a new MetisPartitioner
      */
      static
      std::shared_ptr<NetMPI::MetisPartitioner> makeShared(
                      std::shared_ptr<PanNDE::ArrayFactory<int32_t>> i32_maker,
                      std::shared_ptr<PanNDE::ArrayFactory<int64_t>> i64_maker,
                      std::shared_ptr<PanNDE::ArrayFactory<double>> dbl_maker,
                      std::shared_ptr<PanNDE::Communicator>communicator=nullptr){
        auto partitioner=std::make_shared<NetMPI::MetisPartitioner>(NetMPI::MetisPartitioner(
                            i32_maker,i64_maker,dbl_maker,communicator));
        return std::move(partitioner);
      };

      /*!
      * Partitions a global mesh into multiple domains.
      * 
      * This method divides the provided mesh into the specified number of partitions
      * using the METIS graph partitioning algorithm. The partitioning information is
      * stored internally and used for subsequent distribution operations.
      * 
      * \param Nparts int Number of partitions to create
      * \param global_mesh std::shared_ptr<PanNDE::Mesh> The complete mesh to be partitioned
      * \throws std::runtime_error If the provided mesh is null
      */
      void partitionMesh(int Nparts,std::shared_ptr<PanNDE::Mesh> global_mesh) override {
        if(nullptr==global_mesh){invalidMeshError();};
        gmesh=global_mesh;
        partition_count=Nparts;
        goPartition();
      };

      /*!
      * Distributes mesh partitions to requesting processes.
      * 
      * After partitioning, this method handles the distribution of mesh partitions
      * to the appropriate processes. Each process specifies which partition IDs it
      * wants to receive, and the process containing the global mesh sends the
      * corresponding partitions.
      * 
      * \param partitionIds std::shared_ptr<PanNDE::Array<int>> Partition IDs requested by this process
      * \param maker std::shared_ptr<PanNDE::MeshFactory> Factory for creating mesh objects
      * \param mesh_source_rank int Rank of the process containing the global mesh (default: 0)
      * \return std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> Array of local mesh partitions
      */
      std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> distributeMeshPartitions(
                                  std::shared_ptr<PanNDE::Array<int>> partitionIds,
                                  std::shared_ptr<PanNDE::MeshFactory> maker,
                                  int mesh_source_rank=0) override {
        copyIds(partitionIds);
        if(mesh_source_rank==rank){//send
          for(int kr=0;kr<nranks;kr++){
            if(kr==rank){continue;};
            sendPartitions(kr);
          };
          local_mesh_array=makeLocalPartitions(maker);
        }else{//recv
          local_mesh_array=recvPartitions(maker,mesh_source_rank);
        };
        return local_mesh_array;
      };

      /*!
      * Distributes field partitions based on the previously distributed mesh.
      * 
      * This method distributes a field across the partitioned mesh, ensuring that
      * field values are correctly mapped to the corresponding mesh elements in each
      * partition. The mesh must have been distributed before calling this method.
      * 
      * \param partitionIds std::shared_ptr<PanNDE::Array<int>> Partition IDs requested by this process
      * \param maker std::shared_ptr<PanNDE::FieldFactory> Factory for creating field objects
      * \param global_field std::shared_ptr<PanNDE::Field> The global field to distribute (required on source rank)
      * \param field_source_rank int Rank of the process containing the global field (default: 0)
      * \return std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Field>>> Array of local field partitions
      * \throws std::runtime_error If mesh has not been distributed or if field is invalid
      */
      std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Field>>> distributeFieldPartitions(
                                  std::shared_ptr<PanNDE::Array<int>> partitionIds,
                                  std::shared_ptr<PanNDE::FieldFactory> maker,
                                  std::shared_ptr<PanNDE::Field> global_field=nullptr,
                                  int field_source_rank=0) override {
        checkFieldDistributionReadiness(global_field,field_source_rank);
        //distribute from partitioning rank
        std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Field>>> local_field_array=nullptr;
        if(field_source_rank==rank){//send
          for(int kr=0;kr<nranks;kr++){
            if(kr==rank){continue;};
            sendFieldPartitions(kr,global_field);
          };
          local_field_array=makeLocalFieldPartitions(maker,global_field);
        }else{//recv
          local_field_array=recvFieldPartitions(maker,field_source_rank);
        };
        return local_field_array;
      };

    private:
      /*!
      * Copies partition IDs from the input array to an internal storage array.
      * \param part_ids std::shared_ptr<PanNDE::Array<int32_t>> Source array of partition IDs
      */
      void copyIds(std::shared_ptr<PanNDE::Array<int32_t>> part_ids){
        local_partition_ids=i32mfg->makeManagedArray();
        local_partition_ids->resize(part_ids->size());
        for(int k=0;k<part_ids->size();k++){local_partition_ids->at(k)=part_ids->at(k);};
      };

      /*!
      * Verifies that conditions are met for field distribution.
      * \param global_field std::shared_ptr<PanNDE::Field> The global field to distribute
      * \param src_rank int Source rank for the field distribution
      * \throws std::runtime_error If prerequisites are not met
      */
      void checkFieldDistributionReadiness(std::shared_ptr<PanNDE::Field> global_field,int src_rank){
        if(nullptr==local_mesh_array){throw std::runtime_error("mesh must be distributed");};
        if(rank==src_rank){
          if(nullptr==global_field){throw std::runtime_error("no mesh provided");};
          if(PanNDE::Field::NODE!=global_field->type() && PanNDE::Field::CELL!=global_field->type()){
            throw std::runtime_error("invalid field type for redistribution");
          };
        };
      };

      /*!
      * Creates field partitions for local use on the source process.
      * \param fmaker std::shared_ptr<PanNDE::FieldFactory> Factory for creating fields
      * \param gfield std::shared_ptr<PanNDE::Field> The global field to partition
      * \return std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Field>>> Array of local field partitions
      */
      std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Field>>> makeLocalFieldPartitions(
                            std::shared_ptr<PanNDE::FieldFactory> fmaker,
                            std::shared_ptr<PanNDE::Field> gfield){
        auto field_array=fmaker->makeManagedFieldArray();
        field_array->resize(local_partition_ids->size());
        for(int kp=0;kp<local_partition_ids->size();kp++){
          field_array->at(kp)=makeLocalFieldPartition(kp,fmaker,gfield);
        };
        return std::move(field_array);
      };

      /*!
      * Creates a single field partition for local use on the source process.
      * \param idx int Index of the partition in the local array
      * \param fmaker std::shared_ptr<PanNDE::FieldFactory> Factory for creating fields
      * \param gfield std::shared_ptr<PanNDE::Field> The global field
      * \return std::shared_ptr<PanNDE::Field> The partitioned field
      * \throws std::runtime_error If partition size mismatches occur
      */
      std::shared_ptr<PanNDE::Field> makeLocalFieldPartition(int idx,
                                std::shared_ptr<PanNDE::FieldFactory> fmaker,
                                std::shared_ptr<PanNDE::Field> gfield){
        auto field=fmaker->makeManagedField(local_mesh_array->at(idx),gfield->type());
        auto gids=(PanNDE::Field::NODE==gfield->type())?
                            getNodeGlobalIds(local_partition_ids->at(idx)):
                            getCellGlobalIds(local_partition_ids->at(idx));
        if(gids->size()!=field->size()){
          throw std::runtime_error("partition size mismatch "+
                                    std::to_string(gids->size())+
                                    " vs "+
                                    std::to_string(field->size()));
        };
        for(int k=0;k<gids->size();k++){
          field->at(k)=gfield->at(gids->at(k));
        };
        return std::move(field);
      };

      /*!
      * Receives field partitions from the source process.
      * \param fmaker std::shared_ptr<PanNDE::FieldFactory> Factory for creating fields
      * \param field_source_rank int Rank of the source process
      * \return std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Field>>> Array of received field partitions
      */
      std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Field>>> recvFieldPartitions(
                            std::shared_ptr<PanNDE::FieldFactory> fmaker,
                            int field_source_rank){
        comm->sendArray(local_partition_ids,field_source_rank);
        PanNDE::Field::FieldType ftype;
        comm->recvValue((int*)&ftype,field_source_rank);
        auto field_array=fmaker->makeManagedFieldArray();
        field_array->resize(local_partition_ids->size());
        comm->waitall();
        for(int kp=0;kp<local_partition_ids->size();kp++){
          field_array->at(kp)=recvFieldPartition(kp,fmaker,ftype,field_source_rank);
        };
        return std::move(field_array);
      };

      /*!
      * Sends field partitions to a requesting process.
      * \param destination_rank int Rank of the requesting process
      * \param global_field std::shared_ptr<PanNDE::Field> The global field to partition and send
      */
      void sendFieldPartitions(int destination_rank,
                            std::shared_ptr<PanNDE::Field> global_field){
        auto part_ids=comm->recvArray(i32mfg,destination_rank);
        comm->sendValue(global_field->type(),destination_rank);
        comm->waitall();
        for(int kp=0;kp<part_ids->size();kp++){
          sendFieldPartition(part_ids->at(kp),global_field,destination_rank);
        };
      };

      /*!
      * Receives a single field partition from the source process.
      * \param idx int Index in the local partition array
      * \param fmaker std::shared_ptr<PanNDE::FieldFactory> Factory for creating fields
      * \param ftype PanNDE::Field::FieldType Type of the field being received
      * \param source_rank int Rank of the source process
      * \return std::shared_ptr<PanNDE::Field> The received field partition
      * \throws std::runtime_error If field size mismatches occur
      */
      std::shared_ptr<PanNDE::Field> recvFieldPartition(int idx,
                            std::shared_ptr<PanNDE::FieldFactory> fmaker,
                            PanNDE::Field::FieldType ftype,int source_rank){
        auto field=fmaker->makeManagedField(local_mesh_array->at(idx),ftype);
        //get data, copy data
        auto fdata=comm->recvArray(dblmfg,source_rank);comm->waitall();
        if(fdata->size()!=field->size()){
          throw std::runtime_error("field xfer size mismatch "+
                                    std::to_string(fdata->size())+
                                    " vs "+
                                    std::to_string(field->size()));
        };
        for(int k=0;k<fdata->size();k++){field->at(k)=fdata->at(k);};
        return std::move(field);
      };

      /*!
      * Sends a single field partition to a requesting process.
      * \param part_id int ID of the partition to send
      * \param gfield std::shared_ptr<PanNDE::Field> The global field
      * \param destination_rank int Rank of the requesting process
      */
      void sendFieldPartition(int part_id,std::shared_ptr<PanNDE::Field> gfield,
                              int destination_rank){
        //send data
        auto gids=(PanNDE::Field::NODE==gfield->type())?
                            getNodeGlobalIds(part_id):getCellGlobalIds(part_id);
        auto fdata=dblmfg->makeManagedArray();fdata->resize(gids->size());
        for(int k=0;k<gids->size();k++){fdata->at(k)=gfield->at(gids->at(k));};
        comm->sendArray(fdata,destination_rank);
        comm->waitall();
      };

      /*!
      * Creates mesh partitions for local use on the source process.
      * \param maker std::shared_ptr<PanNDE::MeshFactory> Factory for creating meshes
      * \return std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> Array of local mesh partitions
      */
      std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> makeLocalPartitions(
                                  std::shared_ptr<PanNDE::MeshFactory> maker){
        auto mesharray=maker->makeManagedMeshArray();
        mesharray->resize(local_partition_ids->size());
        for(int k=0;k<local_partition_ids->size();k++){
          mesharray->at(k)=makeLocalPartition(local_partition_ids->at(k),maker);
        };
        return std::move(mesharray);
      };

      /*!
      * Creates a single mesh partition for local use on the source process.
      * \param part_id int ID of the partition to create
      * \param maker std::shared_ptr<PanNDE::MeshFactory> Factory for creating meshes
      * \return std::shared_ptr<PanNDE::Mesh> The partitioned mesh
      */
      std::shared_ptr<PanNDE::Mesh> makeLocalPartition(int part_id,
                            std::shared_ptr<PanNDE::MeshFactory> maker){
        auto ngids=getNodeGlobalIds(part_id);
        auto coords=getNodeCoords(ngids);
        auto nowners=getNodeOwner(ngids);
        auto idx_map=buildPartitionNodesIdxMap(ngids);

        auto cgids=getCellGlobalIds(part_id);
        auto corners=getCellCorners(cgids,idx_map);
        auto cowners=getCellOwner(cgids,part_id);

        auto nodes=buildNodeVector(ngids,coords,nowners,maker);
        auto cells=buildCellVector(cgids,corners,cowners,maker);

        return std::move(maker->makeManagedMesh(nodes.data(),nodes.size(),
                                                cells.data(),cells.size(),part_id));
      };

      /*!
      * Gets the global node IDs for a given partition.
      * \param part_id int Partition ID
      * \return std::shared_ptr<PanNDE::Array<int64_t>> Array of global node IDs
      */
      std::shared_ptr<PanNDE::Array<int64_t>> getNodeGlobalIds(int part_id){
        return std::move(loadArrayFromSet(partitionEngine.getPartitionNodes(part_id)));
      };

      /*!
      * Gets the global cell IDs for a given partition.
      * \param part_id int Partition ID
      * \return std::shared_ptr<PanNDE::Array<int64_t>> Array of global cell IDs
      */
      std::shared_ptr<PanNDE::Array<int64_t>> getCellGlobalIds(int part_id){
        return std::move(loadArrayFromSet(partitionEngine.getPartitionCells(part_id)));
      };

      /*!
      * Converts a set of IDs to an array.
      * \param idxSet std::set<int64_t>* Set of IDs
      * \return std::shared_ptr<PanNDE::Array<int64_t>> Array containing the IDs
      */
      std::shared_ptr<PanNDE::Array<int64_t>> loadArrayFromSet(std::set<int64_t>* idxSet){
        auto array=i64mfg->makeManagedArray();
        array->resize(0);
        array->resize(idxSet->size());
        int cnt=0;
        for(auto it=idxSet->begin();it!=idxSet->end();it++){
          array->at(cnt)=(*it);cnt++;
        };
        return std::move(array);
      };

      /*!
      * Gets the coordinates of nodes specified by their global IDs.
      * \param node_ids std::shared_ptr<PanNDE::Array<int64_t>> Global node IDs
      * \return std::shared_ptr<PanNDE::Array<double>> Flattened array of 3D node coordinates
      */
      std::shared_ptr<PanNDE::Array<double>> getNodeCoords(
                                    std::shared_ptr<PanNDE::Array<int64_t>> node_ids){
        auto coords=dblmfg->makeManagedArray();
        coords->resize(0);coords->resize(3*(node_ids->size()));
        double pt[3];
        for(int k=0;k<node_ids->size();k++){
          gmesh->nodeCoordinate(node_ids->at(k),pt);
          for(int kd=0;kd<3;kd++){
            coords->at(3*k+kd)=pt[kd];
          };
        };
        return std::move(coords);
      };

      /*!
      * Gets the node indices for cells specified by their global IDs.
      * \param cell_ids std::shared_ptr<PanNDE::Array<int64_t>> Global cell IDs
      * \param node_idx_map std::map<int64_t,int32_t>& Mapping from global to local node indices
      * \return std::shared_ptr<PanNDE::Array<int32_t>> Flattened array of cell-node connectivity
      */
      std::shared_ptr<PanNDE::Array<int32_t>> getCellCorners(
                                    std::shared_ptr<PanNDE::Array<int64_t>> cell_ids,
                                    std::map<int64_t,int32_t>& node_idx_map){
        auto corners=i32mfg->makeManagedArray();
        corners->resize(0);corners->resize(8*(cell_ids->size()));
        int32_t box[8];
        for(int k=0;k<cell_ids->size();k++){
          gmesh->cell(cell_ids->at(k),box);
          for(int kb=0;kb<8;kb++){
            int local_nidx=node_idx_map.at(box[kb]);
            if(0>local_nidx){throw std::logic_error("nodal index map fault");};
            corners->at(8*k+kb)=node_idx_map.at(box[kb]);
          };
        };
        return std::move(corners);
      };

      /*!
      * Gets the owner partitions for nodes specified by their global IDs.
      * \param gids std::shared_ptr<PanNDE::Array<int64_t>> Global node IDs
      * \return std::shared_ptr<PanNDE::Array<int32_t>> Array of owner partition IDs
      */
      std::shared_ptr<PanNDE::Array<int32_t>> getNodeOwner(
                            std::shared_ptr<PanNDE::Array<int64_t>> gids){
        auto owners=i32mfg->makeManagedArray();
        auto nodeowners=partitionEngine.getFullDomainNodeOwnership();
        owners->resize(gids->size());
        for(auto k=0;k<gids->size();k++){owners->at(k)=nodeowners->at(gids->at(k));};
        return std::move(owners);
      };

      /*!
      * Gets the owner partitions for cells specified by their global IDs.
      * \param gids std::shared_ptr<PanNDE::Array<int64_t>> Global cell IDs
      * \param partition_id int Partition ID
      * \return std::shared_ptr<PanNDE::Array<int32_t>> Array of owner partition IDs
      */
      std::shared_ptr<PanNDE::Array<int32_t>> getCellOwner(
                            std::shared_ptr<PanNDE::Array<int64_t>> gids,
                            int partition_id){
        auto owners=i32mfg->makeManagedArray();
        auto cellowners=partitionEngine.getFullDomainCellOwnership();
        owners->resize(gids->size());
        for(auto k=0;k<gids->size();k++){owners->at(k)=cellowners->at(gids->at(k));};
        return std::move(owners);
      };

      /*!
      * Builds a mapping from global to local node indices for a partition.
      * \param gids std::shared_ptr<PanNDE::Array<int64_t>> Global node IDs
      * \return std::map<int64_t,int32_t> Map from global to local node indices
      */
      std::map<int64_t,int32_t> buildPartitionNodesIdxMap(
                        std::shared_ptr<PanNDE::Array<int64_t>> gids){
        std::map<int64_t,int32_t> idx_map;idx_map.clear();
        for(int32_t k=0;k<gids->size();k++){
          idx_map.emplace(gids->at(k),k);
        };
        return idx_map;
      };

      /*!
      * Builds a vector of nodes for a partition.
      * \param gids std::shared_ptr<PanNDE::Array<int64_t>> Global node IDs
      * \param coords std::shared_ptr<PanNDE::Array<double>> Node coordinates
      * \param nowners std::shared_ptr<PanNDE::Array<int32_t>> Node owner partitions
      * \param maker std::shared_ptr<PanNDE::MeshFactory> Factory for creating node objects
      * \return std::vector<PanNDE::Mesh::Node> Vector of node objects
      */
      std::vector<PanNDE::Mesh::Node> buildNodeVector(
                          std::shared_ptr<PanNDE::Array<int64_t>> gids,
                          std::shared_ptr<PanNDE::Array<double>> coords,
                          std::shared_ptr<PanNDE::Array<int32_t>> nowners,
                          std::shared_ptr<PanNDE::MeshFactory> maker){
        std::vector<PanNDE::Mesh::Node> nodes;nodes.resize(gids->size());
        double pt[3];
        for(int kn=0;kn<nodes.size();kn++){
          for(int kd=0;kd<3;kd++){pt[kd]=coords->at(3*kn+kd);};
          nodes.at(kn)=maker->makeNode(pt,gids->at(kn),nowners->at(kn));
        };
        return nodes;
      };

      /*!
      * Builds a vector of cells for a partition.
      * \param gids std::shared_ptr<PanNDE::Array<int64_t>> Global cell IDs
      * \param boxes std::shared_ptr<PanNDE::Array<int32_t>> Cell-node connectivity
      * \param cowners std::shared_ptr<PanNDE::Array<int32_t>> Cell owner partitions
      * \param maker std::shared_ptr<PanNDE::MeshFactory> Factory for creating cell objects
      * \return std::vector<PanNDE::Mesh::Cell> Vector of cell objects
      */
      std::vector<PanNDE::Mesh::Cell> buildCellVector(
                          std::shared_ptr<PanNDE::Array<int64_t>> gids,
                          std::shared_ptr<PanNDE::Array<int32_t>> boxes,
                          std::shared_ptr<PanNDE::Array<int32_t>> cowners,
                          std::shared_ptr<PanNDE::MeshFactory> maker){
        std::vector<PanNDE::Mesh::Cell> cells;cells.resize(gids->size());
        int32_t box[8];
        for(int kc=0;kc<cells.size();kc++){
          for(int kb=0;kb<8;kb++){box[kb]=boxes->at(8*kc+kb);};
          cells.at(kc)=maker->makeCell(box,8,gids->at(kc),cowners->at(kc));
        };
        return cells;
      };

      /*!
      * Receives mesh partitions from the source process.
      * \param maker std::shared_ptr<PanNDE::MeshFactory> Factory for creating meshes
      * \param source_rank int Rank of the source process
      * \return std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> Array of received mesh partitions
      */
      std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> recvPartitions(
                            std::shared_ptr<PanNDE::MeshFactory> maker,
                            int source_rank){
        comm->sendArray(local_partition_ids,source_rank);
        auto mesharray=maker->makeManagedMeshArray();
        mesharray->resize(local_partition_ids->size());
        comm->waitall();
        for(int k=0;k<local_partition_ids->size();k++){
          mesharray->at(k)=recvPartition(local_partition_ids->at(k),maker,source_rank);
        };
        return std::move(mesharray);
      };

      /*!
      * Sends mesh partitions to a requesting process.
      * \param destination_rank int Rank of the requesting process
      */
      void sendPartitions(int destination_rank){
        auto remote_part_ids=comm->recvArray(i32mfg,destination_rank);
        comm->waitall();
        for(int kp=0;kp<remote_part_ids->size();kp++){
          sendPartition(remote_part_ids->at(kp),destination_rank);
        };
      };

      /*!
      * Receives a single mesh partition from the source process.
      * \param part_id int ID of the partition to receive
      * \param maker std::shared_ptr<PanNDE::MeshFactory> Factory for creating meshes
      * \param source_rank int Rank of the source process
      * \return std::shared_ptr<PanNDE::Mesh> The received mesh partition
      */
      std::shared_ptr<PanNDE::Mesh> recvPartition(int part_id,
                            std::shared_ptr<PanNDE::MeshFactory> maker,
                            int source_rank){
        auto local_nodes=recvPartitionNodes(part_id,maker,source_rank);
        auto local_cells=recvPartitionCells(part_id,maker,source_rank);
        printf("[%i] data received.\n",comm->getProcessId());
        return std::move(maker->makeManagedMesh(local_nodes.data(),local_nodes.size(),
                                                local_cells.data(),local_cells.size(),part_id));
      };

      /*!
      * Sends a single mesh partition to a requesting process.
      * \param part_id int ID of the partition to send
      * \param destination_rank int Rank of the requesting process
      * \return std::map<int64_t,int32_t> Map from global to local node indices
      */
      void sendPartition(int part_id,int destination_rank){
        auto idx_map=sendPartitionNodes(part_id,destination_rank);
        sendPartitionCells(part_id,destination_rank,idx_map);
      };

      /*!
      * Receives node data for a partition from the source process.
      * \param part_id int ID of the partition being received
      * \param maker std::shared_ptr<PanNDE::MeshFactory> Factory for creating node objects
      * \param source_rank int Rank of the source process
      * \return std::vector<PanNDE::Mesh::Node> Vector of received nodes
      */
      std::vector<PanNDE::Mesh::Node> recvPartitionNodes(int part_id,
                          std::shared_ptr<PanNDE::MeshFactory> maker,
                          int source_rank){
        auto gids=comm->recvArray(i64mfg,source_rank);
        comm->waitall();
        auto coords=comm->recvArray(dblmfg,source_rank);
        auto nowners=comm->recvArray(i32mfg,source_rank);
        comm->waitall();
        return buildNodeVector(gids,coords,nowners,maker);
      };

      /*!
      * Sends node data for a partition to a requesting process.
      * \param part_id int ID of the partition to send
      * \param destination_rank int Rank of the requesting process
      * \return std::map<int64_t,int32_t> Map from global to local node indices
      */
      std::map<int64_t,int32_t> sendPartitionNodes(int part_id,int destination_rank){
        auto gids=getNodeGlobalIds(part_id);
        comm->sendArray(gids,destination_rank);
        comm->waitall();
        auto node_coords_to_send=getNodeCoords(gids);
        comm->sendArray(node_coords_to_send,destination_rank);
        auto node_owners_to_send=getNodeOwner(gids);
        comm->sendArray(node_owners_to_send,destination_rank);
        comm->waitall();
        return buildPartitionNodesIdxMap(gids);
      };

      /*!
      * Receives cell data for a partition from the source process.
      * \param part_id int ID of the partition being received
      * \param maker std::shared_ptr<PanNDE::MeshFactory> Factory for creating cell objects
      * \param source_rank int Rank of the source process
      * \return std::vector<PanNDE::Mesh::Cell> Vector of received cells
      */
      std::vector<PanNDE::Mesh::Cell> recvPartitionCells(int part_id,
                          std::shared_ptr<PanNDE::MeshFactory> maker,
                          int source_rank){
        auto gids=comm->recvArray(i64mfg,source_rank);
        comm->waitall();
        auto boxes=comm->recvArray(i32mfg,source_rank);
        auto cowners=comm->recvArray(i32mfg,source_rank);
        comm->waitall();
        return buildCellVector(gids,boxes,cowners,maker);
      };

      /*!
      * Sends cell data for a partition to a requesting process.
      * \param part_id int ID of the partition to send
      * \param destination_rank int Rank of the requesting process
      * \param idx_map std::map<int64_t,int32_t> Map from global to local node indices
      */
      void sendPartitionCells(int part_id,int destination_rank,
                              std::map<int64_t,int32_t> idx_map){
        auto gids=getCellGlobalIds(part_id);
        comm->sendArray(gids,destination_rank);
        comm->waitall();
        auto corners=getCellCorners(gids,idx_map);
        comm->sendArray(corners,destination_rank);
        auto cell_owners=getCellOwner(gids,part_id);
        comm->sendArray(cell_owners,destination_rank);
        comm->waitall();
      };

      /*!
      * Configures the communicator for the partitioner.
      * \param communicator std::shared_ptr<PanNDE::Communicator> Communicator to use (or create a new one if nullptr)
      */
      void configure(std::shared_ptr<PanNDE::Communicator> communicator){
        if(nullptr==communicator){
          comm=std::make_shared<NetMPI::MPICommunicator>(NetMPI::MPICommunicator());
        }else{comm=communicator;};
        rank=comm->getProcessId();
        nranks=comm->getNumberOfProcesses();
      };

      /*!
      * Throws an error for an invalid mesh.
      * \throws std::runtime_error With message about invalid mesh
      */
      void invalidMeshError(){throw std::runtime_error("invalid mesh");};

      /*!
      * Runs the partitioning algorithm on the global mesh.
      */
      void goPartition(){
        printf("[%i] Executing partitioner\n",rank);
        partitionEngine.partition(gmesh,partition_count);
        printf("[%i] partitioner executed\n",rank);
      };

      //! Process rank
      int rank=0;
      //! Total number of processes
      int nranks=1;
      //! Number of partitions to create
      int partition_count=1;
      //! Rank of the process performing partitioning
      int partitioner_rank=0;
      //! The global mesh to partition
      std::shared_ptr<PanNDE::Mesh> gmesh=nullptr;
      //! Communicator for inter-process data exchange
      std::shared_ptr<PanNDE::Communicator> comm=nullptr;

      //! Internal partitioner implementation
      NetMPI::MetisHexPartitioner partitionEngine;

      //! Factory for 32-bit integer arrays
      std::shared_ptr<PanNDE::ArrayFactory<int32_t>> i32mfg;
      //! Factory for 64-bit integer arrays
      std::shared_ptr<PanNDE::ArrayFactory<int64_t>> i64mfg;
      //! Factory for double-precision arrays
      std::shared_ptr<PanNDE::ArrayFactory<double>> dblmfg;
      
      //! Array of partition IDs assigned to this process
      std::shared_ptr<PanNDE::Array<int32_t>> local_partition_ids=nullptr;
      //! Array of local mesh partitions
      std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> local_mesh_array=nullptr;
  };
};