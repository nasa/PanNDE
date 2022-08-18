/*! \headerfile MetisPartitioner.hpp "modules/NetMPI/include/MetisPartitioner.hpp"
* "MetisPartitioner.hpp" contains the class implementation for partitioning and distributing a mesh 
* and it's associated fields using metis
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
  * Implements the required externally faced operations for partitioning and distributing a mesh 
  * and it's associated fields
  *
  */
  class MetisPartitioner : public Controller::Partitioner {
    public:
      /*!
      * constructor
      * \param i32_maker Array factory for 32bit ints
      * \param i64_maker Array factory for 64bit ints
      * \param dbl_maker Array factory for doubles
      * \param communicator communicator tool for distributing mesh and fields
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
      * partition a mesh
      * \param Nparts number of parts to decompose the mesh into
      * \param global_mesh the monolithic mesh to partition
      */
      void partitionMesh(int Nparts,std::shared_ptr<PanNDE::Mesh> global_mesh)override{
        if(nullptr==global_mesh){invalidMeshError();};
        gmesh=global_mesh;
        partition_count=Nparts;
        goPartition();
      };
      /*!
      * distribute the mesh
      * \param partitionIds which partitions the local process is requesting
      * \param maker the mesh factory
      * \param mesh_source_rank process which has the global mesh and is distributing the mesh partitions
      */
      std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> distributeMeshPartitions(
                                  std::shared_ptr<PanNDE::Array<int>> partitionIds,
                                  std::shared_ptr<PanNDE::MeshFactory> maker,
                                  int mesh_source_rank=0)override{
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
      * distribute a field
      * \param partitionIds which partitions the local process is requesting
      * \param maker the field factory
      * \param global_field the field that must be partitioned that is currently on the field_source_rank process
      * \param field_source_rank process which has the global field and is distributing the field partitions
      */
      std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Field>>> distributeFieldPartitions(
                                  std::shared_ptr<PanNDE::Array<int>> partitionIds,
                                  std::shared_ptr<PanNDE::FieldFactory> maker,
                                  std::shared_ptr<PanNDE::Field> global_field=nullptr,
                                  int field_source_rank=0)override{
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
      void copyIds(std::shared_ptr<PanNDE::Array<int32_t>> part_ids){
        local_partition_ids=i32mfg->makeManagedArray();
        local_partition_ids->resize(part_ids->size());
        for(int k=0;k<part_ids->size();k++){local_partition_ids->at(k)=part_ids->at(k);};
      };

      void checkFieldDistributionReadiness(std::shared_ptr<PanNDE::Field> global_field,int src_rank){
        if(nullptr==local_mesh_array){throw std::runtime_error("mesh must be distributed");};
        if(rank==src_rank){
          if(nullptr==global_field){throw std::runtime_error("no mesh provided");};
          if(PanNDE::Field::NODE!=global_field->type() && PanNDE::Field::CELL!=global_field->type()){
            throw std::runtime_error("invalid field type for redistribution");
          };
        };
      };
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
      void sendFieldPartitions(int destination_rank,
                            std::shared_ptr<PanNDE::Field> global_field){
        auto part_ids=comm->recvArray(i32mfg,destination_rank);
        comm->sendValue(global_field->type(),destination_rank);
        comm->waitall();
        for(int kp=0;kp<part_ids->size();kp++){
          sendFieldPartition(part_ids->at(kp),global_field,destination_rank);
        };
      };
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

      std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> makeLocalPartitions(
                                  std::shared_ptr<PanNDE::MeshFactory> maker){
        auto mesharray=maker->makeManagedMeshArray();
        mesharray->resize(local_partition_ids->size());
        for(int k=0;k<local_partition_ids->size();k++){
          mesharray->at(k)=makeLocalPartition(local_partition_ids->at(k),maker);
        };
        return std::move(mesharray);
      };
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
      std::shared_ptr<PanNDE::Array<int64_t>> getNodeGlobalIds(int part_id){
        return std::move(loadArrayFromSet(partitionEngine.getPartitionNodes(part_id)));
      };
      std::shared_ptr<PanNDE::Array<int64_t>> getCellGlobalIds(int part_id){
        return std::move(loadArrayFromSet(partitionEngine.getPartitionCells(part_id)));
      };
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
      std::shared_ptr<PanNDE::Array<int32_t>> getNodeOwner(
                            std::shared_ptr<PanNDE::Array<int64_t>> gids){
        auto owners=i32mfg->makeManagedArray();
        auto nodeowners=partitionEngine.getFullDomainNodeOwnership();
        owners->resize(gids->size());
        for(auto k=0;k<gids->size();k++){owners->at(k)=nodeowners->at(gids->at(k));};
        return std::move(owners);
      };
      std::shared_ptr<PanNDE::Array<int32_t>> getCellOwner(
                            std::shared_ptr<PanNDE::Array<int64_t>> gids,
                            int partition_id){
        auto owners=i32mfg->makeManagedArray();
        auto cellowners=partitionEngine.getFullDomainCellOwnership();
        owners->resize(gids->size());
        for(auto k=0;k<gids->size();k++){owners->at(k)=cellowners->at(gids->at(k));};
        return std::move(owners);
      };
      std::map<int64_t,int32_t> buildPartitionNodesIdxMap(
                        std::shared_ptr<PanNDE::Array<int64_t>> gids){
        std::map<int64_t,int32_t> idx_map;idx_map.clear();
        for(int32_t k=0;k<gids->size();k++){
          idx_map.emplace(gids->at(k),k);
        };
        return idx_map;
      };
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
      void sendPartitions(int destination_rank){
        auto remote_part_ids=comm->recvArray(i32mfg,destination_rank);
        comm->waitall();
        for(int kp=0;kp<remote_part_ids->size();kp++){
          sendPartition(remote_part_ids->at(kp),destination_rank);
        };
      };


      std::shared_ptr<PanNDE::Mesh> recvPartition(int part_id,
                            std::shared_ptr<PanNDE::MeshFactory> maker,
                            int source_rank){
        auto local_nodes=recvPartitionNodes(part_id,maker,source_rank);
        auto local_cells=recvPartitionCells(part_id,maker,source_rank);
        printf("[%i] data received.\n",comm->getProcessId());
        return std::move(maker->makeManagedMesh(local_nodes.data(),local_nodes.size(),
                                                local_cells.data(),local_cells.size(),part_id));
      };
      void sendPartition(int part_id,int destination_rank){
        auto idx_map=sendPartitionNodes(part_id,destination_rank);
        sendPartitionCells(part_id,destination_rank,idx_map);
      };
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


      void configure(std::shared_ptr<PanNDE::Communicator> communicator){
        if(nullptr==communicator){
          comm=std::make_shared<NetMPI::MPICommunicator>(NetMPI::MPICommunicator());
        }else{comm=communicator;};
        rank=comm->getProcessId();
        nranks=comm->getNumberOfProcesses();
      };
      void invalidMeshError(){throw std::runtime_error("invalid mesh");};
      void goPartition(){
        printf("[%i] Executing partitioner\n",rank);
        partitionEngine.partition(gmesh,partition_count);
        printf("[%i] partitioner executed\n",rank);
      };

      int rank=0;int nranks=1;
      int partition_count=1;
      int partitioner_rank=0;
      std::shared_ptr<PanNDE::Mesh> gmesh=nullptr;
      std::shared_ptr<PanNDE::Communicator> comm=nullptr;

      NetMPI::MetisHexPartitioner partitionEngine;

      std::shared_ptr<PanNDE::ArrayFactory<int32_t>> i32mfg;
      std::shared_ptr<PanNDE::ArrayFactory<int64_t>> i64mfg;
      std::shared_ptr<PanNDE::ArrayFactory<double>> dblmfg;
      
      std::shared_ptr<PanNDE::Array<int32_t>> local_partition_ids=nullptr;
      std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> local_mesh_array=nullptr;
  };
};

