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

#include <memory>
#include <vector>
#include <cmath>

#include "Partitioner.hpp"

#include "gtest/gtest.h"

using testing::Types;

template<class ST>
std::shared_ptr<ST> makeShared(){
  return std::move(std::make_shared<ST>(ST()));
};

template <class I,class C,class MS,class MF,class FF,template<typename T>class AF>
struct PartitionerTypeDefinitions
{
  typedef I Impl;
  typedef C Comm;
  typedef MS MeshStub;
  typedef MF MeshFactory;
  typedef FF FieldFactory;
  typedef AF<int> I32ArrayFactory;
  typedef AF<int64_t> I64ArrayFactory;
  typedef AF<double> DblArrayFactory;
};

template<class T>
class PartitionerTests : public ::testing::Test {
  protected:
    void SetUp() override {
      rank=communicator->getProcessId();
      comm_size=communicator->getNumberOfProcesses();
    };
    void TearDown() override {};

    std::shared_ptr<typename T::Impl> partitioner=nullptr;
    std::shared_ptr<typename T::Comm> communicator=makeShared<typename T::Comm>();
    std::shared_ptr<typename T::MeshStub> global_mesh=nullptr;
    std::shared_ptr<PanNDE::Field> gfield=nullptr;
    std::shared_ptr<typename T::MeshFactory> mesh_maker=makeShared<typename T::MeshFactory>();
    std::shared_ptr<typename T::FieldFactory> field_maker=makeShared<typename T::FieldFactory>();
    std::shared_ptr<typename T::I32ArrayFactory> i32_maker=makeShared<typename T::I32ArrayFactory>();
    std::shared_ptr<typename T::I64ArrayFactory> i64_maker=makeShared<typename T::I64ArrayFactory>();
    std::shared_ptr<typename T::DblArrayFactory> dbl_maker=makeShared<typename T::DblArrayFactory>();
    int rank;int comm_size;

    void makePartitioner(){
      partitioner=std::make_shared<typename T::Impl>(i32_maker,i64_maker,
                                                     dbl_maker,communicator);
      EXPECT_NE(nullptr,partitioner);
    };
    void makeCommunicator(){
      EXPECT_NE(nullptr,communicator);
    };
    void makeFactories(){
      EXPECT_NE(nullptr,mesh_maker);
      EXPECT_NE(nullptr,field_maker);
      EXPECT_NE(nullptr,i32_maker);
    };
    void makeGlobalMesh(){
      if(0==rank){
        global_mesh=std::make_shared<typename T::MeshStub>(0,1);
        EXPECT_NE(nullptr,global_mesh);
      };
      if(0!=rank){EXPECT_EQ(nullptr,global_mesh);};
    };
    void makeInfrastructure(){
      makeCommunicator();
      makePartitioner();
      makeFactories();
      makeGlobalMesh();
    };
    void partition(int Nparts){
      if(0==rank){partitioner->partitionMesh(Nparts,global_mesh);};
    };
    std::shared_ptr<PanNDE::Array<int32_t>> makePartIds(int parts_per_rank=1){
      auto partids=i32_maker->makeManagedArray();
      partids->resize(parts_per_rank);
      for(int k=0;k<parts_per_rank;k++){
        partids->at(k)=parts_per_rank*rank+k;
      };
      return std::move(partids);
    };
    void makeGlobalNodeField(){
      if(0==rank){
        gfield=field_maker->makeManagedField(global_mesh,PanNDE::Field::NODE);
        for(int k=0;k<gfield->size();k++){gfield->at(k)=0.1*(double(k));};
      };
    };
};

TYPED_TEST_SUITE_P(PartitionerTests);

TYPED_TEST_P(PartitionerTests,Nothing){};

TYPED_TEST_P(PartitionerTests,MakePartitioner){
  this->makePartitioner();
};
TYPED_TEST_P(PartitionerTests,MakeCommunicator){
  this->makeCommunicator();
};
TYPED_TEST_P(PartitionerTests,MakeFactories){
  this->makeFactories();
};
TYPED_TEST_P(PartitionerTests,MakeGlobalMesh){
  this->makeGlobalMesh();
};
TYPED_TEST_P(PartitionerTests,MakeMeshPartition){
  this->makeInfrastructure();
  this->partition(this->comm_size);
  auto partids=this->makePartIds(1);
  auto mesharray=this->partitioner->distributeMeshPartitions(partids,this->mesh_maker,0);
  
  ASSERT_EQ(1,mesharray->size());
  EXPECT_NE(nullptr,mesharray->at(0));
  EXPECT_LT(0,mesharray->at(0)->nodeCount());
  EXPECT_LT(0,mesharray->at(0)->cellCount());
  int32_t box[8];
  for(int k=0;k<mesharray->at(0)->cellCount();k++){
    mesharray->at(0)->cell(k,box);
    for(int kb=0;kb<8;kb++){EXPECT_NE(-1,box[kb]);};
  };
  double pt[3];
  for(int k=0;k<mesharray->at(0)->nodeCount();k++){
    mesharray->at(0)->nodeCoordinate(k,pt);
    for(int kb=0;kb<3;kb++){EXPECT_LE(0.0,pt[kb]);};
  };
};
TYPED_TEST_P(PartitionerTests,MakeMeshPartition2){
  this->makeInfrastructure();
  int Nparts=2*(this->comm_size);
  this->partition(Nparts);
  auto partids=this->makePartIds(2);
  auto mesharray=this->partitioner->distributeMeshPartitions(partids,this->mesh_maker,0);
  
  ASSERT_EQ(2,mesharray->size());
  for(int km=0;km<2;km++){
    EXPECT_NE(nullptr,mesharray->at(km));
    EXPECT_LT(0,mesharray->at(km)->nodeCount());
    EXPECT_LT(0,mesharray->at(km)->cellCount());
    int32_t box[8];
    for(int k=0;k<mesharray->at(km)->cellCount();k++){
      mesharray->at(km)->cell(k,box);
      for(int kb=0;kb<8;kb++){EXPECT_NE(-1,box[kb]);};
    };
    double pt[3];
    for(int k=0;k<mesharray->at(km)->nodeCount();k++){
      mesharray->at(km)->nodeCoordinate(k,pt);
      for(int kb=0;kb<3;kb++){EXPECT_LE(0.0,pt[kb]);};
    };
  };
};
TYPED_TEST_P(PartitionerTests,MakeMeshPartitionFields){
  this->makeInfrastructure();
  this->partition(this->comm_size);
  auto partids=this->makePartIds(1);
  auto mesharray=this->partitioner->distributeMeshPartitions(partids,this->mesh_maker,0);
  ASSERT_EQ(1,mesharray->size());
  EXPECT_NE(nullptr,mesharray->at(0));
  this->makeGlobalNodeField();
  auto fieldarray=this->partitioner->distributeFieldPartitions(partids,
                                          this->field_maker,this->gfield,0);

  ASSERT_NE(nullptr,fieldarray);
  ASSERT_EQ(1,fieldarray->size());
  ASSERT_NE(nullptr,fieldarray->at(0));
  for(int k=0;k<fieldarray->at(0)->size();k++){
    EXPECT_EQ(0.1*(double(mesharray->at(0)->globalNodeId(k))),fieldarray->at(0)->at(k));
  };
};
TYPED_TEST_P(PartitionerTests,MakeMeshPartitionFields2){
  this->makeInfrastructure();
  this->partition(2*(this->comm_size));
  auto partids=this->makePartIds(2);
  auto mesharray=this->partitioner->distributeMeshPartitions(partids,this->mesh_maker,0);
  ASSERT_EQ(2,mesharray->size());
  EXPECT_NE(nullptr,mesharray->at(0));
  EXPECT_NE(nullptr,mesharray->at(1));
  this->makeGlobalNodeField();
  auto fieldarray=this->partitioner->distributeFieldPartitions(partids,
                                          this->field_maker,this->gfield,0);

  ASSERT_NE(nullptr,fieldarray);
  ASSERT_EQ(2,fieldarray->size());
  ASSERT_NE(nullptr,fieldarray->at(0));
  ASSERT_NE(nullptr,fieldarray->at(1));
  for(int k=0;k<fieldarray->at(0)->size();k++){
    EXPECT_EQ(0.1*(double(mesharray->at(0)->globalNodeId(k))),fieldarray->at(0)->at(k));
  };
  for(int k=0;k<fieldarray->at(1)->size();k++){
    EXPECT_EQ(0.1*(double(mesharray->at(1)->globalNodeId(k))),fieldarray->at(1)->at(k));
  };
};

REGISTER_TYPED_TEST_SUITE_P(PartitionerTests,
                            Nothing,MakePartitioner,MakeCommunicator,MakeFactories,
                            MakeGlobalMesh,MakeMeshPartition,MakeMeshPartition2,
                            MakeMeshPartitionFields,MakeMeshPartitionFields2);
