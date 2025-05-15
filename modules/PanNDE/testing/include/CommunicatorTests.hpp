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
#include <array>
#include <cmath>

#include "Array.hpp"

#include "gtest/gtest.h"

using testing::Types;

template<class FT>
std::shared_ptr<FT> makeShared(){
  return std::move(std::make_shared<FT>(FT()));
};

template <class I,class MF,class FF,template<typename T>class AF>
struct CommunicatorTypeDefinitions
{
  typedef I Impl;
  typedef MF mesh;
  typedef FF field_factory;
  typedef AF<int64_t> i64array_factory;
  typedef AF<int32_t> i32array_factory;
  typedef AF<double> darray_factory;
  typedef AF<PanNDE::Mesh::Node> nodeArray_factory;
  typedef AF<PanNDE::Mesh::Cell> cellArray_factory;
};

template<class T>
class CommunicatorTests : public ::testing::Test {
  protected:
    void SetUp() override {};
    void TearDown() override {};
    void getRankAndSize(){
      rank=communicator->getProcessId();
      nranks=communicator->getNumberOfProcesses();
      ASSERT_LE(1,nranks);
      ASSERT_LE(0,rank);
      ASSERT_GT(nranks,rank);
    };
    int rank,nranks;
    std::shared_ptr<PanNDE::Communicator> communicator=makeShared<typename T::Impl>();
    std::shared_ptr<PanNDE::ArrayFactory<int64_t>> i64array_maker=makeShared<typename T::i64array_factory>();
    std::shared_ptr<PanNDE::ArrayFactory<int32_t>> i32array_maker=makeShared<typename T::i32array_factory>();
    std::shared_ptr<PanNDE::ArrayFactory<double>> darray_maker=makeShared<typename T::darray_factory>();
};

TYPED_TEST_SUITE_P(CommunicatorTests);

TYPED_TEST_P(CommunicatorTests,Nothing){};

TYPED_TEST_P(CommunicatorTests,GetMyRankAndSize){
  this->getRankAndSize();
};

TYPED_TEST_P(CommunicatorTests,BroadCastRanks){
  this->getRankAndSize();
  std::vector<int> ranks;ranks.resize(this->nranks);
  ranks.at(this->rank)=this->rank;
  for(int k=0;k<this->nranks;k++){
    this->communicator->broadcastValue(&ranks.at(k),k);
  };
  for(int k=0;k<this->nranks;k++){
    EXPECT_EQ(k,ranks.at(k));
  };
};

TYPED_TEST_P(CommunicatorTests,BroadCastRanksAsDoubles){
  this->getRankAndSize();
  std::vector<double> ranks;ranks.resize(this->nranks);
  ranks.at(this->rank)=0.1*(double(this->rank));
  for(int k=0;k<this->nranks;k++){
    this->communicator->broadcastValue(&ranks.at(k),k);
  };
  for(int k=0;k<this->nranks;k++){
    EXPECT_EQ(0.1*(double(k)),ranks.at(k));
  };
};

TYPED_TEST_P(CommunicatorTests,BroadCastRanksAsI64){
  this->getRankAndSize();
  std::vector<int64_t> ranks;ranks.resize(this->nranks);
  ranks.at(this->rank)=10*(this->rank);
  for(int k=0;k<this->nranks;k++){
    this->communicator->broadcastValue(&ranks.at(k),k);
  };
  for(int k=0;k<this->nranks;k++){
    EXPECT_EQ(10*k,ranks.at(k));
  };
};

TYPED_TEST_P(CommunicatorTests,DoubleRoundRobin){
  double x;
  this->getRankAndSize();
  if(0==this->rank){x=234.235;};
  for(int k=1;k<this->nranks;k++){
    if(k-1==this->rank){this->communicator->sendValue(x,k);};
    if(k==this->rank){this->communicator->recvValue(&x,k-1);};
    this->communicator->waitall();
  };
  EXPECT_EQ(234.235,x);
};

TYPED_TEST_P(CommunicatorTests,IntRoundRobin){
  int x;
  this->getRankAndSize();
  if(0==this->rank){x=234;};
  for(int k=1;k<this->nranks;k++){
    if(k-1==this->rank){this->communicator->sendValue(x,k);};
    if(k==this->rank){this->communicator->recvValue(&x,k-1);};
    this->communicator->waitall();
  };
  EXPECT_EQ(234,x);
};

TYPED_TEST_P(CommunicatorTests,Int64RoundRobin){
  int64_t x;
  this->getRankAndSize();
  if(0==this->rank){x=234;};
  for(int k=1;k<this->nranks;k++){
    if(k-1==this->rank){this->communicator->sendValue(x,k);};
    if(k==this->rank){this->communicator->recvValue(&x,k-1);};
    this->communicator->waitall();
  };
  EXPECT_EQ(234,x);
};

TYPED_TEST_P(CommunicatorTests,BroadCastArrayI32){
  this->getRankAndSize();
  int N=10;
  auto array=this->i32array_maker->makeManagedArray();
  array->resize(N);
  for(int kp=0;kp<this->communicator->getNumberOfProcesses();kp++){
    if(kp==this->communicator->getProcessId()){
      for(int k=0;k<N;k++){array->at(k)=N*kp+k+1;};
    };
    this->communicator->broadcastArray(array,kp);
    for(int k=0;k<N;k++){EXPECT_EQ(N*kp+k+1,array->at(k));};
  };
};

TYPED_TEST_P(CommunicatorTests,BroadCastArrayI64){
  this->getRankAndSize();
  int N=10;
  auto array=this->i64array_maker->makeManagedArray();
  array->resize(N);
  for(int kp=0;kp<this->communicator->getNumberOfProcesses();kp++){
    if(kp==this->communicator->getProcessId()){
      for(int k=0;k<N;k++){array->at(k)=int64_t(N*kp+k+1);};
    };
    this->communicator->broadcastArray(array,kp);
    for(int k=0;k<N;k++){EXPECT_EQ(N*kp+k+1,array->at(k));};
  };
};

TYPED_TEST_P(CommunicatorTests,BroadCastArrayF64){
  this->getRankAndSize();
  int N=10;
  auto array=this->darray_maker->makeManagedArray();
  array->resize(N);
  for(int kp=0;kp<this->communicator->getNumberOfProcesses();kp++){
    if(kp==this->communicator->getProcessId()){
      for(int k=0;k<N;k++){array->at(k)=(double(N*kp+k+1))*0.1;};
    };
    this->communicator->broadcastArray(array,kp);
    for(int k=0;k<N;k++){EXPECT_EQ(double(N*kp+k+1)*0.1,array->at(k));};
  };
};

TYPED_TEST_P(CommunicatorTests,AllGatherI32){
  this->getRankAndSize();
  int value=this->communicator->getProcessId()+1;
  auto array=this->communicator->allGatherValue(value,this->i32array_maker);
  for(int k=0;k<this->communicator->getNumberOfProcesses();k++){
    EXPECT_EQ(k+1,array->at(k));
  };
};
TYPED_TEST_P(CommunicatorTests,AllGatherI64){
  this->getRankAndSize();
  int64_t value=10*(this->communicator->getProcessId()+1);
  auto array=this->communicator->allGatherValue(value,this->i64array_maker);
  for(int k=0;k<this->communicator->getNumberOfProcesses();k++){
    EXPECT_EQ(10*(k+1),array->at(k));
  };
};
TYPED_TEST_P(CommunicatorTests,AllGatherF64){
  this->getRankAndSize();
  double value=0.1*(double(this->rank+1));
  auto array=this->communicator->allGatherValue(value,this->darray_maker);
  for(int k=0;k<this->nranks;k++){
    EXPECT_EQ(0.1*(double(k+1)),array->at(k));
  };
};

TYPED_TEST_P(CommunicatorTests,DoubleArrayRoundRobin){
  this->getRankAndSize();
  int N=10;
  auto array=this->darray_maker->makeManagedArray();
  if(0==this->rank){
    array->resize(N);
    for(int k=0;k<N;k++){array->at(k)=(double(k+1))*0.1;};};
  for(int k=1;k<this->nranks;k++){
    if(k-1==this->rank){this->communicator->sendArray(array,k);};
    if(k==this->rank){array=this->communicator->recvArray(this->darray_maker,k-1);};
    this->communicator->waitall();
  };
  for(int k=0;k<N;k++){EXPECT_EQ(double(k+1)*0.1,array->at(k));};
};

TYPED_TEST_P(CommunicatorTests,I32ArrayRoundRobin){
  this->getRankAndSize();
  int N=10;
  auto array=this->i32array_maker->makeManagedArray();
  if(0==this->rank){
    array->resize(N);
    for(int k=0;k<N;k++){array->at(k)=(k+1);};};
  for(int k=1;k<this->nranks;k++){
    if(k-1==this->rank){this->communicator->sendArray(array,k);};
    if(k==this->rank){array=this->communicator->recvArray(this->i32array_maker,k-1);};
    this->communicator->waitall();
  };
  for(int k=0;k<N;k++){EXPECT_EQ((k+1),array->at(k));};
};

TYPED_TEST_P(CommunicatorTests,I64ArrayRoundRobin){
  this->getRankAndSize();
  int N=10;
  auto array=this->i64array_maker->makeManagedArray();
  if(0==this->rank){
    array->resize(N);
    for(int k=0;k<N;k++){array->at(k)=(k+1);};};
  for(int k=1;k<this->nranks;k++){
    if(k-1==this->rank){this->communicator->sendArray(array,k);};
    if(k==this->rank){array=this->communicator->recvArray(this->i64array_maker,k-1);};
    this->communicator->waitall();
  };
  for(int k=0;k<N;k++){EXPECT_EQ((k+1),array->at(k));};
};

TYPED_TEST_P(CommunicatorTests,HaloExchange){
  this->getRankAndSize();
  std::shared_ptr<PanNDE::Mesh> pmesh=
      std::make_shared<typename TypeParam::mesh>(
          typename TypeParam::mesh(this->rank,this->nranks));
  ASSERT_NE(nullptr,pmesh);
  auto field_maker=makeShared<typename TypeParam::field_factory>();
  ASSERT_NE(nullptr,field_maker);
  auto pbundle=field_maker->makeEmptyManagedFieldBundle();
  pbundle->mesh()=pmesh;
  pbundle->emplaceField("node_field",PanNDE::Field::NODE);
  pbundle->emplaceField("cell_field",PanNDE::Field::CELL);
  auto node_field=pbundle->field("node_field");
  auto cell_field=pbundle->field("cell_field");
  ASSERT_NE(nullptr,pbundle);
  ASSERT_NE(nullptr,pbundle->field("node_field"));
  ASSERT_NE(nullptr,pbundle->field("cell_field"));
  EXPECT_EQ(PanNDE::Field::NODE,pbundle->field("node_field")->type());
  EXPECT_EQ(PanNDE::Field::CELL,pbundle->field("cell_field")->type());

  for(int k=0;k<cell_field->size();k++){
    cell_field->at(k)=(double(-0.1*k));
  };
  for(int k=0;k<node_field->size();k++){
    node_field->at(k)=(double(0.1*k));
  };
  this->communicator->setupDataLinks(pbundle);
  for(int k=0;k<pbundle->fieldCount();k++){
    this->communicator->startHaloExchange(pbundle->fieldName(k));
  };
  for(int k=0;k<pbundle->fieldCount();k++){
    this->communicator->waitUntilDone(pbundle->fieldName(k));
  };
  for(int kn=0;kn<pmesh->nodeCount();kn++){
    if(pmesh->nodeHomePartition(kn)!=this->rank){
      EXPECT_NE((double(0.1*kn)),node_field->at(kn));
    };
    if(pmesh->nodeHomePartition(kn)==this->rank){
      EXPECT_DOUBLE_EQ((double(0.1*kn)),node_field->at(kn));
    };
  };
  for(int kc=0;kc<pmesh->cellCount();kc++){
    if(pmesh->cellHomePartition(kc)!=this->rank){
      FAIL();
    };
    if(pmesh->cellHomePartition(kc)==this->rank){
      EXPECT_DOUBLE_EQ((double(-0.1*kc)),cell_field->at(kc));
    };
  };
};

REGISTER_TYPED_TEST_SUITE_P(CommunicatorTests,
                            Nothing,GetMyRankAndSize,BroadCastRanks,
                            BroadCastRanksAsDoubles,BroadCastRanksAsI64,
                            DoubleRoundRobin,IntRoundRobin,Int64RoundRobin,
                            BroadCastArrayI32,BroadCastArrayI64,BroadCastArrayF64,
                            DoubleArrayRoundRobin,I32ArrayRoundRobin,I64ArrayRoundRobin,
                            HaloExchange,AllGatherF64,AllGatherI64,AllGatherI32);
