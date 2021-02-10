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

#include "Gateway.hpp"
#include "Communicator.hpp"
#include "Array.hpp"
#include "Mesh.hpp"
#include "Field.hpp"

#include "PrePartitionedStubHexMesh.hpp"

#include "gtest/gtest.h"

using testing::Types;

template<class ST>
std::shared_ptr<ST> makeShared(){
  return std::move(std::make_shared<ST>(ST()));
};

template <class I,class C,class MF,class FF,template<typename T>class AF>
struct GatewayTypeDefinitions
{
  typedef I Impl;
  typedef C Comm;
  typedef MF MeshFactory;
  typedef FF FieldFactory;
  typedef AF<double> DblFactory;
};

template<class T>
class GatewayTests : public ::testing::Test {
  protected:
    void SetUp() override {};
    void TearDown() override {};
    std::shared_ptr<PanNDE::Communicator> comm=makeShared<typename T::Comm>();
    std::shared_ptr<PanNDE::MeshFactory> meshmfg=makeShared<typename T::MeshFactory>();
    std::shared_ptr<PanNDE::FieldFactory> fieldmfg=makeShared<typename T::FieldFactory>();
    std::shared_ptr<Controller::Gateway> gateway=nullptr;
    std::shared_ptr<PanNDE::ArrayFactory<double>>f64mfg=makeShared<typename T::DblFactory>();
    std::string read_in_file="testing/angle_case_v2.vtu";

    void makeSerialGateway(){
      gateway=makeShared<typename T::Impl>();
    };
    void makeParallelGateway(){
      ASSERT_NE(nullptr,comm);
      gateway=std::make_shared<typename T::Impl>(typename T::Impl(comm));
    };
    void openFileWithGateway(){
      gateway->open(read_in_file);
    }
};

TYPED_TEST_SUITE_P(GatewayTests);

TYPED_TEST_P(GatewayTests,Nothing){};

TYPED_TEST_P(GatewayTests,CreateSerial){
  this->makeSerialGateway();
  EXPECT_NE(nullptr,this->gateway);
};

TYPED_TEST_P(GatewayTests,CreateParallel){
  this->makeParallelGateway();
  EXPECT_NE(nullptr,this->gateway);
};

TYPED_TEST_P(GatewayTests,GetWriteTimes){
  this->makeSerialGateway();
  this->openFileWithGateway();
  auto writetimes=this->gateway->getArray("write_times",this->f64mfg);
  ASSERT_NE(nullptr,writetimes);
  EXPECT_LT(0,writetimes->size());
  for(int k=1;k<writetimes->size();k++){EXPECT_LT(writetimes->at(k-1),writetimes->at(k));};
};

TYPED_TEST_P(GatewayTests,ParallelGetWriteTimes){
  this->makeParallelGateway();
  this->openFileWithGateway();
  auto writetimes=this->gateway->getArray("write_times",this->f64mfg);
  ASSERT_NE(nullptr,writetimes);
  EXPECT_LT(0,writetimes->size());
  for(int k=1;k<writetimes->size();k++){EXPECT_LT(writetimes->at(k-1),writetimes->at(k));};
};

TYPED_TEST_P(GatewayTests,GetTimeStep){
  this->makeSerialGateway();
  this->openFileWithGateway();
  double dt=this->gateway->getValue("dt");
  EXPECT_LT(0,dt);
};

TYPED_TEST_P(GatewayTests,ParallelGetTimeStep){
  this->makeParallelGateway();
  this->openFileWithGateway();
  double dt=this->gateway->getValue("dt");
  EXPECT_LT(0,dt);
};

TYPED_TEST_P(GatewayTests,GetMesh){
  this->makeSerialGateway();
  this->openFileWithGateway();
  auto mesh=this->gateway->getMesh(this->meshmfg);
  ASSERT_NE(nullptr,mesh);
  ASSERT_LT(0,mesh->nodeCount());
  ASSERT_LT(0,mesh->cellCount());
  double aabb[2][3];
  double pt[3];
  mesh->nodeCoordinate(0,pt);
  for(int k=0;k<3;k++){aabb[0][k]=pt[k];aabb[1][k]=pt[k];}
  for(int k=1;k<mesh->nodeCount();k++){
    mesh->nodeCoordinate(k,pt);
    for(int ka=0;ka<3;ka++){
      aabb[0][ka]=std::min(aabb[0][ka],pt[ka]);
      aabb[1][ka]=std::max(aabb[1][ka],pt[ka]);
    };
  };
  for(int ka=0;ka<3;ka++){EXPECT_LT(aabb[0][ka],aabb[1][ka]);};
};

TYPED_TEST_P(GatewayTests,ParallelGetMesh){
  this->makeParallelGateway();
  this->openFileWithGateway();
  auto mesh=this->gateway->getMesh(this->meshmfg);
  if(0==this->comm->getProcessId()){
    ASSERT_NE(nullptr,mesh);
    ASSERT_LT(0,mesh->nodeCount());
    ASSERT_LT(0,mesh->cellCount());
    double aabb[2][3];
    double pt[3];
    mesh->nodeCoordinate(0,pt);
    for(int k=0;k<3;k++){aabb[0][k]=pt[k];aabb[1][k]=pt[k];}
    for(int k=1;k<mesh->nodeCount();k++){
      mesh->nodeCoordinate(k,pt);
      for(int ka=0;ka<3;ka++){
        aabb[0][ka]=std::min(aabb[0][ka],pt[ka]);
        aabb[1][ka]=std::max(aabb[1][ka],pt[ka]);
      };
    };
    for(int ka=0;ka<3;ka++){EXPECT_LT(aabb[0][ka],aabb[1][ka]);};
  }else{ASSERT_EQ(nullptr,mesh);};
};

TYPED_TEST_P(GatewayTests,GetFields){
  this->makeSerialGateway();
  this->openFileWithGateway();
  auto mesh=this->gateway->getMesh(this->meshmfg);
  ASSERT_NE(nullptr,mesh);
  ASSERT_LT(0,mesh->nodeCount());
  ASSERT_LT(0,mesh->cellCount());

  auto density=this->gateway->getField("density",this->fieldmfg);
  auto C11=this->gateway->getField("C11",this->fieldmfg);
  ASSERT_NE(nullptr,density);
  ASSERT_NE(nullptr,C11);
  ASSERT_EQ(PanNDE::Field::NODE,density->type());
  ASSERT_EQ(PanNDE::Field::CELL,C11->type());

  for(int k=0;k<density->size();k++){EXPECT_LT(0.,density->at(k));};
  for(int k=0;k<C11->size();k++){EXPECT_LT(0.,C11->at(k));};
};

TYPED_TEST_P(GatewayTests,ParallelGetFields){
  this->makeParallelGateway();
  this->openFileWithGateway();
  auto mesh=this->gateway->getMesh(this->meshmfg);
  if(0==this->comm->getProcessId()){
    ASSERT_NE(nullptr,mesh);
    ASSERT_LT(0,mesh->nodeCount());
    ASSERT_LT(0,mesh->cellCount());
  };

  auto density=this->gateway->getField("density",this->fieldmfg);
  auto C11=this->gateway->getField("C11",this->fieldmfg);

  if(0==this->comm->getProcessId()){
    ASSERT_NE(nullptr,density);
    ASSERT_NE(nullptr,C11);
    ASSERT_EQ(PanNDE::Field::NODE,density->type());
    ASSERT_EQ(PanNDE::Field::CELL,C11->type());

    for(int k=0;k<density->size();k++){EXPECT_LT(0.,density->at(k));};
    for(int k=0;k<C11->size();k++){EXPECT_LT(0.,C11->at(k));};
  }else{
    EXPECT_EQ(nullptr,density);
    EXPECT_EQ(nullptr,C11);
  };
};

TYPED_TEST_P(GatewayTests,WriteToFile){
  this->makeSerialGateway();
  this->openFileWithGateway();
  auto mesh=this->gateway->getMesh(this->meshmfg);
  if(0==this->comm->getProcessId()){
    ASSERT_NE(nullptr,mesh);
    ASSERT_LT(0,mesh->nodeCount());
    ASSERT_LT(0,mesh->cellCount());
  };

  auto density=this->gateway->getField("density",this->fieldmfg);
  auto C11=this->gateway->getField("C11",this->fieldmfg);
  ASSERT_NE(nullptr,density);
  ASSERT_NE(nullptr,C11);

  double dt=this->gateway->getValue("dt");
  auto writetimes=this->gateway->getArray("write_times",this->f64mfg);
  ASSERT_LT(0,dt);
  ASSERT_NE(nullptr,writetimes);
  ASSERT_LT(0,writetimes->size());
  
  if(0==this->comm->getProcessId()){
    auto meta=this->f64mfg->makeManagedDataBundle();
    auto solution=this->fieldmfg->makeEmptyManagedFieldBundle();

    meta->emplaceScalar("dt",dt);
    meta->emplaceArray("write_times",writetimes);
    solution->mesh()=mesh;
    solution->emplaceField("density",density);
    solution->emplaceField("C11",C11);

    this->gateway->writeSolution("data/serial_testwrite",solution,meta,1);
  };
};

TYPED_TEST_P(GatewayTests,ParallelWriteToFile){
  this->makeParallelGateway();
  auto mesh=std::make_shared<Stubs::PrePartitionedStubHexMesh>(
                Stubs::PrePartitionedStubHexMesh(this->comm->getProcessId(),
                                                 this->comm->getNumberOfProcesses()));
  double dt=0.1;
  auto writetimes=this->f64mfg->makeManagedArray();
  writetimes->resize(4);
  for(int k=0;k<4;k++){writetimes->at(k)=dt*(double(k+1));};

  auto density=this->fieldmfg->makeManagedField(mesh,PanNDE::Field::NODE);
  auto C11=this->fieldmfg->makeManagedField(mesh,PanNDE::Field::CELL);

  for(int k=0;k<density->size();k++){density->at(k)=0.1*(double(k+1));};
  for(int k=0;k<C11->size();k++){C11->at(k)=-0.1*(double(k+1));};

  auto meta=this->f64mfg->makeManagedDataBundle();
  meta->emplaceScalar("dt",dt);
  meta->emplaceArray("write_times",writetimes);
  
  auto solution=this->fieldmfg->makeEmptyManagedFieldBundle();
  solution->mesh()=mesh;
  solution->emplaceField("density",density);
  solution->emplaceField("C11",C11);

  this->gateway->writeSolution("data/parallel_testwrite",solution,meta,1);
};

REGISTER_TYPED_TEST_SUITE_P(GatewayTests,
                            Nothing,CreateSerial,CreateParallel,GetWriteTimes,
                            ParallelGetWriteTimes,GetTimeStep,ParallelGetTimeStep,
                            GetMesh,ParallelGetMesh,GetFields,ParallelGetFields,
                            WriteToFile,ParallelWriteToFile);
