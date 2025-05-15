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
#include <cstdint>
#include <cmath>

#include "Model.hpp"
#include "HostRSGTDUT.hpp"
#include "ModelTests.hpp"

#include "gtest/gtest.h"

//Because the communicator will need to work with many formulations of the PanNDE 
//  configured data, the HostData module will be included, developers should add to the 
//  Types typedef as more *Data modules are developed.
#include "Array.hpp"
#include "Communicator.hpp"

#include "HostArray.hpp"
#include "MPICommunicator.hpp"

#include "PlateBuilder.hpp"
#include "HexMesh.hpp"
#include "HostField.hpp"
#include "internal/RightHex.hpp"


typedef Types<ModelTypeDefinitions<HostSolver::HostRSGTDUT,
                                   NetMPI::MPICommunicator,
                                   HostData::HostArrayFactory>>HostRSGTDUTTypes;

INSTANTIATE_TYPED_TEST_SUITE_P(HostRSGTDUTTests,    // Instance name
                               ModelTests,       // Test case name
                               HostRSGTDUTTypes);   // Type list

class RightHexTests : public ::testing::Test {
  protected:
    void SetUp() override {};
    void TearDown() override {};

    template<class T>
    inline std::shared_ptr<T> makeShared(){
      return std::move(std::make_shared<T>(T()));
    };
    std::shared_ptr<PanNDE::Mesh> makeUnitCell(){
      auto mesh_maker=makeShared<HostData::HexMeshFactory>();
      auto plate_maker=makeShared<Stubs::PlateBuilder>();
      auto mesh=plate_maker->makeMesh(mesh_maker,1.0,1.0,0.25,1,1,false);
      return mesh;
    };
    std::shared_ptr<PanNDE::Mesh> makeCubeAroundNode(){
      auto mesh_maker=makeShared<HostData::HexMeshFactory>();
      auto plate_maker=makeShared<Stubs::PlateBuilder>();
      auto mesh=plate_maker->makeMesh(mesh_maker,1.0,1.0,0.25,2,2,false);
      return mesh;
    };
    std::shared_ptr<PanNDE::Field> makeNodeField(std::shared_ptr<PanNDE::Mesh> mesh){
      auto field_maker=makeShared<HostData::HostFieldFactory>();
      auto field=field_maker->makeManagedField(mesh,PanNDE::Field::NODE);
      return field;
    };
    std::shared_ptr<PanNDE::Field> makeCellField(std::shared_ptr<PanNDE::Mesh> mesh){
      auto field_maker=makeShared<HostData::HostFieldFactory>();
      auto field=field_maker->makeManagedField(mesh,PanNDE::Field::CELL);
      return field;
    };
    std::shared_ptr<HostSolver::RightHex> makeDifferentiator(std::shared_ptr<PanNDE::Mesh> mesh){
      auto ddqi=std::make_shared<HostSolver::RightHex>(HostSolver::RightHex(mesh));
      return ddqi;
    };
};

TEST_F(RightHexTests,Nothing){};

TEST_F(RightHexTests,UnitCellMesh){
  auto mesh=this->makeUnitCell();
  EXPECT_EQ(1,mesh->cellCount());
  EXPECT_EQ(8,mesh->nodeCount());
  double pt[3];
  int box[8];
  mesh->cell(0,box);
  mesh->nodeCoordinate(box[0],pt);
  double expected[3]={0.,0.,0.};
  for(int kd=0;kd<3;kd++){EXPECT_DOUBLE_EQ(expected[kd],pt[kd]) << "box: 0, " << kd;};
  mesh->nodeCoordinate(box[1],pt);
  expected[0]=1.;expected[1]=0.;expected[2]=0.;
  for(int kd=0;kd<3;kd++){EXPECT_DOUBLE_EQ(expected[kd],pt[kd]) << "box: 1, " << kd;};
  mesh->nodeCoordinate(box[2],pt);
  expected[0]=1.;expected[1]=1.;expected[2]=0.;
  for(int kd=0;kd<3;kd++){EXPECT_DOUBLE_EQ(expected[kd],pt[kd]) << "box: 2, " << kd;};
  mesh->nodeCoordinate(box[3],pt);
  expected[0]=0.;expected[1]=1.;expected[2]=0.;
  for(int kd=0;kd<3;kd++){EXPECT_DOUBLE_EQ(expected[kd],pt[kd]) << "box: 3, " << kd;};
  mesh->nodeCoordinate(box[4],pt);
  expected[0]=0.;expected[1]=0.;expected[2]=0.25;
  for(int kd=0;kd<3;kd++){EXPECT_DOUBLE_EQ(expected[kd],pt[kd]) << "box: 4, " << kd;};
  mesh->nodeCoordinate(box[5],pt);
  expected[0]=1.;expected[1]=0.;expected[2]=0.25;
  for(int kd=0;kd<3;kd++){EXPECT_DOUBLE_EQ(expected[kd],pt[kd]) << "box: 5, " << kd;};
  mesh->nodeCoordinate(box[6],pt);
  expected[0]=1.;expected[1]=1.;expected[2]=0.25;
  for(int kd=0;kd<3;kd++){EXPECT_DOUBLE_EQ(expected[kd],pt[kd]) << "box: 6, " << kd;};
  mesh->nodeCoordinate(box[7],pt);
  expected[0]=0.;expected[1]=1.;expected[2]=0.25;
  for(int kd=0;kd<3;kd++){EXPECT_DOUBLE_EQ(expected[kd],pt[kd]) << "box: 7, " << kd;};
};
TEST_F(RightHexTests,NodeFieldOnUnitCellMesh){
  auto mesh=this->makeUnitCell();
  auto nodefield=this->makeNodeField(mesh);
  EXPECT_EQ(8,nodefield->size());
};
TEST_F(RightHexTests,CheckDerivativesOfCoordinates){
  auto mesh=this->makeUnitCell();
  auto nodefield=this->makeNodeField(mesh);
  auto ddqi=this->makeDifferentiator(mesh);
  double pt[3];
  for(int k=0;k<mesh->nodeCount();k++){
    mesh->nodeCoordinate(k,pt);
    nodefield->at(k)=pt[0];
  };
  EXPECT_DOUBLE_EQ(1.,ddqi->DDx_cc(0,nodefield));
  EXPECT_DOUBLE_EQ(0.,ddqi->DDy_cc(0,nodefield));
  EXPECT_DOUBLE_EQ(0.,ddqi->DDz_cc(0,nodefield));
  for(int k=0;k<mesh->nodeCount();k++){
    mesh->nodeCoordinate(k,pt);
    nodefield->at(k)=pt[1];
  };
  EXPECT_DOUBLE_EQ(0.,ddqi->DDx_cc(0,nodefield));
  EXPECT_DOUBLE_EQ(1.,ddqi->DDy_cc(0,nodefield));
  EXPECT_DOUBLE_EQ(0.,ddqi->DDz_cc(0,nodefield));
  for(int k=0;k<mesh->nodeCount();k++){
    mesh->nodeCoordinate(k,pt);
    nodefield->at(k)=pt[2];
  };
  EXPECT_DOUBLE_EQ(0.,ddqi->DDx_cc(0,nodefield));
  EXPECT_DOUBLE_EQ(0.,ddqi->DDy_cc(0,nodefield));
  EXPECT_DOUBLE_EQ(1.,ddqi->DDz_cc(0,nodefield));
};
TEST_F(RightHexTests,CheckDerivativesOfDistance){
  auto mesh=this->makeUnitCell();
  auto nodefield=this->makeNodeField(mesh);
  auto ddqi=this->makeDifferentiator(mesh);
  double pt[3];
  for(int k=0;k<mesh->nodeCount();k++){
    mesh->nodeCoordinate(k,pt);
    nodefield->at(k)=sqrt(pt[0]*pt[0]+pt[1]*pt[1]+pt[2]*pt[2]);
  };
  double ddx=0.25*(sqrt(2.)+(sqrt(33.)/4.)-(1./4.))/1.;
  double ddy=ddx;
  double ddz=0.25*((sqrt(17.)/2.+sqrt(33.)/4.+1./4.)-(2.+sqrt(2.)))/0.25;
  EXPECT_DOUBLE_EQ(ddx,ddqi->DDx_cc(0,nodefield));
  EXPECT_DOUBLE_EQ(ddy,ddqi->DDy_cc(0,nodefield));
  EXPECT_DOUBLE_EQ(ddz,ddqi->DDz_cc(0,nodefield));
};
TEST_F(RightHexTests,CheckDerivativesOfNodeIndex){
  auto mesh=this->makeUnitCell();
  auto nodefield=this->makeNodeField(mesh);
  auto ddqi=this->makeDifferentiator(mesh);
  for(int k=0;k<mesh->nodeCount();k++){
    nodefield->at(k)=double(k);
  };
  double ddx=0.25*((4.+5.+6.+7.)-(0.+1.+2.+3.))/1.;
  double ddy=0.25*((2.+3.+6.+7.)-(0.+1.+4.+5.))/1.;
  double ddz=0.25*((1.+3.+5.+7.)-(0.+2.+4.+6.))/0.25;
  EXPECT_DOUBLE_EQ(ddx,ddqi->DDx_cc(0,nodefield));
  EXPECT_DOUBLE_EQ(ddy,ddqi->DDy_cc(0,nodefield));
  EXPECT_DOUBLE_EQ(ddz,ddqi->DDz_cc(0,nodefield));
};
TEST_F(RightHexTests,Check8pack){
  auto mesh=this->makeCubeAroundNode();
  EXPECT_EQ(8,mesh->cellCount());
  EXPECT_EQ(27,mesh->nodeCount());
  double pt[3];
  mesh->nodeCoordinate(13,pt);
  EXPECT_DOUBLE_EQ(0.5,pt[0]);
  EXPECT_DOUBLE_EQ(0.5,pt[1]);
  EXPECT_DOUBLE_EQ(0.125,pt[2]);
};
TEST_F(RightHexTests,CellFieldOn8Pack){
  auto mesh=this->makeCubeAroundNode();
  auto cellfield=this->makeCellField(mesh);
  EXPECT_EQ(8,cellfield->size());
};
TEST_F(RightHexTests,CheckDerivativesOfCellCoordinates){
  auto mesh=this->makeCubeAroundNode();
  auto cellfield=this->makeCellField(mesh);
  auto ddqi=this->makeDifferentiator(mesh);
  double pt[3];
  int box[8];

  EXPECT_DOUBLE_EQ(0.5,ddqi->getNodeDx(13));
  EXPECT_DOUBLE_EQ(0.5,ddqi->getNodeDy(13));
  EXPECT_DOUBLE_EQ(0.125,ddqi->getNodeDz(13));

  for(int idx=0;idx<3;idx++){
    for(int kc=0;kc<8;kc++){
      double ptc[3]={0.,0.,0.};
      mesh->cell(kc,box);
      for(int kb=0;kb<8;kb++){
        mesh->nodeCoordinate(box[kb],pt);
        //printf("%i-%i: %f %f %f\n",kb,box[kb],pt[0],pt[1],pt[2]);
        for(int kd=0;kd<3;kd++){ptc[kd]+=pt[kd]/8.;};
      };
      //printf("%i: %f %f %f\n",kc,ptc[0],ptc[1],ptc[2]);
      cellfield->at(kc)=ptc[idx];
    };
    EXPECT_DOUBLE_EQ(double(0==idx),ddqi->DDx_nc(13,cellfield));
    EXPECT_DOUBLE_EQ(double(1==idx),ddqi->DDy_nc(13,cellfield));
    EXPECT_DOUBLE_EQ(double(2==idx),ddqi->DDz_nc(13,cellfield));
  };
};

