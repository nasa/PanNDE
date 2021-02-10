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
#include <cmath>
#include <cstdint>

#include "Array.hpp"
#include "Mesh.hpp"

#include "gtest/gtest.h"

using testing::Types;


template<class FT>
std::shared_ptr<FT> makeShared(){
  return std::move(std::make_shared<FT>(FT()));
};

template <class I,class F,template<typename T>class AF>
struct MeshTypeDefinitions
{
  typedef I Impl;
  typedef F Factory;
  typedef AF<int32_t> IArrayMfg;
  typedef AF<double> DArrayMfg;
  typedef AF<PanNDE::Mesh::Node> NArrayMfg;
  typedef AF<PanNDE::Mesh::Cell> CArrayMfg;
};

template<class T>
class MeshTests : public ::testing::Test {
  protected:
    void SetUp() override {
      mesh_factory=std::make_shared<typename T::Factory>(typename T::Factory());
    };
    void TearDown() override {};
    std::shared_ptr<PanNDE::MeshFactory> mesh_factory=nullptr;
};

template<class T>
class SingleCubicHex{
  public:
    std::shared_ptr<PanNDE::Mesh> makeCube(int simpartition=0){
      auto MMFG=makeShared<typename T::Factory>();
      auto NMFG=makeShared<typename T::NArrayMfg>();
      auto nodes=NMFG->makeManagedArray();
      auto CMFG=makeShared<typename T::CArrayMfg>();
      auto cells=CMFG->makeManagedArray();
      nodes->resize(8);
      cells->resize(1);
      double coord[3];
      for(int kb=0;kb<8;kb++){
        coord[0]=X[kb];coord[1]=Y[kb];coord[2]=Z[kb];
        nodes->at(kb)=MMFG->makeNode(coord,7-kb,4+simpartition);
      };
      cells->at(0)=MMFG->makeCell(cell,8,1,5+simpartition);
      auto mesh=MMFG->makeManagedMesh(nodes,cells,simpartition);
      return std::move(mesh);
    };
    void checkIsCube(std::shared_ptr<PanNDE::Mesh> mesh){
      double coord[3];
      int32_t box[8];
      mesh->cell(0,box);
      auto abox=mesh->cell(0);
      for(int kb=0;kb<8;kb++){
        mesh->nodeCoordinate(kb,coord);
        EXPECT_EQ(X[kb],coord[0]);
        EXPECT_EQ(Y[kb],coord[1]);
        EXPECT_EQ(Z[kb],coord[2]);
        auto r=mesh->nodeCoordinate(kb);
        EXPECT_EQ(X[kb],r->at(0));
        EXPECT_EQ(Y[kb],r->at(1));
        EXPECT_EQ(Z[kb],r->at(2));
        EXPECT_EQ(cell[kb],box[kb]);
        EXPECT_EQ(cell[kb],abox->at(kb));
      };
    };
    void checkConnectedCells(std::shared_ptr<PanNDE::Mesh> mesh){
      int32_t box[8];
      mesh->connectedCells(0,box);
      auto abox=mesh->connectedCells(0);
      for(int k=0;k<8;k++){EXPECT_EQ(box[k],(6==k)?0:-1);EXPECT_EQ(abox->at(k),(6==k)?0:-1);};
      mesh->connectedCells(1,box);
      abox=mesh->connectedCells(1);
      for(int k=0;k<8;k++){EXPECT_EQ(box[k],(7==k)?0:-1);EXPECT_EQ(abox->at(k),(7==k)?0:-1);};
      mesh->connectedCells(2,box);
      abox=mesh->connectedCells(2);
      for(int k=0;k<8;k++){EXPECT_EQ(box[k],(4==k)?0:-1);EXPECT_EQ(abox->at(k),(4==k)?0:-1);};
      mesh->connectedCells(3,box);
      abox=mesh->connectedCells(3);
      for(int k=0;k<8;k++){EXPECT_EQ(box[k],(5==k)?0:-1);EXPECT_EQ(abox->at(k),(5==k)?0:-1);};
      mesh->connectedCells(4,box);
      abox=mesh->connectedCells(4);
      for(int k=0;k<8;k++){EXPECT_EQ(box[k],(2==k)?0:-1);EXPECT_EQ(abox->at(k),(2==k)?0:-1);};
      mesh->connectedCells(5,box);
      abox=mesh->connectedCells(5);
      for(int k=0;k<8;k++){EXPECT_EQ(box[k],(3==k)?0:-1);EXPECT_EQ(abox->at(k),(3==k)?0:-1);};
      mesh->connectedCells(6,box);
      abox=mesh->connectedCells(6);
      for(int k=0;k<8;k++){EXPECT_EQ(box[k],(0==k)?0:-1);EXPECT_EQ(abox->at(k),(0==k)?0:-1);};
      mesh->connectedCells(7,box);
      abox=mesh->connectedCells(7);
      for(int k=0;k<8;k++){EXPECT_EQ(box[k],(1==k)?0:-1);EXPECT_EQ(abox->at(k),(1==k)?0:-1);};
    };
    SingleCubicHex(){
      X[0]=0.0;Y[0]=0.0;Z[0]=0.0;
      X[1]=1.0;Y[1]=0.0;Z[1]=0.0;
      X[2]=1.0;Y[2]=1.0;Z[2]=0.0;
      X[3]=0.0;Y[3]=1.0;Z[3]=0.0;
      X[4]=0.0;Y[4]=0.0;Z[4]=1.0;
      X[5]=1.0;Y[5]=0.0;Z[5]=1.0;
      X[6]=1.0;Y[6]=1.0;Z[6]=1.0;
      X[7]=0.0;Y[7]=1.0;Z[7]=1.0;
      for(int kb=0;kb<8;kb++){cell[kb]=kb;};
    };
  private:
    double X[8];
    double Y[8];
    double Z[8];
    int32_t cell[8];
};

template<typename T>
class FlatPlateHexMesh {
  public:
    std::shared_ptr<PanNDE::Mesh> makeFlatPlateHexMesh(std::array<double,3> ds={0.1,0.1,0.1},
                                                       std::array<int,3> Nn={11,11,5},
                                                       int simpartition=0){
      auto MMFG=makeShared<typename T::Factory>();
      auto NMFG=makeShared<typename T::NArrayMfg>();
      auto nodes=NMFG->makeManagedArray();
      auto CMFG=makeShared<typename T::CArrayMfg>();
      auto cells=CMFG->makeManagedArray();
      int Nnt=Nn[0]*Nn[1]*Nn[2];
      int Nct=(Nn[0]-1)*(Nn[1]-1)*(Nn[2]-1);
      int idx;
      nodes->resize(Nnt);
      cells->resize(Nct);
      double coord[3];
      for(int kx=0;kx<Nn[0];kx++){
        for(int ky=0;ky<Nn[1];ky++){
          for(int kz=0;kz<Nn[2];kz++){
            coord[0]=ds[0]*kx;coord[1]=ds[1]*ky;coord[2]=ds[2]*kz;
            idx=kz+Nn[2]*(ky+Nn[1]*kx);
            nodes->at(idx)=MMFG->makeNode(coord,idx,simpartition);
          };
        };
      };
      int32_t box[8];
      for(int kx=0;kx<(Nn[0]-1);kx++){
        for(int ky=0;ky<(Nn[1]-1);ky++){
          for(int kz=0;kz<(Nn[2]-1);kz++){
            idx=kz+(Nn[2]-1)*(ky+(Nn[1]-1)*kx);

            box[0]=(kz)+(Nn[2]*((ky)+Nn[1]*(kx)));
            box[1]=(kz)+(Nn[2]*((ky)+Nn[1]*(kx+1)));
            box[2]=(kz)+(Nn[2]*((ky+1)+Nn[1]*(kx+1)));
            box[3]=(kz)+(Nn[2]*((ky+1)+Nn[1]*(kx)));
            box[4]=(kz+1)+(Nn[2]*((ky)+Nn[1]*(kx)));
            box[5]=(kz+1)+(Nn[2]*((ky)+Nn[1]*(kx+1)));
            box[6]=(kz+1)+(Nn[2]*((ky+1)+Nn[1]*(kx+1)));
            box[7]=(kz+1)+(Nn[2]*((ky+1)+Nn[1]*(kx)));
            

            cells->at(idx)=MMFG->makeCell(box,8,idx,simpartition);
          };
        };
      };
      auto mesh=MMFG->makeManagedMesh(nodes,cells,simpartition);
      return std::move(mesh);
    };
    void checkMeshConnections(std::shared_ptr<PanNDE::Mesh> mesh,std::array<int,3> Nn){
      for(int kx=0;kx<(Nn[0]-1);kx++){
        for(int ky=0;ky<(Nn[1]-1);ky++){
          for(int kz=0;kz<(Nn[2]-1);kz++){
            int cidx=kz+(Nn[2]-1)*(ky+(Nn[1]-1)*kx);
            auto cbox=mesh->cell(cidx);
            EXPECT_EQ((kz)+(Nn[2]*((ky)+Nn[1]*(kx))),cbox->at(0));
            EXPECT_EQ((kz)+(Nn[2]*((ky)+Nn[1]*(kx+1))),cbox->at(1));
            EXPECT_EQ((kz)+(Nn[2]*((ky+1)+Nn[1]*(kx+1))),cbox->at(2));
            EXPECT_EQ((kz)+(Nn[2]*((ky+1)+Nn[1]*(kx))),cbox->at(3));
            EXPECT_EQ((kz+1)+(Nn[2]*((ky)+Nn[1]*(kx))),cbox->at(4));
            EXPECT_EQ((kz+1)+(Nn[2]*((ky)+Nn[1]*(kx+1))),cbox->at(5));
            EXPECT_EQ((kz+1)+(Nn[2]*((ky+1)+Nn[1]*(kx+1))),cbox->at(6));
            EXPECT_EQ((kz+1)+(Nn[2]*((ky+1)+Nn[1]*(kx))),cbox->at(7));
            
            auto nbox=mesh->connectedCells(cbox->at(0));EXPECT_EQ(cidx,nbox->at(6));
            nbox=mesh->connectedCells(cbox->at(1));EXPECT_EQ(cidx,nbox->at(7));
            nbox=mesh->connectedCells(cbox->at(2));EXPECT_EQ(cidx,nbox->at(4));
            nbox=mesh->connectedCells(cbox->at(3));EXPECT_EQ(cidx,nbox->at(5));
            nbox=mesh->connectedCells(cbox->at(4));EXPECT_EQ(cidx,nbox->at(2));
            nbox=mesh->connectedCells(cbox->at(5));EXPECT_EQ(cidx,nbox->at(3));
            nbox=mesh->connectedCells(cbox->at(6));EXPECT_EQ(cidx,nbox->at(0));
            nbox=mesh->connectedCells(cbox->at(7));EXPECT_EQ(cidx,nbox->at(1));
          };
        };
      };
    };
};

TYPED_TEST_SUITE_P(MeshTests);

TYPED_TEST_P(MeshTests,Nothing){};

TYPED_TEST_P(MeshTests,FactoriesExist){
  auto NMFG=makeShared<typename TypeParam::NArrayMfg>();
  auto CMFG=makeShared<typename TypeParam::CArrayMfg>();
  EXPECT_NE(nullptr,CMFG);
  EXPECT_NE(nullptr,NMFG);
  EXPECT_NE(nullptr,this->mesh_factory);
};

TYPED_TEST_P(MeshTests,CubeCheck){
  SingleCubicHex<TypeParam> hex;
  auto mesh=hex.makeCube();
  ASSERT_NE(nullptr,mesh);
  hex.checkIsCube(mesh);
};
TYPED_TEST_P(MeshTests,CubeCheckNodeToCell){
  SingleCubicHex<TypeParam> hex;
  auto mesh=hex.makeCube();
  ASSERT_NE(nullptr,mesh);
  hex.checkIsCube(mesh);
  hex.checkConnectedCells(mesh);
};
TYPED_TEST_P(MeshTests,CheckCubeCounts){
  SingleCubicHex<TypeParam> hex;
  auto mesh=hex.makeCube();
  EXPECT_EQ(1,mesh->cellCount());
  EXPECT_EQ(8,mesh->nodeCount());
};
TYPED_TEST_P(MeshTests,CheckPartition){
  SingleCubicHex<TypeParam> hex;
  auto mesh=hex.makeCube(5);
  EXPECT_EQ(5,mesh->partitionId());
};
TYPED_TEST_P(MeshTests,CheckGlobalIds){
  SingleCubicHex<TypeParam> hex;
  auto mesh=hex.makeCube();
  for(auto k=0;k<8;k++){EXPECT_EQ(7-k,mesh->globalNodeId(k));};
  EXPECT_EQ(1,mesh->globalCellId(0));
};
TYPED_TEST_P(MeshTests,CheckOwnerIds){
  SingleCubicHex<TypeParam> hex;
  auto mesh=hex.makeCube();
  for(auto k=0;k<8;k++){EXPECT_EQ(4,mesh->nodeHomePartition(k));};
  EXPECT_EQ(5,mesh->cellHomePartition(0));
};
TYPED_TEST_P(MeshTests,CopyRawPtr){
  SingleCubicHex<TypeParam> hex;
  auto mesh=hex.makeCube();
  auto newmesh=makeShared<typename TypeParam::Impl>();
  newmesh->copy(mesh.get());
  hex.checkIsCube(newmesh);
  hex.checkConnectedCells(newmesh);
};
TYPED_TEST_P(MeshTests,CopyManagedPtr){
  SingleCubicHex<TypeParam> hex;
  auto mesh=hex.makeCube();
  auto newmesh=makeShared<typename TypeParam::Impl>();
  newmesh->copy(mesh);
  hex.checkIsCube(newmesh);
  hex.checkConnectedCells(newmesh);
};

TYPED_TEST_P(MeshTests,MakePlate){
  FlatPlateHexMesh<TypeParam> plate;
  auto mesh=plate.makeFlatPlateHexMesh();
};
TYPED_TEST_P(MeshTests,CheckPlateConnections){
  FlatPlateHexMesh<TypeParam> plate;
  std::array<double,3> ds={0.1,0.1,0.1};
  std::array<int,3> Nn={11,11,5};
  auto mesh=plate.makeFlatPlateHexMesh(ds,Nn);
  plate.checkMeshConnections(mesh,Nn);
};


REGISTER_TYPED_TEST_SUITE_P(MeshTests,
                            Nothing,FactoriesExist,CubeCheck,
                            CubeCheckNodeToCell,CheckCubeCounts,
                            CheckPartition,CheckGlobalIds,
                            CheckOwnerIds,CopyRawPtr,CopyManagedPtr,
                            MakePlate,CheckPlateConnections);
