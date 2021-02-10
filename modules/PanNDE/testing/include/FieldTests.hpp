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

#include "HostField.hpp"

#include "gtest/gtest.h"

using testing::Types;


template<class FT>
std::shared_ptr<FT> makeShared(){
  return std::move(std::make_shared<FT>(FT()));
};

template <class I,class F,class MF>
struct FieldTypeDefinitions
{
  typedef I Impl;
  typedef F Factory;
  typedef MF MeshFactory;
};

template<class T>
class FieldTests : public ::testing::Test {
  protected:
    void SetUp() override {};
    void TearDown() override {};
    std::shared_ptr<PanNDE::FieldFactory> fieldMfg=makeShared<typename T::Factory>();
    std::shared_ptr<PanNDE::MeshFactory> meshMfg=makeShared<typename T::MeshFactory>();

    void checkFieldFactory(){ASSERT_NE(nullptr,fieldMfg);};
    void checkMeshFactory(){ASSERT_NE(nullptr,meshMfg);};
    std::shared_ptr<PanNDE::Mesh> makeMesh(std::array<double,3> ds={0.1,0.1,0.1},
                                           std::array<int,3> Nn={11,11,5}){
      std::vector<PanNDE::Mesh::Node> nodes;nodes.resize(Nn.at(0)*Nn.at(1)*Nn.at(2));
      std::vector<PanNDE::Mesh::Cell> cells;cells.resize((Nn.at(0)-1)*(Nn.at(1)-1)*(Nn.at(2)-1));
      for(int kx=0;kx<Nn.at(0);kx++){
        for(int ky=0;ky<Nn.at(1);ky++){
          for(int kz=0;kz<Nn.at(2);kz++){
            int idx=kz+Nn.at(2)*(ky+Nn.at(1)*kx);
            double coord[3]={kx*ds.at(0),ky*ds.at(1),kz*ds.at(2)};
            nodes.at(idx)=meshMfg->makeNode(coord,idx,0);
          };
        };
      };
      for(int kx=0;kx<(Nn.at(0)-1);kx++){
        for(int ky=0;ky<(Nn.at(1)-1);ky++){
          for(int kz=0;kz<(Nn.at(2)-1);kz++){
            int idx=kz+(Nn.at(2)-1)*(ky+(Nn.at(1)-1)*kx);
            int32_t box[8]={
              (kz)+(Nn.at(2)*((ky)+Nn.at(1)*(kx))),
              (kz)+(Nn.at(2)*((ky)+Nn.at(1)*(kx+1))),
              (kz)+(Nn.at(2)*((ky+1)+Nn.at(1)*(kx+1))),
              (kz)+(Nn.at(2)*((ky+1)+Nn.at(1)*(kx))),
              (kz+1)+(Nn.at(2)*((ky)+Nn.at(1)*(kx))),
              (kz+1)+(Nn.at(2)*((ky)+Nn.at(1)*(kx+1))),
              (kz+1)+(Nn.at(2)*((ky+1)+Nn.at(1)*(kx+1))),
              (kz+1)+(Nn.at(2)*((ky+1)+Nn.at(1)*(kx)))
            };
            cells.at(idx)=meshMfg->makeCell(box,8,idx,0);
          };
        };
      };
      auto mesh=meshMfg->makeManagedMesh(nodes.data(),nodes.size(),cells.data(),cells.size());
      return std::move(mesh);
    };
    void checkMapper(std::shared_ptr<PanNDE::Field> src_field,PanNDE::Field::FieldType tp){
      auto tgt_field=this->fieldMfg->makeManagedField(src_field->mesh(),tp);
      ASSERT_NE(nullptr,tgt_field);
      tgt_field->mapFrom(src_field);
      if(PanNDE::Field::NODE==tgt_field->type()){
        for(auto k=0;k<src_field->mesh()->nodeCount();k++){
          EXPECT_EQ(src_field->atNode(k),tgt_field->at(k));
        };
      }else{
        for(auto k=0;k<src_field->mesh()->cellCount();k++){
          EXPECT_EQ(src_field->atCell(k),tgt_field->at(k));
        };
      };
    };
};

TYPED_TEST_SUITE_P(FieldTests);

TYPED_TEST_P(FieldTests,Nothing){};

TYPED_TEST_P(FieldTests,CheckSetup){
  this->checkMeshFactory();
  auto mesh=this->makeMesh();
  ASSERT_NE(nullptr,mesh);
  this->checkFieldFactory();
};

TYPED_TEST_P(FieldTests,makeFields){
  auto mesh=this->makeMesh();
  auto cell_field=this->fieldMfg->makeManagedField(mesh,PanNDE::Field::CELL);
  ASSERT_NE(nullptr,cell_field);
  auto node_field=this->fieldMfg->makeManagedField(mesh,PanNDE::Field::NODE);
  ASSERT_NE(nullptr,node_field);
};
TYPED_TEST_P(FieldTests,checkFieldInfo){
  auto mesh=this->makeMesh();
  auto cell_field=this->fieldMfg->makeManagedField(mesh,PanNDE::Field::CELL);
  ASSERT_NE(nullptr,cell_field);
  auto node_field=this->fieldMfg->makeManagedField(mesh,PanNDE::Field::NODE);
  ASSERT_NE(nullptr,node_field);
  auto const_field=this->fieldMfg->makeManagedField(mesh,PanNDE::Field::CONSTANT);
  ASSERT_NE(nullptr,node_field);

  EXPECT_EQ(PanNDE::Field::CELL,cell_field->type());
  EXPECT_EQ(PanNDE::Field::NODE,node_field->type());
  EXPECT_EQ(PanNDE::Field::CONSTANT,const_field->type());

  EXPECT_EQ(mesh,cell_field->mesh());
  EXPECT_EQ(mesh,node_field->mesh());
  EXPECT_EQ(mesh,const_field->mesh());

  EXPECT_EQ(mesh->cellCount(),cell_field->size());
  EXPECT_EQ(mesh->nodeCount(),node_field->size());
  EXPECT_EQ(1,const_field->size());

  EXPECT_NE(nullptr,cell_field->data());
  EXPECT_NE(nullptr,node_field->data());
  EXPECT_NE(nullptr,const_field->data());
};
TYPED_TEST_P(FieldTests,checkFieldActions){
  auto mesh=this->makeMesh();
  auto cell_field=this->fieldMfg->makeManagedField(mesh,PanNDE::Field::CELL);
  ASSERT_NE(nullptr,cell_field);
  auto node_field=this->fieldMfg->makeManagedField(mesh,PanNDE::Field::NODE);
  ASSERT_NE(nullptr,node_field);
  auto const_field=this->fieldMfg->makeManagedField(mesh,PanNDE::Field::CONSTANT);
  ASSERT_NE(nullptr,node_field);

  for(auto k=0;k<cell_field->size();k++){cell_field->at(k)=((double)(k));};
  for(auto k=0;k<cell_field->size();k++){EXPECT_EQ((double)(k),cell_field->at(k));};
  for(auto k=0;k<node_field->size();k++){node_field->at(k)=-((double)(k));};
  for(auto k=0;k<node_field->size();k++){EXPECT_EQ(-(double)(k),node_field->at(k));};
  const_field->at(0)=0.1;
  EXPECT_EQ(0.1,const_field->at(345));
  const_field->at(45)=0.2;
  EXPECT_EQ(0.2,const_field->at(12));

  for(auto k=0;k<mesh->nodeCount();k++){
    auto coord=mesh->nodeCoordinate(k);
    node_field->at(k)=coord->at(0)+coord->at(1)+coord->at(2);
  };
  for(auto k=0;k<mesh->cellCount();k++){
    auto box=mesh->cell(k);
    double pt=0.;
    for(int kb=0;kb<box->size();kb++){
      auto coord=mesh->nodeCoordinate(box->at(kb));
      for(int kd=0;kd<3;kd++){pt+=coord->at(kd);};
    };
    pt=pt/(double(box->size()));
    cell_field->at(k)=pt;
  };
  for(auto k=0;k<mesh->cellCount();k++){
    EXPECT_FLOAT_EQ(cell_field->at(k),node_field->atCell(k));
    EXPECT_EQ(const_field->at(0),const_field->atCell(k));
  };

  for(auto k=0;k<mesh->cellCount();k++){cell_field->at(k)=1.0;};
  for(auto k=0;k<mesh->nodeCount();k++){node_field->at(k)=1.0;};
  const_field->at(0)=1.0;
  for(auto k=0;k<mesh->nodeCount();k++){
    auto box=mesh->connectedCells(k);
    int cnt=0;for(auto kb=0;kb<box->size();kb++){cnt+=(-1==box->at(kb))?0:1;};
    EXPECT_FLOAT_EQ(node_field->at(k)*(double(cnt)/double(box->size())),cell_field->atNode(k));
    EXPECT_FLOAT_EQ(node_field->at(k)*(double(cnt)/double(box->size())),const_field->atNode(k));
  };
};

TYPED_TEST_P(FieldTests,Mapping){
  auto mesh=this->makeMesh();
  auto cell_field=this->fieldMfg->makeManagedField(mesh,PanNDE::Field::CELL);
  ASSERT_NE(nullptr,cell_field);
  auto node_field=this->fieldMfg->makeManagedField(mesh,PanNDE::Field::NODE);
  ASSERT_NE(nullptr,node_field);
  auto const_field=this->fieldMfg->makeManagedField(mesh,PanNDE::Field::CONSTANT);
  ASSERT_NE(nullptr,node_field);
  for(auto k=0;k<mesh->nodeCount();k++){
    auto coord=mesh->nodeCoordinate(k);
    node_field->at(k)=coord->at(0)+coord->at(1)+coord->at(2);
  };
  for(auto k=0;k<mesh->cellCount();k++){
    auto box=mesh->cell(k);
    double pt=0.;
    for(int kb=0;kb<box->size();kb++){
      auto coord=mesh->nodeCoordinate(box->at(kb));
      for(int kd=0;kd<3;kd++){pt+=coord->at(kd);};
    };
    pt=pt/(double(box->size()));
    cell_field->at(k)=pt;
  };

  this->checkMapper(cell_field,PanNDE::Field::NODE);
  this->checkMapper(cell_field,PanNDE::Field::CELL);
  this->checkMapper(cell_field,PanNDE::Field::CONSTANT);
  this->checkMapper(node_field,PanNDE::Field::NODE);
  this->checkMapper(node_field,PanNDE::Field::CELL);
  this->checkMapper(node_field,PanNDE::Field::CONSTANT);
  this->checkMapper(const_field,PanNDE::Field::NODE);
  this->checkMapper(const_field,PanNDE::Field::CELL);
  this->checkMapper(const_field,PanNDE::Field::CONSTANT);
};
TYPED_TEST_P(FieldTests,makeFieldBundle){
  auto mesh=this->makeMesh();
  auto cell_field=this->fieldMfg->makeManagedField(mesh,PanNDE::Field::CELL);
  ASSERT_NE(nullptr,cell_field);
  auto node_field=this->fieldMfg->makeManagedField(mesh,PanNDE::Field::NODE);
  ASSERT_NE(nullptr,node_field);
  auto bundle=this->fieldMfg->makeEmptyManagedFieldBundle();
  bundle->mesh()=mesh;
  bundle->emplaceField("cell",cell_field);
  bundle->emplaceField("node",node_field);

  for(auto k=0;k<mesh->cellCount();k++){
    auto box=mesh->cell(k);
    double pt=0.;
    for(int kb=0;kb<box->size();kb++){
      auto coord=mesh->nodeCoordinate(box->at(kb));
      for(int kd=0;kd<3;kd++){pt+=coord->at(kd);};
    };
    pt=pt/(double(box->size()));
    cell_field->at(k)=pt;
  };
  for(auto k=0;k<mesh->nodeCount();k++){
    auto coord=mesh->nodeCoordinate(k);
    double pt=0.;
    for(int kd=0;kd<3;kd++){pt+=coord->at(kd);};
    node_field->at(k)=pt;
  };

  for(auto k=0;k<mesh->cellCount();k++){
    EXPECT_EQ(cell_field->at(k),bundle->field("cell")->at(k));
  };
  for(auto k=0;k<mesh->nodeCount();k++){
    EXPECT_EQ(node_field->at(k),bundle->field("node")->at(k));
  };
  EXPECT_EQ(2,cell_field.use_count());
  EXPECT_EQ(2,node_field.use_count());

  ASSERT_EQ(2,bundle->fieldCount());
  auto name=bundle->fieldName(0);
  EXPECT_STREQ("cell",name.c_str());
  name=bundle->fieldName(1);
  EXPECT_STREQ("node",name.c_str());
};

REGISTER_TYPED_TEST_SUITE_P(FieldTests,
                            Nothing,CheckSetup,makeFields,checkFieldInfo,
                            checkFieldActions,Mapping,makeFieldBundle);
