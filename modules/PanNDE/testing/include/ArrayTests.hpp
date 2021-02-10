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

#include "Array.hpp"

#include "gtest/gtest.h"

using testing::Types;

template <template<typename T>class I, 
          template<typename T>class B,
          template<typename T>class F,
          typename T>
struct ArrayTypeDefinitions
{
  typedef class I<T> Impl;
  typedef class B<T> Bundle;
  typedef class F<T> Factory;
  typedef T Dtype;
};

template<class T>
class ArrayTests : public ::testing::Test {
  protected:
    void SetUp() override {
      factory=std::make_shared<typename T::Factory>(typename T::Factory());
      array=factory->makeManagedArray();
    };
    void TearDown() override {};
    std::shared_ptr<PanNDE::Array<typename T::Dtype>> array=nullptr;
    std::shared_ptr<PanNDE::DataBundle<typename T::Dtype>> bundle=nullptr;
    std::shared_ptr<PanNDE::ArrayFactory<typename T::Dtype>> factory=nullptr;
    int N=10;
};

TYPED_TEST_SUITE_P(ArrayTests);

TYPED_TEST_P(ArrayTests,Nothing){};

TYPED_TEST_P(ArrayTests,Resize){
  // Inside a test, refer to TypeParam to get the type parameter.
  //TypeParam n = 0;
  //...
  ASSERT_NE(nullptr,this->array);
  this->array->resize(this->N);
  EXPECT_EQ(this->N,this->array->size());
};

TYPED_TEST_P(ArrayTests,Assign){
  // Inside a test, refer to TypeParam to get the type parameter.
  //TypeParam n = 0;
  //...
  ASSERT_NE(nullptr,this->array);
  this->array->resize(this->N);
  EXPECT_EQ(this->N,this->array->size());
  for(auto k=0;k<(this->array->size());k++){
    this->array->at(k)=((typename TypeParam::Dtype)(M_PI*((double)(k+1))));
    EXPECT_NE(0.0,this->array->at(k));
  };
  for(auto k=0;k<(this->array->size());k++){
    EXPECT_EQ(((typename TypeParam::Dtype)(M_PI*(double)(k+1))),this->array->at(k));
  };
};

TYPED_TEST_P(ArrayTests,Copy1){
  ASSERT_NE(nullptr,this->array);
  this->array->resize(this->N);
  EXPECT_EQ(this->N,this->array->size());
  for(auto k=0;k<(this->array->size());k++){
    this->array->at(k)=((typename TypeParam::Dtype)(M_PI*((double)(k+1))));
    EXPECT_NE(0.0,this->array->at(k));
  };
  auto array_copy=this->factory->makeManagedArray();
  array_copy->copy(this->array);
  for(auto k=0;k<(this->array->size());k++){
    EXPECT_EQ(this->array->at(k),array_copy->at(k));
  };
};

TYPED_TEST_P(ArrayTests,Copy2){
  ASSERT_NE(nullptr,this->array);
  this->array->resize(this->N);
  EXPECT_EQ(this->N,this->array->size());
  for(auto k=0;k<(this->array->size());k++){
    this->array->at(k)=((typename TypeParam::Dtype)(M_PI*((double)(k+1))));
    EXPECT_NE(0.0,this->array->at(k));
  };
  auto array_copy=this->factory->newArray();
  array_copy->copy(this->array);
  for(auto k=0;k<(this->array->size());k++){
    EXPECT_EQ(this->array->at(k),array_copy->at(k));
  };
  this->factory->deleteArray(array_copy);
};

TYPED_TEST_P(ArrayTests,Copy3){
  ASSERT_NE(nullptr,this->array);
  std::vector<typename TypeParam::Dtype> vec;
  vec.resize(this->N);
  for(auto k=0;k<vec.size();k++){
    vec.at(k)=((typename TypeParam::Dtype)(M_PI*((double)(k+1))));
    EXPECT_NE(0.0,vec.at(k));
  };
  PanNDE::CArray<typename TypeParam::Dtype> carray;
  carray.data=vec.data();carray.length=vec.size();
  this->array->copy(carray);
  ASSERT_LT(0,this->array->size());
  for(auto k=0;k<this->array->size();k++){
    EXPECT_NE(0.0,carray.data[k]);
    EXPECT_NE(0.0,this->array->at(k));
    EXPECT_EQ(carray.data[k],this->array->at(k));
  };
};

TYPED_TEST_P(ArrayTests,CArrayAccess){
  // Inside a test, refer to TypeParam to get the type parameter.
  //TypeParam n = 0;
  //...
  ASSERT_NE(nullptr,this->array);
  this->array->resize(this->N);
  EXPECT_EQ(this->N,this->array->size());
  for(auto k=0;k<(this->array->size());k++){
    this->array->at(k)=((typename TypeParam::Dtype)(M_PI*((double)(k+1))));
    EXPECT_NE(0.0,this->array->at(k));
  };
  auto carray=this->array->getCArray();
  EXPECT_LT(0,carray.length);
  for(auto k=0;k<carray.length;k++){
    EXPECT_EQ(this->array->at(k),carray.data[k]);
    carray.data[k]=-carray.data[k];
    EXPECT_EQ(this->array->at(k),carray.data[k]);
    EXPECT_NE(0.0,carray.data[k]);
  };
};

TYPED_TEST_P(ArrayTests,CheckPointer){
  // Inside a test, refer to TypeParam to get the type parameter.
  //TypeParam n = 0;
  //...
  ASSERT_NE(nullptr,this->array);
  this->array->resize(this->N);
  EXPECT_EQ(this->N,this->array->size());
  for(auto k=0;k<(this->array->size());k++){
    this->array->at(k)=((typename TypeParam::Dtype)(M_PI*((double)(k+1))));
    EXPECT_NE(0.0,this->array->at(k));
  };
  for(auto k=0;k<(this->array->size());k++){
    EXPECT_EQ(this->array->at(k),this->array->data()[k]);
  };
};

TYPED_TEST_P(ArrayTests,DataBundleTests){
  this->bundle=this->factory->makeManagedDataBundle();
  ASSERT_NE(nullptr,this->bundle);
  auto value=(typename TypeParam::Dtype)(M_PI*((double)(10)));
  this->bundle->emplaceScalar("A",value);
  this->bundle->emplaceScalar("B",-value);
  ASSERT_NE(nullptr,this->array);
  this->array->resize(this->N);
  EXPECT_EQ(this->N,this->array->size());
  for(auto k=0;k<(this->array->size());k++){
    this->array->at(k)=((typename TypeParam::Dtype)(M_PI*((double)(k+1))));
    EXPECT_NE(0.0,this->array->at(k));
  };
  this->bundle->emplaceArray("Aa",this->array);

  EXPECT_EQ(this->array,this->bundle->array("Aa"));
  EXPECT_EQ(value,this->bundle->scalar("A"));
  EXPECT_EQ(-value,this->bundle->scalar("B"));
};

REGISTER_TYPED_TEST_SUITE_P(ArrayTests,
                            Nothing,Resize,Assign,Copy1,Copy2,Copy3,
                            CArrayAccess,CheckPointer,DataBundleTests);
