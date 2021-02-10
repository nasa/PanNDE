// Not the best testing, but enough to check the HostMVSmoother for prototyping
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
#include <vector>

#include "Array.hpp"
#include "MultiVariate.hpp"

#include "gtest/gtest.h"

using testing::Types;


template<class FT>
std::shared_ptr<FT> makeShared(){
  return std::move(std::make_shared<FT>(FT()));
};

template <class I,template<typename T>class AF>
struct MultiVariateTypeDefinitions
{
  typedef I Impl;
  typedef AF<double> DArrayMfg;
};

template<class T>
class MultiVariateTests : public ::testing::Test {
  protected:
    void SetUp() override {};
    void TearDown() override {};

    std::shared_ptr<PanNDE::MultiVariate> smoother=nullptr;
    double dx=0.1;

    void make1DSmoother(double weight,std::vector<double> values){
      auto fx=produce1DSamples(values);
      smoother=std::make_shared<typename T::Impl>(typename T::Impl(1,weight,fx));
    };
    std::vector<std::shared_ptr<PanNDE::Array<double>>> produce1DSamples(
                                      std::vector<double> values){
      int Nsamples=values.size();
      std::vector<std::shared_ptr<PanNDE::Array<double>>> samples;
      samples.resize(Nsamples);
      auto dbl_maker=makeShared<typename T::DArrayMfg>();
      for(int kx=0;kx<Nsamples;kx++){
        samples.at(kx)=dbl_maker->makeManagedArray();
        samples.at(kx)->resize(2);
        samples.at(kx)->at(0)=dx*double(kx+20);
        samples.at(kx)->at(1)=values.at(kx);
      };
      return samples;
    };
    void make2DSmoother(double weight){
      auto fx=produce2DSamples();
      smoother=std::make_shared<typename T::Impl>(typename T::Impl(2,weight,fx));
    };
    std::vector<std::shared_ptr<PanNDE::Array<double>>> produce2DSamples(){
      int Nside=5;
      std::vector<std::shared_ptr<PanNDE::Array<double>>> samples;
      samples.resize(Nside*Nside);
      auto dbl_maker=makeShared<typename T::DArrayMfg>();
      double x0=10.*dx;
      double y0=10.*dx;
      for(int kx=0;kx<Nside;kx++){
        for(int ky=0;ky<Nside;ky++){
          int idx=kx*Nside+ky;
          samples.at(idx)=dbl_maker->makeManagedArray();
          samples.at(idx)->resize(3);
          double x=dx*double(kx+7);
          double y=dx*double(ky+7);
          samples.at(idx)->at(0)=x;
          samples.at(idx)->at(1)=y;
          samples.at(idx)->at(2)=exp(-((x-x0)*(x-x0)+(y-y0)*(y-y0)));
        };
      };
      return samples;
    };
};

TYPED_TEST_SUITE_P(MultiVariateTests);

TYPED_TEST_P(MultiVariateTests,Nothing){};

TYPED_TEST_P(MultiVariateTests,MakeSmoother1D){
  int Ns=25;
  std::vector<double> values;values.resize(Ns);
  for(int k=0;k<Ns;k++){
    values.at(k)=1.-cos(2.0*M_PI*(double(k)/(double(Ns-1))));
  };
  this->make1DSmoother(10./(this->dx),values);
  auto dbl_maker=makeShared<typename TypeParam::DArrayMfg>();
  auto arg=dbl_maker->makeManagedArray();
  arg->resize(1);
  for(int k=0;k<3*Ns;k++){
    arg->at(0)=double(k)*(this->dx);
    double tv=this->smoother->evalAt(arg);
    EXPECT_LE(0.,tv);
  };
};
TYPED_TEST_P(MultiVariateTests,MakeSmoother2D){

  this->make2DSmoother(10./(this->dx));

  auto dbl_maker=makeShared<typename TypeParam::DArrayMfg>();
  auto arg=dbl_maker->makeManagedArray();
  arg->resize(2);
  for(int kx=0;kx<25;kx++){
    for(int ky=0;ky<25;ky++){
      arg->at(0)=double(kx)*(this->dx);
      arg->at(1)=double(ky)*(this->dx);
      double tv=this->smoother->evalAt(arg);
      EXPECT_LE(0.,tv);
      //printf("%g ",tv);
    };
    //printf("\n");
  };
};

REGISTER_TYPED_TEST_SUITE_P(MultiVariateTests,
                            Nothing,MakeSmoother1D,MakeSmoother2D);
