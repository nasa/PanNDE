//

#include <memory>
#include <cmath>

#include "Array.hpp"
#include "Univariate.hpp"
#include "MultiVariate.hpp"

#include "HostNormalXDZ.hpp"

#include "HostArray.hpp"
#include "HostUnivariate.hpp"

#include "gtest/gtest.h"


class HostNormalXDZTests : public ::testing::Test {
  protected:
    void SetUp() override {
      time=makeShared<HostData::HostArray<double>>();
      response=makeShared<HostData::HostArray<double>>();
      int Nt=1001;
      time->resize(Nt);
      response->resize(Nt);
      double dt=0.001;
      for(int k=0;k<Nt;k++){
        time->at(k)=k*dt;
        response->at(k)=0.5*(1.-cos(2.*M_PI*double(k)/double(Nt-1)))*(sin(2.*M_PI*5.*double(k)/double(Nt-1)));
      };
      signal=std::make_shared<HostData::HostUnivariate>(time,response);
    };
    void TearDown() override {};

    template<class T> 
    std::shared_ptr<T> makeShared(){
      return std::make_shared<T>(T());
    };

    double center[3]={1.0,-1.0,0.5};
    double ds=0.01;
    double radius=0.5;
    std::shared_ptr<HostData::HostUnivariate> signal;
    std::shared_ptr<HostData::HostArray<double>> time;
    std::shared_ptr<HostData::HostArray<double>> response;
};


TEST_F(HostNormalXDZTests,Nothing){};

TEST_F(HostNormalXDZTests,Create){
  auto XD=std::make_shared<Stubs::HostNormalXDZ>(Stubs::HostNormalXDZ(center,radius,ds,signal));
};
TEST_F(HostNormalXDZTests,Create2){
  auto XD=Stubs::HostNormalXDZ::makeShared(center,radius,ds,signal);
};

TEST_F(HostNormalXDZTests,CheckVals){
  auto XD=std::make_shared<Stubs::HostNormalXDZ>(Stubs::HostNormalXDZ(center,radius,ds,signal));
  double args[4]={1.,-1.,0.5,1.};
  EXPECT_FLOAT_EQ(XD->evalAt(args),0.0);
  args[0]=1.;args[1]=-1.;args[2]=0.5;args[3]=10.;
  EXPECT_FLOAT_EQ(XD->evalAt(args),0.0);
  args[0]=1.5;args[1]=-1.;args[2]=0.5;args[3]=0.25;
  EXPECT_FLOAT_EQ(XD->evalAt(args),0.0);
  args[0]=1.;args[1]=-1.;args[2]=0.5;args[3]=0.25;
  EXPECT_GT(XD->evalAt(args),0.0);
};
TEST_F(HostNormalXDZTests,CheckVals2){
  auto XD=Stubs::HostNormalXDZ::makeShared(center,radius,ds,signal);
  double args[4]={1.,-1.,0.5,1.};
  EXPECT_FLOAT_EQ(XD->evalAt(args),0.0);
  args[0]=1.;args[1]=-1.;args[2]=0.5;args[3]=10.;
  EXPECT_FLOAT_EQ(XD->evalAt(args),0.0);
  args[0]=1.5;args[1]=-1.;args[2]=0.5;args[3]=0.25;
  EXPECT_FLOAT_EQ(XD->evalAt(args),0.0);
  args[0]=1.;args[1]=-1.;args[2]=0.5;args[3]=0.25;
  EXPECT_GT(XD->evalAt(args),0.0);
};
