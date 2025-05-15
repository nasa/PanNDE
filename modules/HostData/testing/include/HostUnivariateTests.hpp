//

#include <memory>
#include <vector>

#include "Array.hpp"
#include "HostArray.hpp"
#include "Univariate.hpp"
#include "HostUnivariate.hpp"

#include "gtest/gtest.h"


class HostUnivariateTests : public ::testing::Test {
  protected:
    void SetUp() override {
      int N=101;
      time=makeShared<HostData::HostArray<double>>();
      response=makeShared<HostData::HostArray<double>>();
      time->resize(N);
      response->resize(N);
      for(int k=0;k<N;k++){
        time->at(k)=0.01*k;
        response->at(k)=-10.*k;
      };
    };
    void TearDown() override {};

    template<class T> 
    std::shared_ptr<T> makeShared(){
      return std::make_shared<T>(T());
    };

    std::shared_ptr<HostData::HostUnivariate> MakeUV(){
      return std::make_shared<HostData::HostUnivariate>(HostData::HostUnivariate(time,response));    
    };
    std::shared_ptr<HostData::HostArray<double>> time;
    std::shared_ptr<HostData::HostArray<double>> response;
};


TEST_F(HostUnivariateTests,Nothing){};

TEST_F(HostUnivariateTests,Create){
  auto tseries=MakeUV();
};

TEST_F(HostUnivariateTests,Create2){
  auto tseries=HostData::HostUnivariate::makeShared(time,response);
};

TEST_F(HostUnivariateTests,CheckLinear){
  auto tseries=MakeUV();
  EXPECT_FLOAT_EQ(tseries->evalAt(0.01),-10.);
  EXPECT_FLOAT_EQ(tseries->evalAt(0.001),-1.);
  EXPECT_FLOAT_EQ(tseries->evalAt(0.0001),-0.1);
  EXPECT_FLOAT_EQ(tseries->evalAt(0.015),-15.);
  EXPECT_FLOAT_EQ(tseries->evalAt(1.),-1000.);
  EXPECT_FLOAT_EQ(tseries->evalAt(1.1),-1000.);
  EXPECT_FLOAT_EQ(tseries->evalAt(11.),-1000.);
};

TEST_F(HostUnivariateTests,CheckLinear2){
  std::shared_ptr<PanNDE::Univariate> tseries=HostData::HostUnivariate::makeShared(time,response);
  EXPECT_FLOAT_EQ(tseries->evalAt(0.01),-10.);
  EXPECT_FLOAT_EQ(tseries->evalAt(0.001),-1.);
  EXPECT_FLOAT_EQ(tseries->evalAt(0.0001),-0.1);
  EXPECT_FLOAT_EQ(tseries->evalAt(0.015),-15.);
  EXPECT_FLOAT_EQ(tseries->evalAt(1.),-1000.);
  EXPECT_FLOAT_EQ(tseries->evalAt(1.1),-1000.);
  EXPECT_FLOAT_EQ(tseries->evalAt(11.),-1000.);
};

