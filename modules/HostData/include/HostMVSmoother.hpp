/*! \headerfile HostMVSmoother.hpp "modules/HostData/include/HostMVSmoother.hpp"
* "HostMVSmoother.hpp" contains the class definition for creating a generic scalar multivariate object.
* Performance is bad, currently here for development reasons, use as beta product if one must do so.
*
*/

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

#include <cmath>
#include <vector>
#include <memory>

#include "MultiVariate.hpp"
#include "Array.hpp"

#include "HostArray.hpp"

namespace HostData {
  class HostMVSmoother : public PanNDE::MultiVariate {
    public:
      HostMVSmoother(std::shared_ptr<PanNDE::Array<double>> dim_weights){
        Nargs=dim_weights->size();
        weights=makeArray();
        weights->copy(dim_weights);
      };
      HostMVSmoother(int Narguments,double weight){
        Nargs=Narguments;
        weights=makeArray();
        weights->resize(Nargs);
        for(int k=0;k<Nargs;k++){weights->at(k)=weight;};
      };
      HostMVSmoother(std::shared_ptr<PanNDE::Array<double>> dim_weights,
                     std::vector<std::shared_ptr<PanNDE::Array<double>>>& samps){
        Nargs=dim_weights->size();
        weights=makeArray();
        weights->copy(dim_weights);
        loadSamples(samps);
      };
      HostMVSmoother(int Narguments,double weight,
                     std::vector<std::shared_ptr<PanNDE::Array<double>>>& samps){
        Nargs=Narguments;
        weights=makeArray();
        weights->resize(Nargs);
        for(int k=0;k<Nargs;k++){weights->at(k)=weight;};
        loadSamples(samps);
      };

      double evalAt(std::shared_ptr<PanNDE::Array<double>> arguments)override{
        if(arguments->size()!=Nargs){throw std::runtime_error("invalid arguments provided");};
        return kernel(arguments->data());
      };
      //not safe pointer
      double evalAt(double* arguments)override{
        return kernel(arguments);
      };

      //necessary utilities
      int argumentCount()override{
        return Nargs;
      };
      std::shared_ptr<PanNDE::Array<double>> sampleArgument(int arg_index)override{
        return samples.at(arg_index);
      };
    private:
      std::shared_ptr<PanNDE::Array<double>> makeArray(){
        return std::move(std::make_shared<HostData::HostArray<double>>(HostData::HostArray<double>()));
      };
      inline double sqrDist(int idx, double* args){
        double retval=0.;
        for(int k=0;k<Nargs;k++){
          retval+=sqr(weights->at(k)*(args[k]-samples.at(idx)->at(k)));
        };
        return retval;
      };
      inline double sqr(double x){return x*x;};
      inline double window(double dist2){
        if(dist2<1.){return (0.5*(1.0+cos(M_PI*sqrt(dist2))));};
        return 0.;
      };
      inline double kernel(double* args){
        double sumval=0.;
        double retval=0.;
        for(int k=0;k<samples.size();k++){
          double winval=window(sqrDist(k,args));
          sumval+=winval;
          retval+=winval*samples.at(k)->at(Nargs);
          //printf("winval(%g-%g)=%g\n",args[0],samples.at(k)->at(0),winval);
        };
        //printf("at %g rv %g sv %g\n",args[0],retval,sumval);
        retval=retval/sumval;
        if(sumval>0.){return retval;};
        return 0.;
      };
      void loadSamples(std::vector<std::shared_ptr<PanNDE::Array<double>>>& samps){
        samples.resize(0);
        samples.reserve(samps.size());
        for(int k=0;k<samps.size();k++){loadSample(samps.at(k));};
      };
      void loadSample(std::shared_ptr<PanNDE::Array<double>>& samp){
        if(samp->size()!=(Nargs+1)){throw std::runtime_error("invalid sample provided");};
        samples.push_back(makeArray());
        samples.back()->copy(samp);
      };

      int Nargs=0;
      std::shared_ptr<PanNDE::Array<double>> weights;//must be of length of Nargs
      std::vector<std::shared_ptr<PanNDE::Array<double>>> samples;//Array contents must be of length Nargs+1
  };
};