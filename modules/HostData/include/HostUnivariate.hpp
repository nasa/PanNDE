/*! \headerfile HostUnivariate.hpp "modules/HostData/include/HostUnivariate.hpp"
* "HostUnivariate.hpp" contains a class implementation for representing and evaluating
* univariate functions (functions of a single variable). It supports evaluation at arbitrary
* points through linear interpolation between known data points.
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

#include <memory>
#include <map>
#include <stdexcept>

#include "Univariate.hpp"
#include "Array.hpp"

namespace HostData {
  /*! \class HostUnivariate HostUnivariate.hpp "modules/HostData/include/HostUnivariate.hpp"
  * 
  * Implements a univariate function defined by discrete data points.
  * This class represents a function of a single variable (typically time) through
  * a set of known (x,y) data points. It provides linear interpolation between points
  * and constant extrapolation beyond the defined range. The implementation uses an
  * ordered map for efficient lookup and caches the last evaluation for performance.
  *
  */
  class HostUnivariate : public PanNDE::Univariate {
    public:
      /*!
      * Constructs a univariate function from arrays of independent and dependent values.
      * \param time_vals std::shared_ptr<PanNDE::Array<double>> Array of x-values (independent variable)
      * \param func_vals std::shared_ptr<PanNDE::Array<double>> Array of y-values (dependent variable)
      * \throw std::runtime_error If the arrays have different lengths
      */
      HostUnivariate(std::shared_ptr<PanNDE::Array<double>> time_vals,
                     std::shared_ptr<PanNDE::Array<double>> func_vals){
        if(time_vals->size()!=func_vals->size()){throw std::runtime_error("mismatched array lengths");};
        for(int k=0;k<time_vals->size();k++){
          tseries.emplace(time_vals->at(k),func_vals->at(k));
        };
      };

      /*!
      * Creates a shared pointer to a new HostUnivariate instance.
      * \param time_vals std::shared_ptr<PanNDE::Array<double>> Array of x-values (independent variable)
      * \param func_vals std::shared_ptr<PanNDE::Array<double>> Array of y-values (dependent variable)
      * \return std::shared_ptr<HostData::HostUnivariate> Shared pointer to the new univariate function
      * \throw std::runtime_error If the arrays have different lengths
      */
      static
      std::shared_ptr<HostData::HostUnivariate> makeShared(std::shared_ptr<PanNDE::Array<double>> time_vals,
                                                           std::shared_ptr<PanNDE::Array<double>> func_vals){
        auto uv=std::make_shared<HostData::HostUnivariate>(HostData::HostUnivariate(time_vals,func_vals));
        return std::move(uv);
      };

      /*!
      * Evaluates the function at a specified point.
      * Uses linear interpolation between known points. For values beyond the
      * defined range, returns the function value at the closest endpoint.
      * \param t double The point at which to evaluate the function
      * \return double The interpolated function value
      *
      * \note Implements a caching strategy for performance when the same value
      *       is evaluated repeatedly
      */
      double evalAt(double t) override {
        if(t==last_eval_time){
          return last_eval;
        };
        if(t<(tseries.rbegin()->first)){
          last_eval_time=t;
          last_eval=interpolate(t);
          return last_eval;
        };
        last_eval_time=tseries.rbegin()->first;
        last_eval=tseries.rbegin()->second;
        return last_eval;
      };

    private:
      /*!
      * Performs linear interpolation between known data points.
      * \param t double The value at which to interpolate
      * \return double The interpolated function value
      *
      * \note For t values less than the smallest defined point, this will return
      *       an extrapolated value based on the first segment.
      */
      double interpolate(double t){
        for(auto it=tseries.begin();it!=tseries.end();it++){
          if(t<(it->first)){
            auto prev=std::prev(it);
            double Dt=(it->first)-(prev->first);
            double Dv=(it->second)-(prev->second);
            return (Dv/Dt)*(t-(prev->first))+(prev->second);
          };
        };
        return 0.;
      };

      //! Last computed function value (for caching)
      double last_eval=-1;
      
      //! Last evaluation point (for caching)
      double last_eval_time=-1;
      
      //! Map of known (x,y) function points, automatically sorted by x value
      std::map<double,double> tseries;
  };
};