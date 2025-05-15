/*! \headerfile HostNormalXDZ.hpp "modules/Stubs/include/HostNormalXDZ.hpp"
* "HostNormalXDZ.hpp" contains a class that implements a circular transducer with normal
* excitation in the Z-direction, using an external time-domain signal.
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
#include <memory>

#include "MultiVariate.hpp"
#include "Array.hpp"
#include "Univariate.hpp"

namespace Stubs {
  /*! \class HostNormalXDZ HostNormalXDZ.hpp "modules/Stubs/include/HostNormalXDZ.hpp"
  *
  * Implements a circular transducer that applies a provided time-domain signal
  * over a defined circular region, with normal excitation in the Z-direction.
  * Uses an externally provided univariate time signal function, applying it with
  * a spatial window to create a realistic transducer excitation pattern.
  *
  */
  class HostNormalXDZ : public PanNDE::MultiVariate {
    public:
      /*!
      * Constructor for a normal-excitation transducer.
      * \param center double[3] The center coordinates of the transducer
      * \param radius double The radius of the transducer
      * \param ds double The spatial discretization size
      * \param signal std::shared_ptr<PanNDE::Univariate> Time-domain signal function
      */
      HostNormalXDZ(double center[3],double radius,double ds,
                    std::shared_ptr<PanNDE::Univariate> signal){
        tseries=signal;
        this->center[0]=center[0];this->center[1]=center[1];this->center[2]=center[2];
        this->radius=radius;
        this->ds=ds;
      };
      
      /*!
      * Creates a shared pointer to a HostNormalXDZ object.
      * \param center double[3] The center coordinates of the transducer
      * \param radius double The radius of the transducer
      * \param ds double The spatial discretization size
      * \param signal std::shared_ptr<PanNDE::Univariate> Time-domain signal function
      * \return std::shared_ptr<HostNormalXDZ> Shared pointer to a newly created object
      */
      static
      std::shared_ptr<HostNormalXDZ> makeShared(double center[3],double radius,double ds,
                                                std::shared_ptr<PanNDE::Univariate> signal){
        auto xd=std::make_shared<HostNormalXDZ>(HostNormalXDZ(center,radius,ds,signal));
        return std::move(xd);
      };

      /*!
      * Evaluates the transducer function at a given point and time using an array.
      * \param arguments std::shared_ptr<PanNDE::Array<double>> Array containing (x,y,z,t)
      * \return double The value of the transducer function
      */
      double evalAt(std::shared_ptr<PanNDE::Array<double>> arguments) override {
        //x=args[0],y=args[1],z=args[2],t=args[3]
        return evalAt(arguments->data());
      };

      /*!
      * Evaluates the transducer function at a given point and time.
      * \param arguments double[] Array containing (x,y,z,t)
      * \return double The value of the transducer function
      * 
      * Applies both spatial windowing (circular profile with z-constraint) and 
      * the provided time signal to simulate the transducer output.
      */
      double evalAt(double arguments[]) override {
        double dist=distance(arguments,center);
        double scale=tseries->evalAt(arguments[3]);
        double distz=fabs(arguments[2]-center[2]);
        double window=(distz<=0.75*ds)*(dist<radius)*(dist>=0.); 
        return scale*window;
      };

      /*!
      * Checks if a point is within the active region of the transducer.
      * \param arguments std::shared_ptr<PanNDE::Array<double>> Array containing (x,y,z,t)
      * \return bool True if the point is within the active region
      */
      bool checkAt(std::shared_ptr<PanNDE::Array<double>> arguments) override {
        //x=args[0],y=args[1],z=args[2],t=args[3]
        return checkAt(arguments->data());
      };

      /*!
      * Checks if a point is within the active region of the transducer.
      * \param arguments double[] Array containing (x,y,z,t)
      * \return bool True if the point is within the active region
      */
      bool checkAt(double arguments[]) override {
        double dist=distance(arguments,center);
        double distz=fabs(arguments[2]-center[2]);
        bool check=(distz<=0.75*ds)*(dist<radius)*(dist>=0.); 
        return check;
      };

      /*!
      * Gets the number of arguments required for evaluation.
      * \return int Always returns 4 (x,y,z,t)
      */
      int argumentCount() override {return 4;};
      
      /*!
      * Sample argument method (not implemented).
      * \param arg_index int The argument index
      * \return std::shared_ptr<PanNDE::Array<double>> Always returns nullptr
      */
      std::shared_ptr<PanNDE::Array<double>> sampleArgument(int arg_index)override{return nullptr;};

    private:
      /*!
      * Calculates the square of a value.
      * \param x double Input value
      * \return double The squared value
      */
      inline double sqr(double x){return x*x;};
      
      /*!
      * Calculates the Euclidean distance between two 3D points.
      * \param pti double[3] First point
      * \param ptr double[3] Second point
      * \return double The distance between the points
      */
      inline double distance(double pti[3],double ptr[3]){
        double res=0.0;for(int k=0;k<3;k++){res+=sqr(pti[k]-ptr[k]);};
        return sqrt(res);
      };
      
      //! Time-domain signal to apply
      std::shared_ptr<PanNDE::Univariate> tseries;
      //! Radius of the transducer
      double radius;
      //! Spatial discretization size
      double ds;
      //! Center coordinates of the transducer
      double center[3];
  };
};