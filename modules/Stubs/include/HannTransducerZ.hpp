/*! \headerfile HannTransducerZ.hpp "modules/Stubs/include/HannTransducerZ.hpp"
* "HannTransducerZ.hpp" contains a class that implements a Hann-windowed ultrasonic transducer
* operating along the Z-axis for testing and simulation purposes.
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
#include <cstdlib>
#include <memory>

#include "MultiVariate.hpp"
#include "Array.hpp"

namespace Stubs {
  /*! \class HannTransducerZ HannTransducerZ.hpp "modules/Stubs/include/HannTransducerZ.hpp"
  *
  * Implements a Hann-windowed ultrasonic transducer operating along the Z-axis.
  * This class provides a mathematical representation of an ultrasonic transducer
  * that generates a sinusoidal signal modulated by a Hann window function in time,
  * with a spatial profile defined by a circular region in the XY plane.
  *
  */
  class HannTransducerZ : public PanNDE::MultiVariate {
    public:
      /*!
      * Constructor for a Hann-windowed transducer.
      * \param center double[3] The center coordinates of the transducer
      * \param radius double The radius of the transducer
      * \param ds double The spatial discretization size
      */
      HannTransducerZ(double center[3],double radius,double ds){
        for(int k=0;k<3;k++){
          this->center[k]=center[k];
        };
        this->radius=radius;
        this->ds=ds;
      };

      /*!
      * Gets the number of arguments required for evaluation.
      * \return int Always returns 4 (x,y,z,t)
      */
      int argumentCount()override{return 4;};
      
      /*!
      * Sample argument method (not implemented).
      * \param arg_index int The argument index
      * \return std::shared_ptr<PanNDE::Array<double>> Always returns nullptr
      */
      std::shared_ptr<PanNDE::Array<double>> sampleArgument(int arg_index)override{return nullptr;};

      /*!
      * Evaluates the transducer function at a given point and time using an array.
      * \param arguments std::shared_ptr<PanNDE::Array<double>> Array containing (x,y,z,t)
      * \return double The value of the transducer function
      */
      double evalAt(std::shared_ptr<PanNDE::Array<double>> arguments) override {
        return (evalAt(arguments->data()));
      };
      
      /*!
      * Evaluates the transducer function at a given point and time.
      * \param arguments double[] Array containing (x,y,z,t)
      * \return double The value of the transducer function
      * 
      * Applies both spatial windowing (circular profile) and temporal windowing
      * (Hann-windowed sinusoid) to simulate the transducer output.
      */
      double evalAt(double arguments[])override{
        double t=arguments[3];
        double dist=distance(arguments,center);
        double scale=(t>tMax || t<0.0)?0.0:hannWindowedSin(t);
        double spatial_window=double((fabs(arguments[2]-center[2]))<=(0.75*ds));
        //double window=0.5*(cos(M_PI*dist/radius)+1.)*(dist<radius)*spatial_window*(dist>=0.);
        double window=spatial_window*(dist<radius)*(dist>=0.);
        return scale*window;
      };

      /*!
      * Configures the transducer parameters.
      * \param Ncycle int32_t Number of cycles in the pulse (default: 1)
      * \param frequency double Operating frequency in Hz (default: 200 kHz)
      * \param phaseOffset double Phase offset in radians (default: 0.0)
      */
      void configure(int32_t Ncycle=1,double frequency=200.0e3,double phaseOffset=0.0){
        this->Ncycle=Ncycle;
        this->frequency=frequency;
        this->tMax=(double(Ncycle))/frequency;
        this->phaseOffset=phaseOffset;
      };

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
      
      /*!
      * Calculates a Hann-windowed sinusoid at a given time.
      * \param t double The time
      * \return double The value of the Hann-windowed sinusoid
      */
      double hannWindowedSin(double t)
      {
        if(Ncycle==0){return 0.;};
        return 0.5*sin(2.0*M_PI*frequency*t+phaseOffset)*(1.0-cos(2.0*M_PI*t/tMax));
      };
      
      //! Number of cycles in the pulse
      int32_t Ncycle=1;
      //! Operating frequency in Hz
      double frequency=200.0e3;
      //! Maximum time duration of the pulse
      double tMax=(1.0/200.0e3);
      //! Phase offset in radians
      double phaseOffset=0.0;
      //! Center coordinates of the transducer
      double center[3];
      //! Spatial discretization size
      double ds;
      //! Radius of the transducer
      double radius;
  };
};