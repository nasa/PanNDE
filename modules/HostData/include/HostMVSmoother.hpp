/*! \headerfile HostMVSmoother.hpp "modules/HostData/include/HostMVSmoother.hpp"
* "HostMVSmoother.hpp" contains the class implementation for a multivariate smoothing interpolator.
* This class provides kernel-based smoothing interpolation for scattered data points in 
* N-dimensional space. It uses a cosine window function and weighted distance metrics.
*
* Note: This implementation prioritizes flexibility over performance and should be considered
* a beta product with potential performance limitations for large datasets.
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
  /*! \class HostMVSmoother HostMVSmoother.hpp "modules/HostData/include/HostMVSmoother.hpp"
  * 
  * Implements a kernel-based multivariate smoothing interpolator.
  * This class provides scattered data interpolation in N-dimensional space using
  * a cosine window function and dimension-specific weighting. It computes output
  * values as weighted averages of known sample points, with closer samples having
  * more influence on the result.
  * 
  * \warning This implementation is computationally intensive for large datasets
  * and should be considered a beta product with potential performance limitations.
  *
  */
  class HostMVSmoother : public PanNDE::MultiVariate {
    public:
      /*!
      * Constructs a smoother with specified dimension weights but no samples.
      * \param dim_weights std::shared_ptr<PanNDE::Array<double>> Weights for each dimension
      *        (higher weights increase the importance of that dimension in distance calculations)
      */
      HostMVSmoother(std::shared_ptr<PanNDE::Array<double>> dim_weights){
        Nargs=dim_weights->size();
        weights=makeArray();
        weights->copy(dim_weights);
      };
      
      /*!
      * Constructs a smoother with uniform weights across all dimensions but no samples.
      * \param Narguments int Number of input dimensions
      * \param weight double Weight value to use for all dimensions
      */
      HostMVSmoother(int Narguments,double weight){
        Nargs=Narguments;
        weights=makeArray();
        weights->resize(Nargs);
        for(int k=0;k<Nargs;k++){weights->at(k)=weight;};
      };
      
      /*!
      * Constructs a smoother with specified dimension weights and sample points.
      * \param dim_weights std::shared_ptr<PanNDE::Array<double>> Weights for each dimension
      * \param samps std::vector<std::shared_ptr<PanNDE::Array<double>>>& Sample points, each
      *        containing N input values followed by 1 output value
      * \throw std::runtime_error If any sample has incorrect dimensions
      */
      HostMVSmoother(std::shared_ptr<PanNDE::Array<double>> dim_weights,
                     std::vector<std::shared_ptr<PanNDE::Array<double>>>& samps){
        Nargs=dim_weights->size();
        weights=makeArray();
        weights->copy(dim_weights);
        loadSamples(samps);
      };
      
      /*!
      * Constructs a smoother with uniform weights and sample points.
      * \param Narguments int Number of input dimensions
      * \param weight double Weight value to use for all dimensions
      * \param samps std::vector<std::shared_ptr<PanNDE::Array<double>>>& Sample points, each
      *        containing N input values followed by 1 output value
      * \throw std::runtime_error If any sample has incorrect dimensions
      */
      HostMVSmoother(int Narguments,double weight,
                     std::vector<std::shared_ptr<PanNDE::Array<double>>>& samps){
        Nargs=Narguments;
        weights=makeArray();
        weights->resize(Nargs);
        for(int k=0;k<Nargs;k++){weights->at(k)=weight;};
        loadSamples(samps);
      };

      /*!
      * Evaluates the smoothing function at the specified point.
      * Computes a weighted average of sample values, where the weight for each sample
      * depends on its distance from the evaluation point.
      * 
      * \param arguments std::shared_ptr<PanNDE::Array<double>> Coordinates of the point to evaluate
      * \return double Interpolated function value at the point
      * \throw std::runtime_error If the arguments array has incorrect dimensions
      */
      double evalAt(std::shared_ptr<PanNDE::Array<double>> arguments)override{
        if(arguments->size()!=Nargs){throw std::runtime_error("invalid arguments provided");};
        return kernel(arguments->data());
      };
      
      /*!
      * Evaluates the smoothing function at the specified point using raw pointer access.
      * \param arguments double* Pointer to array of coordinates
      * \return double Interpolated function value at the point
      * \warning No size checking is performed; caller must ensure array has correct dimensions
      */
      double evalAt(double* arguments)override{
        return kernel(arguments);
      };

      /*!
      * Gets the number of input dimensions for this function.
      * \return int The number of dimensions (arguments)
      */
      int argumentCount()override{
        return Nargs;
      };
      
      /*!
      * Retrieves a specific sample point.
      * \param arg_index int Index of the sample to retrieve
      * \return std::shared_ptr<PanNDE::Array<double>> The sample point (N input values + 1 output value)
      * \throw std::out_of_range If index is out of bounds
      */
      std::shared_ptr<PanNDE::Array<double>> sampleArgument(int arg_index)override{
        return samples.at(arg_index);
      };
      
    private:
      /*!
      * Creates a new empty array.
      * \return std::shared_ptr<PanNDE::Array<double>> New empty array
      */
      std::shared_ptr<PanNDE::Array<double>> makeArray(){
        return std::move(std::make_shared<HostData::HostArray<double>>(HostData::HostArray<double>()));
      };
      
      /*!
      * Computes the weighted squared distance between a point and a sample.
      * \param idx int Index of the sample
      * \param args double* Coordinates of the point
      * \return double Weighted squared distance
      */
      inline double sqrDist(int idx, double* args){
        double retval=0.;
        for(int k=0;k<Nargs;k++){
          retval+=sqr(weights->at(k)*(args[k]-samples.at(idx)->at(k)));
        };
        return retval;
      };
      
      /*!
      * Computes the square of a value.
      * \param x double Input value
      * \return double Square of the input
      */
      inline double sqr(double x){return x*x;};
      
      /*!
      * Window function that determines the influence of samples based on distance.
      * Uses a cosine-based function that smoothly decreases from 1 to 0 as distance
      * increases, reaching 0 at a distance of 1.
      * 
      * \param dist2 double Squared distance value
      * \return double Window function value (between 0 and 1)
      */
      inline double window(double dist2){
        if(dist2<1.){return (0.5*(1.0+cos(M_PI*sqrt(dist2))));};
        return 0.;
      };
      
      /*!
      * Kernel function that computes the smoothed value at a point.
      * Calculates a weighted average of all sample values based on their
      * distance from the evaluation point.
      * 
      * \param args double* Coordinates of the evaluation point
      * \return double Smoothed function value
      */
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
      
      /*!
      * Loads a collection of sample points.
      * \param samps std::vector<std::shared_ptr<PanNDE::Array<double>>>& Vector of sample points
      * \throw std::runtime_error If any sample has incorrect dimensions
      */
      void loadSamples(std::vector<std::shared_ptr<PanNDE::Array<double>>>& samps){
        samples.resize(0);
        samples.reserve(samps.size());
        for(int k=0;k<samps.size();k++){loadSample(samps.at(k));};
      };
      
      /*!
      * Loads a single sample point.
      * \param samp std::shared_ptr<PanNDE::Array<double>>& The sample to load
      * \throw std::runtime_error If sample has incorrect dimensions
      */
      void loadSample(std::shared_ptr<PanNDE::Array<double>>& samp){
        if(samp->size()!=(Nargs+1)){throw std::runtime_error("invalid sample provided");};
        samples.push_back(makeArray());
        samples.back()->copy(samp);
      };

      //! Number of input dimensions
      int Nargs=0;
      
      //! Weighting factors for each dimension
      std::shared_ptr<PanNDE::Array<double>> weights;
      
      //! Collection of sample points, each containing N input coordinates plus 1 output value
      std::vector<std::shared_ptr<PanNDE::Array<double>>> samples;
  };
};