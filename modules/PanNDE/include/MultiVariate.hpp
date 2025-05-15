/*! \headerfile MultiVariate.hpp "modules/PanNDE/include/MultiVariate.hpp"
* "MultiVariate.hpp" defines the interface for scalar functions of multiple variables
* used throughout the PanNDE framework to represent spatial and temporal distributions
* of physical quantities.
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

#include <memory>

#include "Array.hpp"

namespace PanNDE {
  /*! \class MultiVariate MultiVariate.hpp "modules/PanNDE/include/MultiVariate.hpp"
  *
  * Defines an interface for representing scalar functions of multiple variables.
  *
  * The MultiVariate class provides a generic way to represent mathematical
  * functions that take multiple input arguments and produce a single scalar output.
  * These functions are essential for defining spatial and temporal distributions
  * in physical simulations, including:
  *
  * - Source terms (e.g., heat sources, force distributions)
  * - Boundary conditions (e.g., pressure profiles, displacement fields)
  * - Material property distributions (e.g., spatially varying conductivity)
  * - Initial conditions (e.g., temperature profiles, velocity distributions)
  * - Analytical solutions for validation
  *
  * Vector-valued functions can be composed by using multiple MultiVariate
  * objects, one for each component.
  *
  */
  class MultiVariate {
    public:
      /*!
      * Evaluates the function at a specified point in the input space.
      *
      * This version accepts input arguments as a shared array.
      *
      * \param arguments std::shared_ptr<PanNDE::Array<double>> Array of input values
      * \return double The scalar function value at the specified point
      */
      virtual double evalAt(std::shared_ptr<PanNDE::Array<double>> arguments) =0;

      /*!
      * Evaluates the function at a specified point in the input space.
      *
      * This version accepts input arguments as a C-style array.
      *
      * \param arguments double[] Array of input values
      * \return double The scalar function value at the specified point
      */
      virtual double evalAt(double arguments[]) =0;

      /*!
      * Checks if a point is within the domain of the function.
      *
      * Implementations can use this to define restricted domains where the
      * function is valid. The default implementation always returns true,
      * indicating an unrestricted domain.
      *
      * \param arguments std::shared_ptr<PanNDE::Array<double>> Array of input values to check
      * \return bool True if the point is within the function domain
      */
      virtual bool checkAt(std::shared_ptr<PanNDE::Array<double>> arguments){return true;};

      /*!
      * Checks if a point is within the domain of the function.
      *
      * C-style array version of the domain check function.
      *
      * \param arguments double[] Array of input values to check
      * \return bool True if the point is within the function domain
      */
      virtual bool checkAt(double arguments[]){return true;};

      /*!
      * Gets the number of independent variables (input arguments) required.
      *
      * For example, a function of spatial coordinates and time f(x,y,z,t)
      * would return 4.
      *
      * \return int The number of input arguments
      */
      virtual int argumentCount() =0;

      /*!
      * Provides a sample input for a specific argument position.
      *
      * This can be used to obtain representative input values for testing,
      * initialization, or determining the expected range of each argument.
      *
      * \param arg_index int The index of the argument to sample (0-based)
      * \return std::shared_ptr<PanNDE::Array<double>> An array of sample values, or nullptr if not implemented
      */
      virtual std::shared_ptr<PanNDE::Array<double>> sampleArgument(int arg_index) =0;
  };
};