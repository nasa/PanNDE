/*! \headerfile Model.hpp "modules/PanNDE/include/Model.hpp"
* "Model.hpp" defines the minimal interface for physics models and solvers in PanNDE.
* This intentionally sparse interface provides the essential execution and output methods
* that all models must implement, while allowing specific physics implementations to 
* define their own configuration and parameter interfaces as appropriate.
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

#include "Mesh.hpp"
#include "Field.hpp"

namespace PanNDE{
  /*! \class Model Model.hpp "modules/PanNDE/include/Model.hpp"
  *
  * Defines the core interface for computational physics models and numerical solvers.
  * 
  * This deliberately minimal interface provides only the essential methods required
  * for executing simulations and retrieving results. The sparse design accommodates
  * diverse physics domains (ultrasonic, thermal, electromagnetic, etc.) without 
  * imposing unnecessary constraints on implementation details.
  * 
  * Concrete model implementations are expected to:
  * 1. Define their own domain-specific configuration methods
  * 2. Handle appropriate internal state management
  * 3. Process boundary and initial conditions as needed
  * 4. Implement appropriate numerical algorithms
  * 
  * The interface follows the Open-Closed Principle - it can be extended for new
  * functionality without modifying the core interface itself.
  * 
  */
  class Model{
    public:
      /*!
      * Executes the model simulation for a specified number of steps.
      * 
      * This method advances the internal state of the simulation according to
      * the model's physics and numerical algorithms. It should handle time stepping,
      * convergence checks, and any required state updates.
      * 
      * \param nsteps int Number of model steps to execute (default: 1)
      * \return double A meaningful indicator of model progress, such as:
      *    - Convergence residual for implicit solvers
      *    - Error estimate for adaptive schemes
      *    - Time step size for variable time stepping methods
      *    - Iteration count for iterative solvers
      */
      virtual double solve(int nsteps=1) =0;
      
      /*!
      * Gets the current internal state variables of the simulation.
      * 
      * State variables represent the complete computational state of the model,
      * including primary field variables (e.g., displacement, temperature) and
      * any auxiliary variables needed to advance the solution (e.g., velocities,
      * intermediate values). They fully define the current state of the simulation.
      * 
      * \return std::shared_ptr<PanNDE::FieldBundle> Bundle containing all state fields
      */
      virtual std::shared_ptr<PanNDE::FieldBundle> getStates() =0;
  };
};