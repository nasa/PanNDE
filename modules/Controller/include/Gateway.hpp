/*! \headerfile Gateway.hpp "modules/Controller/include/Gateway.hpp"
* "Gateway.hpp" contains the class definition encapsulating the file input/output management
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

#include <string>
#include <memory>

#include "Mesh.hpp"
#include "Field.hpp"
#include "Array.hpp"
#include "Communicator.hpp"


namespace Controller {
  /*! \class Gateway Gateway.hpp "modules/Controller/include/Gateway.hpp"
  *
  * Defines the methods required to describe the data for file I/O 
  *
  */
  class Gateway {
    public:
      /*!
      * open file for read
      * \param filename the name of the file (using relative addr)
      */
      virtual void open(std::string filename) =0;
      /*!
      * close files
      */
      virtual void close(){};

      /*!
      * get array by keyname from file
      * \param keyname name of variable to be pulled from input file
      * \param maker array factory
      */
      virtual std::shared_ptr<PanNDE::Array<double>> getArray(std::string keyname,
                              std::shared_ptr<PanNDE::ArrayFactory<double>> maker) =0;
      /*!
      * get value by keyname from file
      * \param keyname name of variable to be pulled from input file
      */
      virtual double getValue(std::string keyname) =0;

      /*!
      * get mesh from file
      * \param maker mesh factory
      * \param host_id root process to read mesh from file
      */
      virtual std::shared_ptr<PanNDE::Mesh> getMesh(
                              std::shared_ptr<PanNDE::MeshFactory> maker,int host_id=0) =0;
      /*!
      * get field by keyname from file
      * \param maker field factory
      * \param host_id root process to read mesh from file
      */
      virtual std::shared_ptr<PanNDE::Field> getField(std::string keyname,
                              std::shared_ptr<PanNDE::FieldFactory> maker,int host_id=0) =0;

      /*!
      * write solution variables to file
      * \param filename_base base name for file write
      * \param solution bundled solution data
      * \param write_index time index to be written
      */
      virtual void writeSolution(std::string filename_base,
                                 std::shared_ptr<PanNDE::FieldBundle> solution,
                                 int write_index) =0;
      /*!
      * write solution variables to file
      * \param filename_base base name for file write
      * \param solution bundled solution data
      * \param meta additional metadata to be written to file
      * \param write_index time index to be written
      */
      virtual void writeSolution(std::string filename_base,
                                 std::shared_ptr<PanNDE::FieldBundle> solution,
                                 std::shared_ptr<PanNDE::DataBundle<double>> meta,
                                 int write_index) =0;
      /*!
      * write solution variables to file
      * \param filename_base base name for file write
      * \param solution bundled solution data
      * \param meta additional metadata to be written to file
      */
      virtual void writeSolution(std::string filename_base,
                         std::shared_ptr<PanNDE::FieldBundle> solution,
                         std::shared_ptr<PanNDE::DataBundle<double>> meta) =0;
      /*!
      * write solution variables to file
      * \param filename_base base name for file write
      * \param solution bundled solution data
      */
      virtual void writeSolution(std::string filename_base,
                         std::shared_ptr<PanNDE::FieldBundle> solution) =0;
  };
};