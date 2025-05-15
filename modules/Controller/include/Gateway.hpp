/*! \headerfile Gateway.hpp "modules/Controller/include/Gateway.hpp"
* "Gateway.hpp" contains the abstract interface for file input/output operations.
* It defines methods for reading and writing simulation data, meshes, fields, and
* tabular data from/to external files in various formats through a consistent interface.
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
  * Defines an abstract interface for file input/output operations.
  * 
  * The Gateway class provides a unified interface for reading and writing simulation
  * data across different file formats. It abstracts the details of specific file formats
  * and provides methods for accessing meshes, fields, arrays, and scalar values from files.
  * It also supports writing solution data and tabular data with associated metadata.
  * 
  * Concrete implementations of this interface handle specific file formats like VTK,
  * HDF5, or custom formats, allowing the rest of the code to remain agnostic to the
  * actual storage format used.
  *
  */
  class Gateway {
    public:
      /*!
      * Opens a file for reading.
      * \param filename std::string The name of the file (using relative address)
      */
      virtual void open(std::string filename) =0;
      
      /*!
      * Closes currently open file(s).
      * Default implementation does nothing; derived classes may override.
      */
      virtual void close(){};

      /*!
      * Retrieves an array by name from the open file.
      * \param keyname std::string Name of the variable to retrieve
      * \param maker std::shared_ptr<PanNDE::ArrayFactory<double>> Factory for creating the array
      * \return std::shared_ptr<PanNDE::Array<double>> The retrieved array data
      */
      virtual std::shared_ptr<PanNDE::Array<double>> getArray(std::string keyname,
                              std::shared_ptr<PanNDE::ArrayFactory<double>> maker) =0;
      
      /*!
      * Retrieves a scalar value by name from the open file.
      * \param keyname std::string Name of the variable to retrieve
      * \return double The value from the file
      */
      virtual double getValue(std::string keyname) =0;

      /*!
      * Retrieves a mesh from the open file.
      * \param maker std::shared_ptr<PanNDE::MeshFactory> Factory for creating the mesh
      * \param host_id int Process ID that should read the mesh (in parallel contexts, default: 0)
      * \return std::shared_ptr<PanNDE::Mesh> The retrieved mesh
      */
      virtual std::shared_ptr<PanNDE::Mesh> getMesh(
                              std::shared_ptr<PanNDE::MeshFactory> maker,int host_id=0) =0;
      
      /*!
      * Retrieves a field by name from the open file.
      * \param keyname std::string Name of the field to retrieve
      * \param maker std::shared_ptr<PanNDE::FieldFactory> Factory for creating the field
      * \param host_id int Process ID that should read the field (in parallel contexts, default: 0)
      * \return std::shared_ptr<PanNDE::Field> The retrieved field
      */
      virtual std::shared_ptr<PanNDE::Field> getField(std::string keyname,
                              std::shared_ptr<PanNDE::FieldFactory> maker,int host_id=0) =0;

      /*!
      * Writes solution data to file with a time index.
      * \param filename_base std::string Base name for the output file
      * \param solution std::shared_ptr<PanNDE::FieldBundle> The solution data to write
      * \param write_index int Time index to associate with this output
      */
      virtual void writeSolution(std::string filename_base,
                                 std::shared_ptr<PanNDE::FieldBundle> solution,
                                 int write_index) =0;
      
      /*!
      * Writes solution data and metadata to file with a time index.
      * \param filename_base std::string Base name for the output file
      * \param solution std::shared_ptr<PanNDE::FieldBundle> The solution data to write
      * \param meta std::shared_ptr<PanNDE::DataBundle<double>> Additional metadata to include
      * \param write_index int Time index to associate with this output
      */
      virtual void writeSolution(std::string filename_base,
                                 std::shared_ptr<PanNDE::FieldBundle> solution,
                                 std::shared_ptr<PanNDE::DataBundle<double>> meta,
                                 int write_index) =0;
      
      /*!
      * Writes solution data and metadata to file.
      * \param filename_base std::string Base name for the output file
      * \param solution std::shared_ptr<PanNDE::FieldBundle> The solution data to write
      * \param meta std::shared_ptr<PanNDE::DataBundle<double>> Additional metadata to include
      */
      virtual void writeSolution(std::string filename_base,
                         std::shared_ptr<PanNDE::FieldBundle> solution,
                         std::shared_ptr<PanNDE::DataBundle<double>> meta) =0;
      
      /*!
      * Writes solution data to file.
      * \param filename_base std::string Base name for the output file
      * \param solution std::shared_ptr<PanNDE::FieldBundle> The solution data to write
      */
      virtual void writeSolution(std::string filename_base,
                         std::shared_ptr<PanNDE::FieldBundle> solution) =0;

      /*!
      * Writes tabular data and metadata to file with a time index.
      * \param filename_base std::string Base name for the output file
      * \param data std::shared_ptr<PanNDE::DataBundle<double>> Tabular data to write
      * \param meta std::shared_ptr<PanNDE::DataBundle<double>> Additional metadata to include
      * \param write_index int Time index to associate with this output
      */
      virtual void writeTable(std::string filename_base,
                      std::shared_ptr<PanNDE::DataBundle<double>> data,
                      std::shared_ptr<PanNDE::DataBundle<double>> meta,
                      int write_index) =0;
      
      /*!
      * Writes tabular data and metadata to file.
      * \param filename_base std::string Base name for the output file
      * \param data std::shared_ptr<PanNDE::DataBundle<double>> Tabular data to write
      * \param meta std::shared_ptr<PanNDE::DataBundle<double>> Additional metadata to include
      */
      virtual void writeTable(std::string filename_base,
                      std::shared_ptr<PanNDE::DataBundle<double>> data,
                      std::shared_ptr<PanNDE::DataBundle<double>> meta) =0;
      
      /*!
      * Writes tabular data to file with a time index.
      * \param filename_base std::string Base name for the output file
      * \param data std::shared_ptr<PanNDE::DataBundle<double>> Tabular data to write
      * \param write_index int Time index to associate with this output
      */
      virtual void writeTable(std::string filename_base,
                      std::shared_ptr<PanNDE::DataBundle<double>> data,
                      int write_index) =0;
      
      /*!
      * Writes tabular data to file.
      * \param filename_base std::string Base name for the output file
      * \param data std::shared_ptr<PanNDE::DataBundle<double>> Tabular data to write
      */
      virtual void writeTable(std::string filename_base,
                      std::shared_ptr<PanNDE::DataBundle<double>> data) =0;
  };
};