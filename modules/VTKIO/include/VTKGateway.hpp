/*! \headerfile VTKGateway.hpp "modules/VTKIO/include/VTKGateway.hpp"
* "VTKGateway.hpp" contains the class implementation encapsulating the file input/output 
* management using VTK.
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
#include <string>

#include "Gateway.hpp"
#include "internal/VTKRead.hpp"
#include "internal/VTKWrite.hpp"
#include "internal/VTKPWrite.hpp"
#include "internal/VTKSWrite.hpp"

#include "Array.hpp"
#include "Mesh.hpp"
#include "Field.hpp"

namespace VTKIO {
  /*! \class VTKGateway VTKGateway.hpp "modules/VTKIO/include/VTKGateway.hpp"
  *
  * Implements the methods required to describe the data for file I/O using VTK.
  * Provides a unified interface for both serial and parallel read/write operations.
  *
  */
  class VTKGateway : public Controller::Gateway {
    public:
      /*!
      * Constructor for VTKGateway.
      * \param communicator std::shared_ptr<PanNDE::Communicator> Provides communication required for parallel file writing
      */
      VTKGateway(std::shared_ptr<PanNDE::Communicator> communicator=nullptr){
        reader=VTKIO::VTKRead::makeShared();
        comm=communicator;
        //setWriter();
      };

      /*!
      * Creates a shared pointer to a VTKGateway object.
      * \param communicator std::shared_ptr<PanNDE::Communicator> Communicator for parallel operations (default: nullptr)
      * \return std::shared_ptr<VTKIO::VTKGateway> Shared pointer to a newly created VTKGateway
      */
      static inline
      std::shared_ptr<VTKIO::VTKGateway> makeShared(std::shared_ptr<PanNDE::Communicator> communicator=nullptr){
        auto file_IO=std::make_shared<VTKIO::VTKGateway>(VTKIO::VTKGateway(communicator));
        return std::move(file_IO);
      };

      /*!
      * Opens a file for reading.
      * \param filename std::string The name of the file (using relative address)
      */
      void open(std::string filename)override{
        (nullptr==comm)?(reader->open(filename)):masterOpen(filename);
        return;
      };

      /*!
      * Gets array by keyname from file.
      * \param keyname std::string Name of variable to be pulled from input file
      * \param maker std::shared_ptr<PanNDE::ArrayFactory<double>> Array factory for creating the result
      * \return std::shared_ptr<PanNDE::Array<double>> The array data from the file
      */
      std::shared_ptr<PanNDE::Array<double>> getArray(std::string keyname,
                              std::shared_ptr<PanNDE::ArrayFactory<double>> maker)override{
        if(nullptr==comm){
          return std::move(getSerialArray(keyname,maker));
        }else{return std::move(getParallelArray(keyname,maker));};
      };

      /*!
      * Gets value by keyname from file.
      * \param keyname std::string Name of variable to be pulled from input file
      * \return double The value from the file
      */
      double getValue(std::string keyname)override{
        if(nullptr==comm){
          return std::move(getSerialValue(keyname));
        }else{return std::move(getParallelValue(keyname));};
      };

      /*!
      * Gets mesh from file.
      * \param maker std::shared_ptr<PanNDE::MeshFactory> Mesh factory for creating the result
      * \param host_id int Root process to read mesh from file (default: 0)
      * \return std::shared_ptr<PanNDE::Mesh> The mesh from the file, or nullptr if not on host_id
      */
      std::shared_ptr<PanNDE::Mesh> getMesh(
                              std::shared_ptr<PanNDE::MeshFactory> maker,int host_id=0)override{
        if(nullptr==comm){
          return readMesh(maker);
        }else if(comm->getProcessId()==host_id){
          return readMesh(maker);
        }else{return nullptr;};
      };

      /*!
      * Gets field by keyname from file.
      * \param keyname std::string Name of field to be pulled from input file
      * \param maker std::shared_ptr<PanNDE::FieldFactory> Field factory for creating the result
      * \param host_id int Root process to read field from file (default: 0)
      * \return std::shared_ptr<PanNDE::Field> The field from the file, or nullptr if not on host_id
      */
      std::shared_ptr<PanNDE::Field> getField(std::string keyname,
                              std::shared_ptr<PanNDE::FieldFactory> maker,int host_id=0)override{
        if(nullptr==comm){
          return std::move(readField(keyname,maker));
        }else if(comm->getProcessId()==host_id){
          return std::move(readField(keyname,maker));
        }else{return nullptr;};
      };
      
      /*!
      * Writes solution variables to file with time index.
      * \param filename_base std::string Base name for file write
      * \param solution std::shared_ptr<PanNDE::FieldBundle> Bundled solution data
      * \param write_index int Time index to be written
      */
      void writeSolution(std::string filename_base,
                         std::shared_ptr<PanNDE::FieldBundle> solution,
                         int write_index)override{
        auto writer=setWriter();
        writer->write(filename_base,solution,write_index);
      };

      /*!
      * Writes solution variables to file with metadata and time index.
      * \param filename_base std::string Base name for file write
      * \param solution std::shared_ptr<PanNDE::FieldBundle> Bundled solution data
      * \param meta std::shared_ptr<PanNDE::DataBundle<double>> Additional metadata to be written
      * \param write_index int Time index to be written
      */
      void writeSolution(std::string filename_base,
                         std::shared_ptr<PanNDE::FieldBundle> solution,
                         std::shared_ptr<PanNDE::DataBundle<double>> meta,
                         int write_index)override{
        auto writer=setWriter();
        writer->write(filename_base,solution,meta,write_index);
      };

      /*!
      * Writes solution variables to file with metadata.
      * \param filename_base std::string Base name for file write
      * \param solution std::shared_ptr<PanNDE::FieldBundle> Bundled solution data
      * \param meta std::shared_ptr<PanNDE::DataBundle<double>> Additional metadata to be written
      */
      void writeSolution(std::string filename_base,
                         std::shared_ptr<PanNDE::FieldBundle> solution,
                         std::shared_ptr<PanNDE::DataBundle<double>> meta)override{
        auto writer=setWriter();
        writer->write(filename_base,solution,meta);
      };

      /*!
      * Writes solution variables to file.
      * \param filename_base std::string Base name for file write
      * \param solution std::shared_ptr<PanNDE::FieldBundle> Bundled solution data
      */
      void writeSolution(std::string filename_base,
                         std::shared_ptr<PanNDE::FieldBundle> solution)override{
        auto writer=setWriter();
        writer->write(filename_base,solution);
      };

      /*!
      * Writes tabular data to file with metadata and time index.
      * \param filename_base std::string Base name for file write
      * \param data std::shared_ptr<PanNDE::DataBundle<double>> Data to be written
      * \param meta std::shared_ptr<PanNDE::DataBundle<double>> Additional metadata to be written
      * \param write_index int Time index to be written
      */
      void writeTable(std::string filename_base,
                      std::shared_ptr<PanNDE::DataBundle<double>> data,
                      std::shared_ptr<PanNDE::DataBundle<double>> meta,
                      int write_index)override{
        auto writer=setWriter();
        writer->write(filename_base,data,meta,write_index);
      };

      /*!
      * Writes tabular data to file with metadata.
      * \param filename_base std::string Base name for file write
      * \param data std::shared_ptr<PanNDE::DataBundle<double>> Data to be written
      * \param meta std::shared_ptr<PanNDE::DataBundle<double>> Additional metadata to be written
      */
      void writeTable(std::string filename_base,
                      std::shared_ptr<PanNDE::DataBundle<double>> data,
                      std::shared_ptr<PanNDE::DataBundle<double>> meta)override{
        auto writer=setWriter();
        writer->write(filename_base,data,meta);
      };

      /*!
      * Writes tabular data to file with time index.
      * \param filename_base std::string Base name for file write
      * \param data std::shared_ptr<PanNDE::DataBundle<double>> Data to be written
      * \param write_index int Time index to be written
      */
      void writeTable(std::string filename_base,
                      std::shared_ptr<PanNDE::DataBundle<double>> data,
                      int write_index)override{
        auto writer=setWriter();
        writer->write(filename_base,data,write_index);
      };

      /*!
      * Writes tabular data to file.
      * \param filename_base std::string Base name for file write
      * \param data std::shared_ptr<PanNDE::DataBundle<double>> Data to be written
      */
      void writeTable(std::string filename_base,
                      std::shared_ptr<PanNDE::DataBundle<double>> data)override{
        auto writer=setWriter();
        writer->write(filename_base,data);
      };

    private:
      /*!
      * Opens a file on the master process only.
      * \param filename std::string The name of the file to open
      */
      inline void masterOpen(std::string filename){
        if(0==comm->getProcessId()){reader->open(filename);};
      };

      /*!
      * Creates an appropriate writer based on the communicator.
      * \return std::shared_ptr<VTKIO::VTKWrite> The selected writer implementation
      */
      std::shared_ptr<VTKIO::VTKWrite> setWriter(){
        std::shared_ptr<VTKIO::VTKWrite> writer=nullptr;
        if(nullptr==comm or comm->getNumberOfProcesses()==1){
          writer=VTKIO::VTKSWrite::makeShared();
        }else{
          writer=VTKIO::VTKPWrite::makeShared(comm);
        };
        return std::move(writer);
      };

      /*!
      * Gets array from file in serial mode.
      * \param keyname std::string Name of array to retrieve
      * \param maker std::shared_ptr<PanNDE::ArrayFactory<double>> Factory for creating the array
      * \return std::shared_ptr<PanNDE::Array<double>> The retrieved array
      */
      std::shared_ptr<PanNDE::Array<double>> getSerialArray(std::string keyname,
                              std::shared_ptr<PanNDE::ArrayFactory<double>> maker){
        return reader->getArray(keyname,maker);
      };

      /*!
      * Gets array from file in parallel mode.
      * \param keyname std::string Name of array to retrieve
      * \param maker std::shared_ptr<PanNDE::ArrayFactory<double>> Factory for creating the array
      * \return std::shared_ptr<PanNDE::Array<double>> The retrieved array, broadcast to all processes
      */
      std::shared_ptr<PanNDE::Array<double>> getParallelArray(std::string keyname,
                              std::shared_ptr<PanNDE::ArrayFactory<double>> maker){
        std::shared_ptr<PanNDE::Array<double>> array=nullptr;
        if(0==comm->getProcessId()){
          array=reader->getArray(keyname,maker);
        }else{array=maker->makeManagedArray();};
        comm->broadcastArray(array,0);
        return std::move(array);
      };

      /*!
      * Gets value from file in serial mode.
      * \param keyname std::string Name of value to retrieve
      * \return double The retrieved value
      */
      double getSerialValue(std::string keyname){
        return reader->getValue(keyname);
      };

      /*!
      * Gets value from file in parallel mode.
      * \param keyname std::string Name of value to retrieve
      * \return double The retrieved value, broadcast to all processes
      */
      double getParallelValue(std::string keyname){
        double value;
        if(0==comm->getProcessId()){value=reader->getValue(keyname);};
        comm->broadcastValue(&value,0);
        return std::move(value);
      };

      /*!
      * Reads mesh from file.
      * \param maker std::shared_ptr<PanNDE::MeshFactory> Factory for creating the mesh
      * \return std::shared_ptr<PanNDE::Mesh> The retrieved mesh
      */
      std::shared_ptr<PanNDE::Mesh> readMesh(std::shared_ptr<PanNDE::MeshFactory> maker){
        auto mesh=reader->getMesh(maker);
        return mesh;
      };

      /*!
      * Reads field from file.
      * \param keyname std::string Name of field to retrieve
      * \param maker std::shared_ptr<PanNDE::FieldFactory> Factory for creating the field
      * \return std::shared_ptr<PanNDE::Field> The retrieved field
      */
      std::shared_ptr<PanNDE::Field> readField(std::string keyname,
                                               std::shared_ptr<PanNDE::FieldFactory> maker){
        auto field=reader->getField(keyname,maker);
        return std::move(field);
      };

      //! Communicator for parallel operations
      std::shared_ptr<PanNDE::Communicator> comm=nullptr;

      //! Reader for VTK file operations
      std::shared_ptr<VTKIO::VTKRead> reader=nullptr;
  };
};
