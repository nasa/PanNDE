/*! \headerfile VTKGateway.hpp "modules/VTKIO/include/VTKGateway.hpp"
* "VTKGateway.hpp" contains the class implementation encapsulating the file input/output management using VTK
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

#include "Array.hpp"
#include "Mesh.hpp"
#include "Field.hpp"

namespace VTKIO {
  /*! \class VTKGateway VTKGateway.hpp "modules/VTKIO/include/VTKGateway.hpp"
  *
  * Implements the methods required to describe the data for file I/O 
  *
  */
  class VTKGateway : public Controller::Gateway {
    public:
      /*!
      * constructor
      * \param communcator provides communication required for parallel file writing
      */
      VTKGateway(std::shared_ptr<PanNDE::Communicator> communicator=nullptr){
        reader=makeShared<VTKIO::VTKRead>();
        comm=communicator;
        //setWriter();
      };

      /*!
      * open file for read
      * \param filename the name of the file (using relative addr)
      */
      void open(std::string filename)override{
        reader->open(filename);
      };

      /*!
      * get array by keyname from file
      * \param keyname name of variable to be pulled from input file
      * \param maker array factory
      */
      std::shared_ptr<PanNDE::Array<double>> getArray(std::string keyname,
                              std::shared_ptr<PanNDE::ArrayFactory<double>> maker)override{
        if(nullptr==comm){
          return std::move(getSerialArray(keyname,maker));
        }else{return std::move(getParallelArray(keyname,maker));};
      };
      /*!
      * get value by keyname from file
      * \param keyname name of variable to be pulled from input file
      */
      double getValue(std::string keyname)override{if(nullptr==comm){
          return std::move(getSerialValue(keyname));
        }else{return std::move(getParallelValue(keyname));};
      };

      /*!
      * get mesh from file
      * \param maker mesh factory
      * \param host_id root process to read mesh from file
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
      * get field by keyname from file
      * \param maker field factory
      * \param host_id root process to read mesh from file
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
      * write solution variables to file
      * \param filename_base base name for file write
      * \param solution bundled solution data
      * \param write_index time index to be written
      */
      void writeSolution(std::string filename_base,
                         std::shared_ptr<PanNDE::FieldBundle> solution,
                         int write_index)override{
        auto writer=setWriter();
        writer->write(filename_base,solution,write_index);
      };
      /*!
      * write solution variables to file
      * \param filename_base base name for file write
      * \param solution bundled solution data
      * \param meta additional metadata to be written to file
      * \param write_index time index to be written
      */
      void writeSolution(std::string filename_base,
                         std::shared_ptr<PanNDE::FieldBundle> solution,
                         std::shared_ptr<PanNDE::DataBundle<double>> meta,
                         int write_index)override{
        auto writer=setWriter();
        writer->write(filename_base,solution,meta,write_index);
      };
      /*!
      * write solution variables to file
      * \param filename_base base name for file write
      * \param solution bundled solution data
      * \param meta additional metadata to be written to file
      */
      void writeSolution(std::string filename_base,
                         std::shared_ptr<PanNDE::FieldBundle> solution,
                         std::shared_ptr<PanNDE::DataBundle<double>> meta)override{
        auto writer=setWriter();
        writer->write(filename_base,solution,meta);
      };
      /*!
      * write solution variables to file
      * \param filename_base base name for file write
      * \param solution bundled solution data
      */
      void writeSolution(std::string filename_base,
                         std::shared_ptr<PanNDE::FieldBundle> solution)override{
        auto writer=setWriter();
        writer->write(filename_base,solution);
      };

    private:
      std::shared_ptr<VTKIO::VTKWrite> setWriter(){
        std::shared_ptr<VTKIO::VTKWrite> writer=nullptr;
        if(nullptr==comm){
          writer=makeShared<VTKIO::VTKWrite>();
        }else{
          writer=std::make_shared<VTKIO::VTKPWrite>(VTKIO::VTKPWrite(comm));
        };
        return std::move(writer);
      };
      template<class T>
      static inline std::shared_ptr<T> makeShared(){
        return std::move(std::make_shared<T>(T()));
      };
      std::shared_ptr<PanNDE::Array<double>> getSerialArray(std::string keyname,
                              std::shared_ptr<PanNDE::ArrayFactory<double>> maker){
        return reader->getArray(keyname,maker);
      };
      std::shared_ptr<PanNDE::Array<double>> getParallelArray(std::string keyname,
                              std::shared_ptr<PanNDE::ArrayFactory<double>> maker){
        std::shared_ptr<PanNDE::Array<double>> array=nullptr;
        if(0==comm->getProcessId()){
          array=reader->getArray(keyname,maker);
        }else{array=maker->makeManagedArray();};
        comm->broadcastArray(array,0);
        return std::move(array);
      };
      double getSerialValue(std::string keyname){
        return reader->getValue(keyname);
      };
      double getParallelValue(std::string keyname){
        double value;
        if(0==comm->getProcessId()){value=reader->getValue(keyname);};
        comm->broadcastValue(&value,0);
        return std::move(value);
      };
      std::shared_ptr<PanNDE::Mesh> readMesh(std::shared_ptr<PanNDE::MeshFactory> maker){
        auto mesh=reader->getMesh(maker);
        return mesh;
      };
      std::shared_ptr<PanNDE::Field> readField(std::string keyname,
                                               std::shared_ptr<PanNDE::FieldFactory> maker){
        auto field=reader->getField(keyname,maker);
        return std::move(field);
      };

      std::shared_ptr<PanNDE::Communicator> comm=nullptr;

      std::shared_ptr<VTKIO::VTKRead> reader=nullptr;
  };
};
