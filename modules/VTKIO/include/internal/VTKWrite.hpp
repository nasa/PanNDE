/*! \headerfile VTKWrite.hpp "modules/VTKIO/include/internal/VTKWrite.hpp"
* "VTKWrite.hpp" contains the abstract base class for writing PanNDE data to VTK files.
* Defines interfaces for writing both solution data and tabular data with various options.
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

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLTableWriter.h>
#include <vtkDenseArray.h>
#include <vtkArrayData.h>
#include <vtkXMLDataSetWriter.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkMultiBlockDataSet.h>

#include "Mesh.hpp"
#include "Field.hpp"
#include "Array.hpp"

#include "VTKGrid.hpp"
#include "VTKTable.hpp"

namespace VTKIO {
  /*! \class VTKWrite VTKWrite.hpp "modules/VTKIO/include/internal/VTKWrite.hpp"
  *
  * Abstract base class that defines the interface for writing PanNDE data to VTK files.
  * Provides methods for writing both solution data (as VTK unstructured grids) and 
  * tabular data (as VTK tables), with options for including metadata and time indices.
  * Concrete implementations must provide specific writers for grid and table data.
  *
  */
  class VTKWrite {
    public:
      /*static std::shared_ptr<VTKWrite> makeShared(){
        auto retobj=std::make_shared<VTKWrite>(VTKWrite());
        return std::move(retobj);
      }*/

      /*!
      * Writes solution data with metadata to a file with time index.
      * \param filename_base std::string Base name for the output file
      * \param solution std::shared_ptr<PanNDE::FieldBundle> Solution data to write
      * \param meta std::shared_ptr<PanNDE::DataBundle<double>> Metadata to include in the file
      * \param write_index int Time index to append to the filename
      */
      void write(std::string filename_base,
                 std::shared_ptr<PanNDE::FieldBundle> solution,
                 std::shared_ptr<PanNDE::DataBundle<double>> meta,
                 int write_index){
        grid->setMeta(meta);
        write(filename_base,solution,write_index);
      };

      /*!
      * Writes solution data to a file with time index.
      * \param filename_base std::string Base name for the output file
      * \param solution std::shared_ptr<PanNDE::FieldBundle> Solution data to write
      * \param write_index int Time index to append to the filename
      */
      void write(std::string filename_base,
                 std::shared_ptr<PanNDE::FieldBundle> solution,
                 int write_index){
        auto fname=assembleFilename(filename_base,write_index);
        write(fname,solution);
      };

      /*!
      * Writes solution data with optional metadata to a file.
      * \param filename_base std::string Base name for the output file
      * \param solution std::shared_ptr<PanNDE::FieldBundle> Solution data to write
      * \param meta std::shared_ptr<PanNDE::DataBundle<double>> Optional metadata to include (default: nullptr)
      */
      void write(std::string filename_base,
                 std::shared_ptr<PanNDE::FieldBundle> solution,
                 std::shared_ptr<PanNDE::DataBundle<double>> meta=nullptr){
        if(nullptr!=meta){grid->setMeta(meta);};
        grid->setSolution(solution);
        auto vtkgrid=grid->getVTKObject();
        auto writer=getGridWriter();
        commitToFile(filename_base,writer,vtkgrid);
      };

      /*!
      * Writes tabular data with metadata to a file with time index.
      * \param filename_base std::string Base name for the output file
      * \param data std::shared_ptr<PanNDE::DataBundle<double>> Tabular data to write
      * \param meta std::shared_ptr<PanNDE::DataBundle<double>> Metadata to include in the file
      * \param write_index int Time index to append to the filename
      */
      void write(std::string filename_base,
                 std::shared_ptr<PanNDE::DataBundle<double>> data,
                 std::shared_ptr<PanNDE::DataBundle<double>> meta,
                 int write_index){
        table->setMeta(meta);
        write(filename_base,data,write_index);
      };

      /*!
      * Writes tabular data to a file with time index.
      * \param filename_base std::string Base name for the output file
      * \param data std::shared_ptr<PanNDE::DataBundle<double>> Tabular data to write
      * \param write_index int Time index to append to the filename
      */
      void write(std::string filename_base,
                 std::shared_ptr<PanNDE::DataBundle<double>> data,
                 int write_index){
        auto fname=assembleFilename(filename_base,write_index);
        write(fname,data);
      };

      /*!
      * Writes tabular data with optional metadata to a file.
      * \param filename_base std::string Base name for the output file
      * \param data std::shared_ptr<PanNDE::DataBundle<double>> Tabular data to write
      * \param meta std::shared_ptr<PanNDE::DataBundle<double>> Optional metadata to include (default: nullptr)
      */
      void write(std::string filename_base,
                 std::shared_ptr<PanNDE::DataBundle<double>> data,
                 std::shared_ptr<PanNDE::DataBundle<double>> meta=nullptr){
        if(nullptr!=meta){table->setMeta(meta);};
        table->setData(data);
        auto vtktable=table->getVTKObject();
        auto writer=getTableWriter();
        commitToFile(filename_base,writer,vtktable);
      };

    private:
      /*!
      * Pure virtual method to get a VTK table writer.
      * \return vtkSmartPointer<vtkXMLWriter> Writer for VTK table data
      */
      virtual vtkSmartPointer<vtkXMLWriter> getTableWriter()=0;
      
      /*!
      * Pure virtual method to get a VTK grid writer.
      * \return vtkSmartPointer<vtkXMLWriter> Writer for VTK grid data
      */
      virtual vtkSmartPointer<vtkXMLWriter> getGridWriter()=0;

      /*!
      * Commits data to a file using the specified writer.
      * \param filename_base std::string Base filename for output
      * \param writer vtkSmartPointer<vtkXMLWriter> Writer to use for output
      * \param dataobj vtkSmartPointer<vtkDataObject> Data object to write
      */
      void commitToFile(std::string filename_base,
                        vtkSmartPointer<vtkXMLWriter> writer,
                        vtkSmartPointer<vtkDataObject> dataobj){
        auto filename=appendExtension(filename_base,writer);
        writer->SetFileName(filename.c_str());
        writer->SetInputData(0,dataobj);
        writer->Update();
        writer->Write();
      };

      /*!
      * Assembles a filename with time index.
      * \param filename_base std::string Base filename
      * \param write_index int Time index to append
      * \return std::string The assembled filename
      */
      inline std::string assembleFilename(std::string filename_base,int write_index){
        auto fname=filename_base + "_" + std::to_string(write_index);
        return fname;
      }

      /*!
      * Appends the appropriate file extension based on writer type.
      * \param filename_base std::string Base filename
      * \param writer vtkSmartPointer<vtkXMLWriter> Writer that determines the extension
      * \return std::string Filename with appropriate extension
      */
      inline std::string appendExtension(std::string filename_base,
                                         vtkSmartPointer<vtkXMLWriter> writer){
        auto ext=writer->GetDefaultFileExtension();
        auto fname=filename_base + "." + std::string(ext);
        return fname;
      }

      //! VTKGrid for handling mesh and field data conversion
      std::shared_ptr<VTKGrid> grid=VTKGrid::makeShared();
      
      //! VTKTable for handling tabular data conversion
      std::shared_ptr<VTKTable> table=VTKTable::makeShared();
  };
};