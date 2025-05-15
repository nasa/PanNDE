/*! \headerfile VTKRead.hpp "modules/VTKIO/include/internal/VTKRead.hpp"
* "VTKRead.hpp" contains the class implementation for reading VTK unstructured grid
* files and extracting data elements such as arrays, meshes, and fields.
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
#include <vtkCell.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

#include "Array.hpp"
#include "Mesh.hpp"
#include "Field.hpp"

namespace VTKIO {
  /*! \class VTKRead VTKRead.hpp "modules/VTKIO/include/internal/VTKRead.hpp"
  *
  * Implements functionality for reading VTK unstructured grid files (.vtu).
  * Provides capabilities to extract arrays, scalar values, meshes, and fields
  * from VTK files into the corresponding PanNDE data structures.
  *
  */
  class VTKRead {
    public:
      /*!
      * Creates a shared pointer to a VTKRead object.
      * \return std::shared_ptr<VTKRead> Shared pointer to a newly created VTKRead object
      */
      static std::shared_ptr<VTKRead> makeShared(){
        auto retobj=std::make_shared<VTKRead>(VTKRead());
        return std::move(retobj);
      };

      /*!
      * Opens a VTK file for reading.
      * \param filename std::string Path to the VTK file to open
      */
      void open(std::string filename){
        current_file=filename;
        loadFile();
      };

      /*!
      * Gets an array from the file by name.
      * \param keyname std::string Name of the array to retrieve
      * \param maker std::shared_ptr<PanNDE::ArrayFactory<double>> Factory for creating the array
      * \return std::shared_ptr<PanNDE::Array<double>> The retrieved array data
      * \throws std::runtime_error If no file is open or the named array doesn't exist
      */
      std::shared_ptr<PanNDE::Array<double>> getArray(std::string keyname,
                              std::shared_ptr<PanNDE::ArrayFactory<double>> maker){
        ifNoFileError();
        return std::move(get1DArray(keyname,maker));
      };

      /*!
      * Gets a scalar value from the file by name.
      * \param keyname std::string Name of the value to retrieve
      * \return double The retrieved scalar value
      * \throws std::runtime_error If no file is open, the named value doesn't exist, or is not a scalar
      */
      double getValue(std::string keyname){
        ifNoFileError();
        auto vec=get1Dvector(keyname);
        if(1!=vec.size()){throw std::runtime_error("not scalar quantity");};
        return vec.at(0);
      };

      /*!
      * Gets mesh data from the file.
      * \param maker std::shared_ptr<PanNDE::MeshFactory> Factory for creating the mesh
      * \return std::shared_ptr<PanNDE::Mesh> The mesh constructed from the file data
      * \throws std::runtime_error If no file is open
      */
      std::shared_ptr<PanNDE::Mesh> getMesh(std::shared_ptr<PanNDE::MeshFactory> maker){
        ifNoFileError();
        if(nullptr!=mesh){
          return mesh;
        }else{
          auto nodes=getMeshPoints(maker);
          auto cells=getMeshCells(maker);
          mesh=maker->makeManagedMesh(nodes.data(),nodes.size(),cells.data(),cells.size());
          return mesh;
        };
      };

      /*!
      * Gets a field from the file by name.
      * \param keyname std::string Name of the field to retrieve
      * \param maker std::shared_ptr<PanNDE::FieldFactory> Factory for creating the field
      * \return std::shared_ptr<PanNDE::Field> The field constructed from the file data
      * \throws std::runtime_error If no file is open, no mesh is loaded, or the field doesn't exist
      */
      std::shared_ptr<PanNDE::Field> getField(std::string keyname,
                                              std::shared_ptr<PanNDE::FieldFactory> maker){
        ifNoFileError();
        ifNoMeshError();
        if(ugrid->GetPointData()->HasArray(keyname.c_str())){
          return std::move(getNodeField(keyname,maker));
        }else if(ugrid->GetCellData()->HasArray(keyname.c_str())){
          return std::move(getCellField(keyname,maker));
        }else{
          throwNoName();return nullptr;
        };
      };

      /*!
      * Gets the currently active file path.
      * \return std::string Path of the currently open file
      */
      std::string currentlyActiveFile(){return current_file;};

    private:
      /*!
      * Throws an error if no file is currently open.
      * \throws std::runtime_error If no file is open
      */
      void ifNoFileError(){if(current_file.empty()){throw std::runtime_error("no file provided");};};
      
      /*!
      * Throws an error if no mesh is currently loaded.
      * \throws std::runtime_error If no mesh is loaded
      */
      void ifNoMeshError(){if(nullptr==mesh){throw std::runtime_error("no mesh loaded");};};
      
      /*!
      * Throws an error for non-existent data names.
      * \throws std::runtime_error With message about missing data
      */
      void throwNoName(){throw std::runtime_error("no data by this name present");};
      
      /*!
      * Loads the file and updates the internal data structures.
      */
      void loadFile(){
        vtkreader->SetFileName(current_file.c_str());
        vtkreader->Update();
        ugrid=vtkreader->GetOutput();
        mesh=nullptr;
      };
      
      /*!
      * Gets a 1D array from the file as a PanNDE::Array.
      * \param name std::string Name of the array to retrieve
      * \param maker std::shared_ptr<PanNDE::ArrayFactory<double>> Factory for creating the array
      * \return std::shared_ptr<PanNDE::Array<double>> The retrieved array data
      */
      std::shared_ptr<PanNDE::Array<double>> get1DArray(std::string name,
                                     std::shared_ptr<PanNDE::ArrayFactory<double>> maker){
        auto array=maker->makeManagedArray();
        auto adata=get1Dvector(name);
        array->resize(adata.size());
        for(int k=0;k<array->size();k++){array->at(k)=adata.at(k);};
        return std::move(array);
      };
      
      /*!
      * Gets a 1D array from the file as a vector.
      * \param name std::string Name of the array to retrieve
      * \return std::vector<double> Vector containing the array data
      * \throws std::runtime_error If the named array doesn't exist
      */
      std::vector<double> get1Dvector(std::string name){
        std::vector<double> adata;
        if(ugrid->GetFieldData()->HasArray(name.c_str())){
          vtkSmartPointer<vtkDataArray> readdata=
                          ugrid->GetFieldData()->GetArray(name.c_str());
          adata.resize(readdata->GetNumberOfTuples());
          for(int k=0;k<adata.size();k++){readdata->GetTuple(k,&adata.at(k));};
          return std::move(adata);
        }else{throwNoName();return std::move(adata);};
      };

      /*!
      * Extracts mesh points (nodes) from the VTK file.
      * \param maker std::shared_ptr<PanNDE::MeshFactory> Factory for creating nodes
      * \return std::vector<PanNDE::Mesh::Node> Vector of nodes from the file
      */
      std::vector<PanNDE::Mesh::Node> getMeshPoints(std::shared_ptr<PanNDE::MeshFactory> maker){
        std::vector<PanNDE::Mesh::Node> nodes;
        nodes.reserve(ugrid->GetNumberOfPoints());
        double pt[3];
        for(int kp=0;kp<ugrid->GetNumberOfPoints();kp++){
          ugrid->GetPoint(kp,pt);
          nodes.push_back(maker->makeNode(pt,kp,0));
        };
        return std::move(nodes);
      };
      
      /*!
      * Extracts mesh cells from the VTK file.
      * \param maker std::shared_ptr<PanNDE::MeshFactory> Factory for creating cells
      * \return std::vector<PanNDE::Mesh::Cell> Vector of cells from the file
      */
      std::vector<PanNDE::Mesh::Cell> getMeshCells(std::shared_ptr<PanNDE::MeshFactory> maker){
        std::vector<PanNDE::Mesh::Cell> cells;
        int32_t box[8];
        cells.reserve(ugrid->GetNumberOfCells());
        for(int kc=0;kc<ugrid->GetNumberOfCells();kc++){
          for(int kb=0;kb<8;kb++){
            box[kb]=ugrid->GetCell(kc)->GetPointId(kb);
          };
          cells.push_back(maker->makeCell(box,8,kc,0));
        };
        return std::move(cells);
      };

      /*!
      * Creates a node-based field from VTK point data.
      * \param keyname std::string Name of the field to retrieve
      * \param maker std::shared_ptr<PanNDE::FieldFactory> Factory for creating the field
      * \return std::shared_ptr<PanNDE::Field> The node-based field
      */
      std::shared_ptr<PanNDE::Field> getNodeField(std::string keyname,
                                                  std::shared_ptr<PanNDE::FieldFactory> maker){
        vtkSmartPointer<vtkDataArray> readdata=
                        ugrid->GetPointData()->GetArray(keyname.c_str());
        auto field=maker->makeManagedField(mesh,PanNDE::Field::NODE);
        double value;
        for(int k=0;k<mesh->nodeCount();k++){
          readdata->GetTuple(k,&value);
          field->at(k)=value;
        };
        return std::move(field);
      };
      
      /*!
      * Creates a cell-based field from VTK cell data.
      * \param keyname std::string Name of the field to retrieve
      * \param maker std::shared_ptr<PanNDE::FieldFactory> Factory for creating the field
      * \return std::shared_ptr<PanNDE::Field> The cell-based field
      */
      std::shared_ptr<PanNDE::Field> getCellField(std::string keyname,
                                                  std::shared_ptr<PanNDE::FieldFactory> maker){
        vtkSmartPointer<vtkDataArray> readdata=
                        ugrid->GetCellData()->GetArray(keyname.c_str());
        auto field=maker->makeManagedField(mesh,PanNDE::Field::CELL);
        double value;
        for(int k=0;k<mesh->cellCount();k++){
          readdata->GetTuple(k,&value);
          field->at(k)=value;
        };
        return std::move(field);
      };

      //! Current mesh extracted from the file
      std::shared_ptr<PanNDE::Mesh> mesh=nullptr;

      //! Path of the currently open file
      std::string current_file;
      
      //! VTK XML unstructured grid reader
      vtkSmartPointer<vtkXMLUnstructuredGridReader> vtkreader=
                      vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
      
      //! VTK unstructured grid containing file data
      vtkSmartPointer<vtkUnstructuredGrid> ugrid=
                      vtkSmartPointer<vtkUnstructuredGrid>::New();
  };
};