/*! \headerfile VTKGrid.hpp "modules/VTKIO/include/internal/VTKGrid.hpp"
* "VTKGrid.hpp" contains the class implementation for converting PanNDE mesh
* and field data to VTK unstructured grid format for visualization and analysis.
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
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>

//#include <vtkDataObject.h>

#include "Mesh.hpp"
#include "Field.hpp"
#include "Array.hpp"

#include "VTKMeta.hpp"

namespace VTKIO {
  /*! \class VTKGrid VTKGrid.hpp "modules/VTKIO/include/internal/VTKGrid.hpp"
  *
  * Implements conversion from PanNDE mesh and field data to VTK unstructured grid format
  * for visualization and data analysis.
  *
  */
  class VTKGrid{
    public:
      /*!
      * Creates a shared pointer to a VTKGrid object.
      * \return std::shared_ptr<VTKGrid> Shared pointer to a newly created VTKGrid object
      */
      static std::shared_ptr<VTKGrid> makeShared(){
        auto grid=std::make_shared<VTKGrid>(VTKGrid());
        return std::move(grid);
      };

      /*!
      * Sets the solution data to be represented in the VTK grid.
      * \param solution std::shared_ptr<PanNDE::FieldBundle> The solution data to convert
      */
      void setSolution(std::shared_ptr<PanNDE::FieldBundle> solution){
        setMesh(solution->mesh());
        setFieldData(solution);
      };

      /*!
      * Sets metadata for the VTK grid.
      * \param meta std::shared_ptr<PanNDE::DataBundle<double>> The metadata to set
      */
      void setMeta(std::shared_ptr<PanNDE::DataBundle<double>> meta){
        meta_ops->setMetaData(uGridobj,meta);
      };

      /*!
      * Gets the VTK unstructured grid object.
      * \return vtkSmartPointer<vtkUnstructuredGrid> The VTK unstructured grid representation
      */
      vtkSmartPointer<vtkUnstructuredGrid> getVTKObject(){return uGridobj;};

    private:
      /*!
      * Validates that the mesh has valid dimensions.
      * \param mesh std::shared_ptr<PanNDE::Mesh> The mesh to check
      * \throws std::runtime_error If the mesh has invalid dimensions
      */
      void meshCheck(std::shared_ptr<PanNDE::Mesh> mesh){
        if(0>=mesh->nodeCount()){throw std::runtime_error("invalid mesh size (node count)");};
        if(0>=mesh->cellCount()){throw std::runtime_error("invalid mesh size (cell count)");};
      };

      /*!
      * Adds field data from a FieldBundle to the VTK grid.
      * \param solution std::shared_ptr<PanNDE::FieldBundle> The solution with field data
      * \throws std::runtime_error If there's a field-mesh size mismatch or unsupported field type
      */
      void setFieldData(std::shared_ptr<PanNDE::FieldBundle> solution){
        std::string keyname;
        for(int k=0;k<solution->fieldCount();k++){
          keyname=solution->fieldName(k);
          auto field=solution->field(keyname);
          auto data=extractFieldData(keyname,field);
          meshCheck(solution->mesh());
          if(PanNDE::Field::NODE==field->type()){
            if(field->size()!=solution->mesh()->nodeCount()){
              throw std::runtime_error("field-mesh (node) mismatch: " +
                                        std::to_string(field->size()) +
                                       " vs " +
                                        std::to_string(solution->mesh()->nodeCount()));
            };
            uGridobj->GetPointData()->AddArray(data);
          }else if(PanNDE::Field::CELL==field->type()){
            if(field->size()!=solution->mesh()->cellCount()){
              throw std::runtime_error("field-mesh (cell) mismatch: " +
                                        std::to_string(field->size()) +
                                       " vs " +
                                        std::to_string(solution->mesh()->cellCount()));
            };
            uGridobj->GetCellData()->AddArray(data);
          }else{throw std::runtime_error("Unsupported Field Type (not node/cell)");};
        };
      };

      /*!
      * Extracts field data into a VTK double array.
      * \param keyname std::string The name of the field
      * \param field std::shared_ptr<PanNDE::Field> The field containing the data
      * \return vtkSmartPointer<vtkDoubleArray> VTK array containing the field data
      */
      vtkSmartPointer<vtkDoubleArray> extractFieldData(std::string keyname,
                                                       std::shared_ptr<PanNDE::Field> field){
        vtkSmartPointer<vtkDoubleArray> data=vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfComponents(1);
        data->SetName(keyname.c_str());
        for(int k=0;k<field->size();k++){
          data->InsertNextValue(field->at(k));
        };
        return std::move(data);
      };

      /*!
      * Sets the mesh for the VTK grid.
      * \param mesh std::shared_ptr<PanNDE::Mesh> The mesh to set
      */
      void setMesh(std::shared_ptr<PanNDE::Mesh> mesh){
        uGridobj->Reset();
        setVTKPoints(mesh);
        setVTKCells(mesh);
      };

      /*!
      * Sets the node points for the VTK grid.
      * \param mesh std::shared_ptr<PanNDE::Mesh> The mesh containing node coordinates
      */
      void setVTKPoints(std::shared_ptr<PanNDE::Mesh> mesh){
        double coords[3];
        for(int32_t kn=0;kn<mesh->nodeCount();kn++){
          mesh->nodeCoordinate(kn,coords);
          nodeLocations->InsertNextPoint(coords);
        };
        uGridobj->SetPoints(nodeLocations);
      };

      /*!
      * Sets the cells for the VTK grid.
      * \param mesh std::shared_ptr<PanNDE::Mesh> The mesh containing cell connectivity
      */
      void setVTKCells(std::shared_ptr<PanNDE::Mesh> mesh){
        int32_t cell_nodes[8];
        for(int32_t kc=0;kc<mesh->cellCount();kc++){
          mesh->cell(kc,cell_nodes);
          for(int k=0;k<8;k++){
            hexCell->GetPointIds()->SetId(k,cell_nodes[k]);
          };
          uGridobj->InsertNextCell(hexCell->GetCellType(),
                                   hexCell->GetPointIds());
        };
        uGridobj->BuildLinks();
      };

      //! Helper object for handling VTK metadata operations
      std::shared_ptr<VTKMeta> meta_ops=VTKMeta::makeShared();

      //! Storage for node coordinates in VTK format
      vtkSmartPointer<vtkPoints> nodeLocations=
                      vtkSmartPointer<vtkPoints>::New();
      //! Template for hexahedral cell
      vtkSmartPointer<vtkHexahedron> hexCell =
                      vtkSmartPointer<vtkHexahedron>::New();
      //! The VTK unstructured grid representation
      vtkSmartPointer<vtkUnstructuredGrid> uGridobj =
                      vtkSmartPointer<vtkUnstructuredGrid>::New();
  };
};