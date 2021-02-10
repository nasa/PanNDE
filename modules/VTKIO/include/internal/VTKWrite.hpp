//
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
#include <vtkDenseArray.h>
#include <vtkArrayData.h>
#include <vtkXMLDataSetWriter.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkMultiBlockDataSet.h>

#include "Mesh.hpp"
#include "Field.hpp"
#include "Array.hpp"

namespace VTKIO {
  class VTKWrite {
    public:
      void write(std::string filename_base,
                 std::shared_ptr<PanNDE::FieldBundle> solution,
                 int write_index){
        setFileName(filename_base,write_index);
        setSolution(solution);
        commitToFile();
      };
      void write(std::string filename_base,
                 std::shared_ptr<PanNDE::FieldBundle> solution,
                 std::shared_ptr<PanNDE::DataBundle<double>> meta,
                 int write_index){
        setMeta(meta);
        write(filename_base,solution,write_index);
      };
      void write(std::string filename_no_ext,
                 std::shared_ptr<PanNDE::FieldBundle> solution,
                 std::shared_ptr<PanNDE::DataBundle<double>> meta=nullptr){
        if(nullptr!=meta){setMeta(meta);};
        setSolution(solution);
        setFileName(filename_no_ext);
        commitToFile();
      };

    protected:
      virtual void commitToFile(){
        writer->SetFileName(filename.c_str());
        writer->SetInputData(0,uGridVTK);
        writer->Write();
      };
      virtual void setFileName(std::string filename_base,int write_index){
        filename=filename_base + std::to_string(write_index) + ".vtu";
      };
      virtual void setFileName(std::string filename_base){
        filename=filename_base + ".vtu";
      };
      void setSolution(std::shared_ptr<PanNDE::FieldBundle> solution){
        setMesh(solution->mesh());
        setFieldData(solution);
      };
      void setMeta(std::shared_ptr<PanNDE::DataBundle<double>> meta){
        setScalars(meta);
        setArrays(meta);
      };
      void setScalars(std::shared_ptr<PanNDE::DataBundle<double>> meta){
        double value;
        for(int k=0;k<meta->scalarCount();k++){
          auto keyname=meta->scalarName(k);
          value=meta->scalar(keyname);
          set1DArray(&value,1,keyname);
        };
      };
      void setArrays(std::shared_ptr<PanNDE::DataBundle<double>> meta){
        for(int k=0;k<meta->arrayCount();k++){
          auto keyname=meta->arrayName(k);
          auto array=meta->array(keyname);
          set1DArray(array->data(),array->size(),keyname);
        };
      };
      void set1DArray(double* array,int N,std::string name){
        vtkSmartPointer<vtkDoubleArray> series = vtkSmartPointer<vtkDoubleArray>::New();
        series->SetNumberOfComponents(1);
        series->SetNumberOfTuples(N);
        double value;
        for(int k=0;k<N;k++){
          value=array[k];
          series->SetTuple(k,&value);
        }
        series->SetName(name.c_str());
        uGridVTK->GetFieldData()->AddArray(series);
      };
      void meshCheck(std::shared_ptr<PanNDE::Mesh> mesh){
        if(0>=mesh->nodeCount()){throw std::runtime_error("invalid mesh size");};
        if(0>=mesh->cellCount()){throw std::runtime_error("invalid mesh size");};
      };

      void setFieldData(std::shared_ptr<PanNDE::FieldBundle> solution){
        std::string keyname;
        for(int k=0;k<solution->fieldCount();k++){
          keyname=solution->fieldName(k);
          auto field=solution->field(keyname);
          auto data=extractFieldData(keyname,field);
          meshCheck(solution->mesh());
          if(PanNDE::Field::NODE==field->type()){
            if(field->size()!=solution->mesh()->nodeCount()){
              throw std::runtime_error("field-mesh mismatch: " +
                                        std::to_string(field->size()) +
                                       " vs " +
                                        std::to_string(solution->mesh()->nodeCount()));
            };
            uGridVTK->GetPointData()->AddArray(data);
          }else if(PanNDE::Field::CELL==field->type()){
            if(field->size()!=solution->mesh()->cellCount()){
              throw std::runtime_error("field-mesh mismatch: " +
                                        std::to_string(field->size()) +
                                       " vs " +
                                        std::to_string(solution->mesh()->cellCount()));
            };
            uGridVTK->GetCellData()->AddArray(data);
          }else{throw std::runtime_error("Unsupported Field Type for Write");};
        };
      };
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
      void setMesh(std::shared_ptr<PanNDE::Mesh> mesh){
        uGridVTK->Reset();
        setVTKPoints(mesh);
        setVTKCells(mesh);
      };
      void setVTKPoints(std::shared_ptr<PanNDE::Mesh> mesh){
        double coords[3];
        for(int32_t kn=0;kn<mesh->nodeCount();kn++){
          mesh->nodeCoordinate(kn,coords);
          nodeLocations->InsertNextPoint(coords);
        };
        uGridVTK->SetPoints(nodeLocations);
      };
      void setVTKCells(std::shared_ptr<PanNDE::Mesh> mesh){
        int32_t cell_nodes[8];
        for(int32_t kc=0;kc<mesh->cellCount();kc++){
          mesh->cell(kc,cell_nodes);
          for(int k=0;k<8;k++){
            hexCell->GetPointIds()->SetId(k,cell_nodes[k]);
          };
          uGridVTK->InsertNextCell(hexCell->GetCellType(),
                                   hexCell->GetPointIds());
        };
        uGridVTK->BuildLinks();
      };

      std::string filename;

      vtkSmartPointer<vtkPoints> nodeLocations=
                      vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkHexahedron> hexCell =
                      vtkSmartPointer<vtkHexahedron>::New();
      vtkSmartPointer<vtkUnstructuredGrid> uGridVTK =
                      vtkSmartPointer<vtkUnstructuredGrid>::New();

    private:
      vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = 
                      vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  };
};