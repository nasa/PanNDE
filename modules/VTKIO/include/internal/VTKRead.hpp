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
#include <vtkCell.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

#include "Array.hpp"
#include "Mesh.hpp"
#include "Field.hpp"

namespace VTKIO {
  class VTKRead {
    public:
      void open(std::string filename){
        current_file=filename;
        loadFile();
      };
      std::shared_ptr<PanNDE::Array<double>> getArray(std::string keyname,
                              std::shared_ptr<PanNDE::ArrayFactory<double>> maker){
        ifNoFileError();
        return std::move(get1DArray(keyname,maker));
      };
      double getValue(std::string keyname){
        ifNoFileError();
        auto vec=get1Dvector(keyname);
        if(1!=vec.size()){throw std::runtime_error("not scalar quantity");};
        return vec.at(0);
      };

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

      std::string currentlyActiveFile(){return current_file;};

    private:
      void ifNoFileError(){if(current_file.empty()){throw std::runtime_error("no file provided");};};
      void ifNoMeshError(){if(nullptr==mesh){throw std::runtime_error("no mesh loaded");};};
      void throwNoName(){throw std::runtime_error("no data by this name present");};
      void loadFile(){
        vtkreader->SetFileName(current_file.c_str());
        vtkreader->Update();
        ugrid=vtkreader->GetOutput();
        mesh=nullptr;
      };
      std::shared_ptr<PanNDE::Array<double>> get1DArray(std::string name,
                                     std::shared_ptr<PanNDE::ArrayFactory<double>> maker){
        auto array=maker->makeManagedArray();
        auto adata=get1Dvector(name);
        array->resize(adata.size());
        for(int k=0;k<array->size();k++){array->at(k)=adata.at(k);};
        return std::move(array);
      };
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

      std::shared_ptr<PanNDE::Mesh> mesh=nullptr;

      std::string current_file;
      vtkSmartPointer<vtkXMLUnstructuredGridReader> vtkreader=
                      vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
      vtkSmartPointer<vtkUnstructuredGrid> ugrid=
                      vtkSmartPointer<vtkUnstructuredGrid>::New();
  };
};