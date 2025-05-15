/*! \headerfile VTKMeta.hpp "modules/VTKIO/include/internal/VTKMeta.hpp"
* "VTKMeta.hpp" contains the class implementation for handling metadata operations
* with VTK data objects.
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
#include <vtkDoubleArray.h>

#include <vtkDataObject.h>

#include "Mesh.hpp"
#include "Field.hpp"
#include "Array.hpp"


namespace VTKIO {
  /*! \class VTKMeta VTKMeta.hpp "modules/VTKIO/include/internal/VTKMeta.hpp"
  *
  * Implements functionality to transfer scalar and array metadata from
  * PanNDE DataBundle objects to VTK data objects for visualization and processing.
  *
  */
  class VTKMeta {
    public:
      /*!
      * Creates a shared pointer to a VTKMeta object.
      * \return std::shared_ptr<VTKMeta> Shared pointer to a newly created VTKMeta object
      */
      static std::shared_ptr<VTKMeta> makeShared(){
        auto newobj=std::make_shared<VTKMeta>(VTKMeta());
        return std::move(newobj);
      };

      /*!
      * Sets metadata from a PanNDE DataBundle to a VTK data object.
      * \param dataobj vtkSmartPointer<vtkDataObject> The VTK object to receive metadata
      * \param meta std::shared_ptr<PanNDE::DataBundle<double>> The source metadata
      */
      void setMetaData(vtkSmartPointer<vtkDataObject> dataobj,
                       std::shared_ptr<PanNDE::DataBundle<double>> meta){
        printf("set scalar metadata\n");
        setScalars(dataobj,meta);
        printf("set array metadata\n");
        setArrays(dataobj,meta);
        printf("meta assigned\n");
      };

    private:
      /*!
      * Sets scalar values from the DataBundle to the VTK data object.
      * \param dataobj vtkSmartPointer<vtkDataObject> The VTK object to receive scalar data
      * \param meta std::shared_ptr<PanNDE::DataBundle<double>> The source metadata
      */
      void setScalars(vtkSmartPointer<vtkDataObject> dataobj,
                      std::shared_ptr<PanNDE::DataBundle<double>> meta){
        double value;
        for(int k=0;k<meta->scalarCount();k++){
          auto keyname=meta->scalarName(k);
          value=meta->scalar(keyname);
          set1DArray(dataobj,&value,1,keyname);
        };
      };

      /*!
      * Sets array values from the DataBundle to the VTK data object.
      * \param dataobj vtkSmartPointer<vtkDataObject> The VTK object to receive array data
      * \param meta std::shared_ptr<PanNDE::DataBundle<double>> The source metadata
      */
      void setArrays(vtkSmartPointer<vtkDataObject> dataobj,
                     std::shared_ptr<PanNDE::DataBundle<double>> meta){
        for(int k=0;k<meta->arrayCount();k++){
          auto keyname=meta->arrayName(k);
          auto array=meta->array(keyname);
          set1DArray(dataobj,array->data(),array->size(),keyname);
        };
      };

      /*!
      * Creates and adds a one-dimensional array to a VTK data object's field data.
      * \param dataobj vtkSmartPointer<vtkDataObject> The VTK object to receive the array
      * \param array double* Pointer to the source data array
      * \param N int The number of elements in the array
      * \param name std::string The name identifier for the array
      */
      void set1DArray(vtkSmartPointer<vtkDataObject> dataobj,
                      double* array,int N,std::string name){
        vtkSmartPointer<vtkDoubleArray> series = vtkSmartPointer<vtkDoubleArray>::New();
        series->SetNumberOfComponents(1);
        series->SetNumberOfTuples(N);
        double value;
        for(int k=0;k<N;k++){
          value=array[k];
          series->SetTuple(k,&value);
        }
        series->SetName(name.c_str());
        dataobj->GetFieldData()->AddArray(series);
      };
  };
};