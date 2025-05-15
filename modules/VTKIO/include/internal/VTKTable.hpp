/*! \headerfile VTKTable.hpp "modules/VTKIO/include/internal/VTKTable.hpp"
* "VTKTable.hpp" contains the class implementation for converting PanNDE tabular data
* to VTK table format for visualization and analysis.
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
#include <cstdio>

#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkTable.h>

#include "Array.hpp"

#include "VTKMeta.hpp"

namespace VTKIO {
  /*! \class VTKTable VTKTable.hpp "modules/VTKIO/include/internal/VTKTable.hpp"
  *
  * Implements conversion from PanNDE data bundles to VTK table format.
  * Provides functionality to organize array data as columns in a VTK table
  * and attach metadata for visualization and analysis.
  * Inherits from VTKMeta to leverage metadata handling capabilities.
  *
  */
  class VTKTable : public VTKMeta {
    public:
      /*!
      * Creates a shared pointer to a VTKTable object.
      * \return std::shared_ptr<VTKTable> Shared pointer to a newly created VTKTable object
      */
      static std::shared_ptr<VTKTable> makeShared(){
        auto newobj=std::make_shared<VTKTable>(VTKTable());
        return std::move(newobj);
      };

      /*!
      * Sets tabular data from a PanNDE DataBundle.
      * \param data std::shared_ptr<PanNDE::DataBundle<double>> Data bundle containing arrays to use as table columns
      * 
      * Validates that all arrays have the same length before creating the table.
      * Each array becomes a column in the table with the array name as the column name.
      */
      void setData(std::shared_ptr<PanNDE::DataBundle<double>> data){
        //TODO: check column lengths are all equal, throw error if not
        int Nprev;
        for(int k=0;k<data->arrayCount();k++){
          auto keyname=data->arrayName(k);
          auto array=data->array(keyname);
          auto N=array->size();
          if(k==0){Nprev=N;};
          if(N!=Nprev){
            std::string error_str="column " + keyname + 
                                  "has mismatched length. Table write aborted.\n";
            printf("%s\n",error_str.c_str());
            return;
          };
          setColumns(data);
          return;
        };
      };

      /*!
      * Sets metadata for the VTK table.
      * \param meta std::shared_ptr<PanNDE::DataBundle<double>> Metadata to assign to the table
      */
      void setMeta(std::shared_ptr<PanNDE::DataBundle<double>> meta){
        printf("assign metadata to table\n");
        tbl_meta_ops->setMetaData(tableobj,meta);
      };
      
      /*!
      * Gets the VTK table object.
      * \return vtkSmartPointer<vtkTable> The VTK table representation
      */
      vtkSmartPointer<vtkTable> getVTKObject(){return tableobj;};

    private:
      /*!
      * Sets all columns in the table from the data bundle.
      * \param data std::shared_ptr<PanNDE::DataBundle<double>> Data bundle containing the column arrays
      */
      void setColumns(std::shared_ptr<PanNDE::DataBundle<double>> data){
        for(int k=0;k<data->arrayCount();k++){
          auto keyname=data->arrayName(k);
          auto array=data->array(keyname);
          buildVTKColumn(array->data(),array->size(),keyname);
        };
      };

      /*!
      * Creates a VTK column from an array and adds it to the table.
      * \param array double* Pointer to the source data array
      * \param N int Number of elements in the array
      * \param name std::string Name for the column
      */
      void buildVTKColumn(double* array,int N,std::string name){
        vtkSmartPointer<vtkDoubleArray> series = vtkSmartPointer<vtkDoubleArray>::New();
        series->SetNumberOfComponents(1);
        series->SetNumberOfTuples(N);
        double value;
        for(int k=0;k<N;k++){
          value=array[k];
          series->SetTuple(k,&value);
        }
        series->SetName(name.c_str());
        tableobj->AddColumn(series);
      };

      //! Helper object for handling VTK metadata operations
      std::shared_ptr<VTKMeta> tbl_meta_ops=VTKMeta::makeShared();

      //! The VTK table representation
      vtkSmartPointer<vtkTable> tableobj=vtkSmartPointer<vtkTable>::New();
  };
};