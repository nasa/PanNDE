/*! \headerfile VTKSWrite.hpp "modules/VTKIO/include/internal/VTKSWrite.hpp"
* "VTKSWrite.hpp" contains the class implementation for serial (non-parallel) writing
* of PanNDE data to VTK files.
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
  /*! \class VTKSWrite VTKSWrite.hpp "modules/VTKIO/include/internal/VTKSWrite.hpp"
  *
  * Concrete implementation of VTKWrite for serial (non-parallel) writing of VTK files.
  * Provides standard XML writers for both grid and table data without MPI parallelism.
  * Used for single-process or serial output of PanNDE data to VTK format.
  *
  */
  class VTKSWrite : public VTKWrite{
    public:
      /*!
      * Creates a shared pointer to a VTKSWrite object.
      * \return std::shared_ptr<VTKSWrite> Shared pointer to a newly created VTKSWrite object
      */
      static std::shared_ptr<VTKSWrite> makeShared(){
        auto retobj=std::make_shared<VTKSWrite>(VTKSWrite());
        return std::move(retobj);
      };
    private:
      /*!
      * Gets the VTK table writer for serial output.
      * \return vtkSmartPointer<vtkXMLWriter> The table writer
      */
      vtkSmartPointer<vtkXMLWriter> getTableWriter()override{return table_writer;};
      
      /*!
      * Gets the VTK grid writer for serial output.
      * \return vtkSmartPointer<vtkXMLWriter> The grid writer
      */
      vtkSmartPointer<vtkXMLWriter> getGridWriter()override{return grid_writer;};

      //! Writer for unstructured grid (mesh and field) data
      vtkSmartPointer<vtkXMLUnstructuredGridWriter> grid_writer=
                      vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
      
      //! Writer for tabular data
      vtkSmartPointer<vtkXMLTableWriter> table_writer=
                      vtkSmartPointer<vtkXMLTableWriter>::New();
  };
};