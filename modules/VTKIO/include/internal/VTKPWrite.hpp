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
#include <vtkXMLPUnstructuredGridWriter.h>
//#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkMPIController.h>
#include <vtkMPICommunicator.h>
#include <vtkProcessGroup.h>

#include "Mesh.hpp"
#include "Field.hpp"
#include "Array.hpp"
#include "Communicator.hpp"

#include "VTKWrite.hpp"

namespace VTKIO {
  /*
  Only Intended for a single partition per rank at this time
  */
  class VTKPWrite : public VTKIO::VTKWrite {
    public:
      VTKPWrite(std::shared_ptr<PanNDE::Communicator> communicator){
        comm=communicator;
        Nparts=comm->getNumberOfProcesses();
        part_id=comm->getProcessId();
        int gid=comm->getGlobalId();
        std::vector<int32_t> granks;granks.resize(Nparts);
        granks.at(part_id)=gid;
        //stupid, hacky allgather since I don't have the factory handy
        for(int k=0;k<Nparts;k++){comm->broadcastValue(&granks.at(k),k);};

        vtkgroup->Initialize(vtkcomm->GetWorldCommunicator());
        vtkgroup->RemoveAllProcessIds();
        for(auto it=granks.begin();it!=granks.end();it++){vtkgroup->AddProcessId(*it);};
        vtkcomm->Initialize(vtkgroup);
        //ctrlr->Initialize(nullptr,nullptr,1);
        ctrlr->SetCommunicator(vtkcomm);
        writer->SetController(ctrlr);
      };


    protected:
      void commitToFile()override{
        writer->SetFileName(filename.c_str());
        writer->SetInputData(uGridVTK);
        writer->SetNumberOfPieces(Nparts);
        writer->SetStartPiece(part_id);
        writer->SetEndPiece(part_id);
        writer->Update();
        writer->Write();
      };
      void setFileName(std::string filename_base,int write_index)override{
        filename=filename_base + std::to_string(write_index) + ".pvtu";
      };
      void setFileName(std::string filename_base)override{
        filename=filename_base + ".pvtu";
      };

    private:
      int part_id;
      int Nparts;
      std::shared_ptr<PanNDE::Communicator> comm=nullptr;

      vtkSmartPointer<vtkMPICommunicator> vtkcomm =
                      vtkSmartPointer<vtkMPICommunicator>::New();
      vtkSmartPointer<vtkProcessGroup> vtkgroup =
                      vtkSmartPointer<vtkProcessGroup>::New();
      vtkSmartPointer<vtkMPIController> ctrlr =
                      vtkSmartPointer<vtkMPIController>::New();
      vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer =
                      vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
  };
};