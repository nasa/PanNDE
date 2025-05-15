/*! \headerfile VTKParallel.hpp "modules/VTKIO/include/internal/VTKParallel.hpp"
* "VTKParallel.hpp" contains the class implementation for handling parallel VTK operations
* using MPI communication.
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

#include <vtkMPIController.h>
#include <vtkMPICommunicator.h>
#include <vtkProcessGroup.h>
#include <vtkXMLPDataObjectWriter.h>

#include "Mesh.hpp"
#include "Field.hpp"
#include "Array.hpp"


namespace VTKIO {
  /*! \class VTKParallel VTKParallel.hpp "modules/VTKIO/include/internal/VTKParallel.hpp"
  *
  * Implements functionality for VTK parallel operations using MPI communication.
  * Sets up the necessary VTK parallel infrastructure based on a PanNDE communicator
  * to enable parallel writing of VTK data across multiple processes.
  *
  */
  class VTKParallel {
    public:
      /*!
      * Creates a shared pointer to a VTKParallel object.
      * \param comm std::shared_ptr<PanNDE::Communicator> The PanNDE communicator
      * \return std::shared_ptr<VTKParallel> Shared pointer to a newly created VTKParallel object
      */
      static std::shared_ptr<VTKParallel> makeShared(std::shared_ptr<PanNDE::Communicator> comm){
        auto retobj=std::make_shared<VTKParallel>(VTKParallel(comm));
        return std::move(retobj);
      };

      /*!
      * Constructor that initializes the VTK parallel environment.
      * \param communicator std::shared_ptr<PanNDE::Communicator> The PanNDE communicator
      * 
      * Sets up VTK's MPI communication infrastructure based on the PanNDE communicator.
      * Performs rank mapping between PanNDE and VTK MPI processes.
      */
      VTKParallel(std::shared_ptr<PanNDE::Communicator> communicator):comm(communicator){
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
      };

      /*!
      * Configures a parallel VTK writer for this process.
      * \param writer vtkSmartPointer<vtkXMLPDataObjectWriter> The parallel writer to configure
      * 
      * Sets up the writer with the MPI controller and appropriate piece information
      * for the current process to enable coordinated parallel writing.
      */
      void configPWriter(vtkSmartPointer<vtkXMLPDataObjectWriter> writer){
        writer->SetController(ctrlr);
        writer->SetNumberOfPieces(Nparts);
        writer->SetStartPiece(part_id);
        writer->SetEndPiece(part_id);
      };

    private:
      //! Process ID of the current process within the communicator
      int part_id;
      //! Total number of processes in the communicator
      int Nparts;
      //! PanNDE communicator used for MPI operations
      std::shared_ptr<PanNDE::Communicator> comm=nullptr;

      //! VTK MPI communicator for parallel VTK operations
      vtkSmartPointer<vtkMPICommunicator> vtkcomm =
                      vtkSmartPointer<vtkMPICommunicator>::New();
      //! Process group defining which MPI ranks participate in VTK operations
      vtkSmartPointer<vtkProcessGroup> vtkgroup =
                      vtkSmartPointer<vtkProcessGroup>::New();
      //! VTK MPI controller that coordinates parallel VTK operations
      vtkSmartPointer<vtkMPIController> ctrlr =
                      vtkSmartPointer<vtkMPIController>::New();
  };
};