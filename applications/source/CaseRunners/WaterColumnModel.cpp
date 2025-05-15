/*! WaterColumn.cpp "applications/source/WaterColumn.cpp"
* "WaterColumn.cpp" contains The core program for the Watercolumn UT simulation.
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
#include <memory>
#include <cstdint>
#include <vector>
#include <cstdarg>

#include "StandardNames.hpp"

#include "Array.hpp"
#include "Mesh.hpp"
#include "Field.hpp"
#include "MultiVariate.hpp"
#include "Model.hpp"
#include "Communicator.hpp"

#include "MPICommunicator.hpp"
#include "HostUTModel.hpp"

#include "Boss.hpp"
/*! \class Main WaterColumn.cpp "applications/source/WaterColumn.cpp"
* The main class for WaterColumnModel.cpp.   
*/
class Main {
  public:
    Main(int& argc, char**& argv){
      comm=Boss::makeSharedObject<NetMPI::MPICommunicator>();
      comm->Init(argc,argv);

      factories=Boss::HostFactoryManager::makeShared();
      file_manager=Boss::VTKIOManager::makeShared(argc,argv,factories,comm);
      partitioner_manager=Boss::METISManager::makeShared(factories,comm);
      data_manager=Boss::DataManager::makeShared(factories,comm);
      write_times=factories->makeManagedF64Array();
      xd_manager=Boss::RoundWaterExcitationManager::makeShared(file_manager,data_manager,factories,comm);
    };
    ~Main(){comm->Finalize();};

    void makeModel(/*std::string input_file*/){
      printf("[%i] commence model creation\n",comm->getProcessId());
      file_manager->openInputFile();
      partition();
      loadExcitationData();
      loadTimingData();
      setZeroInitialCondition();
      buildModel();
      printf("[%i] model creation complete\n",comm->getProcessId());
    };

    void runModel(){
      double time=0;
      double xdvals[4]={center->at(0),center->at(1),center->at(2),0};
      for(int ktw=0;ktw<write_times->size();ktw++){
        mainReport("advancing to t=%g\n",write_times->at(ktw));
        while(time<write_times->at(ktw)){
          xdvals[3]=time;
          mainReport("excitation: f(%g)=%e\n",time,moment_excitation.at(0).at(0)->evalAt(xdvals));
          time=model->solve();
          comm->broadcastValue(&time);
        };
        comm->barrier();
        mainReport("t=%e saving write index %i\n",time,ktw);
        file_manager->writeSolution(model->getStates(),ktw);
        mainReport("write index %i saved\n",ktw);
      };
      comm->barrier();
      mainReport("run complete\n");
    };

  private:
    void mainReport(const char *fmt, ...)
    {
      if(0==comm->getProcessId()){
        printf("[0] ");
        va_list argp;
        va_start(argp, fmt);
        vprintf(fmt, argp);
        va_end(argp);
      };
    };
    void buildModel(){
      mainReport("build model\n");
      auto utmodel=std::make_shared<HostSolver::HostUTModel>(
              HostSolver::HostUTModel(dt,data_manager->getDomain(),data_manager->getSolution(),comm));
      for(int k=0;k<moment_excitation.size();k++){
        utmodel->applyMomentExcitation(moment_excitation.at(k).data());
      };
      model=utmodel;
    };
    void setZeroInitialCondition(){
      mainReport("zeroing initial conditions\n");
      PanNDE::ElasticStateNames state_names;
      for(int kd=0;kd<3;kd++){data_manager->emplaceNewNodeSolutionField(state_names.V[kd]);};
      for(int ki=0;ki<3;ki++){
        for(int kj=ki;kj<3;kj++){
          data_manager->emplaceNewCellSolutionField(state_names.S[ki][kj]);
        };
      };
    };
    void partition(){
      partitionAndGetLocalMesh();
      loadElasticFields();
    };
    void partitionAndGetLocalMesh(){
      auto gmesh=file_manager->readMesh();
      if(nullptr!=gmesh){presentBoundingBox(gmesh);};
      auto lmesh=partitioner_manager->partitionAndDistribute(gmesh);
      data_manager->assignMesh(lmesh);
    };
    void loadElasticFields(){
      mainReport("load fields\n");
      PanNDE::ElasticMaterialNames property_names;
      distributeAndEmplaceDomainField(property_names.density);
      for(int ki=0;ki<6;ki++){
        for(int kj=ki;kj<6;kj++){
          distributeAndEmplaceDomainField(property_names.CIJ[ki][kj]);
        };
      };
    };
    void distributeAndEmplaceDomainField(std::string field_name){
      auto gfield=file_manager->getField(field_name);
      auto lfield=partitioner_manager->distributeFieldPartitions(gfield);
      data_manager->emplaceDomainField(field_name,lfield);
    };
    void loadExcitationData(){
      mainReport("load Excitation Data\n");
      int Nxd=xd_manager->getExcitationCount();
      moment_excitation.reserve(Nxd);
      printf("[%i] Loading %i Excitation:\n",comm->getProcessId(),Nxd);
      for(int kxd=0;kxd<Nxd;kxd++){
        moment_excitation.push_back({xd_manager->extractXXExcitation(kxd),xd_manager->extractYYExcitation(kxd),xd_manager->extractZZExcitation(kxd),nullptr,nullptr,nullptr});
      };
      center=xd_manager->getExcitationCenter(0);
    };
    void loadTimingData(){
      mainReport("load timing data\n");
      PanNDE::TimeDomainMetadataNames timing_names;
      dt=file_manager->getValue(timing_names.dt);
      write_times=file_manager->getArray(timing_names.write_times);
      mainReport("  dt=%e\n",dt);
    };

    void presentBoundingBox(std::shared_ptr<PanNDE::Mesh> amesh){
      double aabb[6]={0.,0.,0.,0.,0.,0.};
      double pt[3];
      amesh->nodeCoordinate(0,pt);
      for(int kp=0;kp<3;kp++){aabb[kp]=pt[kp];aabb[kp+3]=pt[kp];};
      for(int k=1;k<amesh->nodeCount();k++){
        amesh->nodeCoordinate(k,pt);
        for(int kp=0;kp<3;kp++){
          aabb[kp]=std::min(pt[kp],aabb[kp]);
          aabb[kp+3]=std::max(pt[kp],aabb[kp+3]);
        };
      };
      printf("[%i] bounding box:\n  %e %e %e\n  %e %e %e\n",
              comm->getProcessId(),aabb[0],aabb[1],aabb[2],
              aabb[3],aabb[4],aabb[5]);
    };

    std::shared_ptr<PanNDE::Communicator> comm=nullptr;
    std::shared_ptr<Boss::VTKIOManager> file_manager=nullptr;
    std::shared_ptr<Boss::HostFactoryManager> factories=nullptr;
    std::shared_ptr<Boss::METISManager> partitioner_manager=nullptr;
    std::shared_ptr<Boss::DataManager> data_manager=nullptr;
    std::shared_ptr<Boss::RoundWaterExcitationManager> xd_manager=nullptr;

    std::shared_ptr<PanNDE::Array<double>> center=nullptr;
    std::vector<std::array<std::shared_ptr<PanNDE::MultiVariate>,6>> moment_excitation;

    double dt;
    std::shared_ptr<PanNDE::Array<double>> write_times=nullptr;
    std::shared_ptr<PanNDE::Model> model=nullptr;
};
/*!
* Executes the model
*/
int main(int argc, char** argv){
  Main driver(argc,argv);
  driver.makeModel();
  driver.runModel();
  return 0;
};