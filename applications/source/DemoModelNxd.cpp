//DemoModelNxd.cpp
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
#include <cstdarg>
#include <cstdint>
#include <vector>
#include <getopt.h>

#include "StandardNames.hpp"

#include "Array.hpp"
#include "Mesh.hpp"
#include "Field.hpp"
#include "MultiVariate.hpp"
#include "HostArray.hpp"
#include "HexMesh.hpp"
#include "HostField.hpp"
#include "HostMVSmoother.hpp"

#include "HannTransducerZ.hpp"

#include "Model.hpp"
#include "HostRSGTDUT.hpp"

#include "Communicator.hpp"
#include "Partitioner.hpp"
#include "MPICommunicator.hpp"
#include "MetisPartitioner.hpp"

#include "Gateway.hpp"
#include "VTKGateway.hpp"


inline std::string ParseInputFile(int& argc, char**& argv,int rank=0){
  int c=getopt(argc, argv, "f:");
  if(-1==c){
    if(0==rank){
      printf("required -f <filename> argument\n");
    };
    exit(0);
  };
  return std::string(optarg);
};
inline std::string ParseOutputBase(int& argc, char**& argv,int rank=0){
  int c=getopt(argc, argv, "o:");
  if(-1==c){
    if(0==rank){
      printf("required -o <output_base> argument\n");
    };
    exit(0);
  };
  return std::string(optarg);
};

class Main {
  public:
    Main(int& argc, char**& argv){
      comm=makeShared<NetMPI::MPICommunicator>();
      comm->Init(argc,argv);
      input_file=ParseInputFile(argc,argv,comm->getProcessId());
      output_base=ParseOutputBase(argc,argv,comm->getProcessId());
      gateway=std::make_shared<VTKIO::VTKGateway>(VTKIO::VTKGateway(comm));
      mesh_maker=makeShared<HostData::HexMeshFactory>();
      field_maker=makeShared<HostData::HostFieldFactory>();
      dbl_maker=makeShared<HostData::HostArrayFactory<double>>();
      i32_maker=makeShared<HostData::HostArrayFactory<int32_t>>();
      i64_maker=makeShared<HostData::HostArrayFactory<int64_t>>();
      partitioner=std::make_shared<NetMPI::MetisPartitioner>(
                    NetMPI::MetisPartitioner(i32_maker,i64_maker,dbl_maker,comm));
      domain=field_maker->makeEmptyManagedFieldBundle();
      solution=field_maker->makeEmptyManagedFieldBundle();
      sim_inputs=dbl_maker->makeManagedDataBundle();
      p_ids=i32_maker->makeManagedArray();p_ids->resize(1);p_ids->at(0)=comm->getProcessId();
    };
    ~Main(){
      gateway->close();
      comm->Finalize();
    };

    void testprint(){printf("[%i] hello.\n",comm->getProcessId());};
    void makeModel(/*std::string input_file*/){
      printf("[%i] commence model creation\n",comm->getProcessId());
      comm->barrier();
      gateway->open(input_file);
      partitionAndGetLocalMesh();
      loadElasticFields();
      loadStubTransducerData();
      //makeStubTransducers();
      loadTimingData();
      setZeroInitialCondition();
      buildModel();
      if(solution->mesh()->nodeCount()!=solution->field("Vx")->size()){
        throw std::runtime_error("invalid data size: " + 
                                  std::to_string(solution->mesh()->nodeCount()) +
                                 " vs " +
                                  std::to_string(solution->field("Vx")->size()));
      };
      printf("[%i] model creation complete\n",comm->getProcessId());
    };

    void ownerCheck(){
      for(int kn=0;kn<solution->mesh()->nodeCount();kn++){
        solution->field("Vz")->at(kn)=solution->mesh()->nodeHomePartition(kn);
      };
      gateway->writeSolution(output_base,model->getStates(),0);
    };
    void runModel(){
      double time=0;
      for(int ktw=0;ktw<write_times->size();ktw++){
        if(0==comm->getProcessId()){printf("[%i] advancing to t=%g\n",comm->getProcessId(),write_times->at(ktw));};
        while(time<write_times->at(ktw)){//(0.5*dt)<(write_times->at(kt)-time)
          if(0==comm->getProcessId()){
            double xdvals[4]={center0[0],center0[1],center0[2],time};
            printf("[%i] excitation: f(%g)=%e\n",
                    comm->getProcessId(),time,
                    xducers.at(2)->evalAt(xdvals));
          };
          time=model->solve();
          comm->broadcastValue(&time);
        };
        comm->barrier();
        if(0==comm->getProcessId()){printf("[%i] t=%e saving write index %i\n",comm->getProcessId(),time,ktw);};
        gateway->writeSolution(output_base,model->getStates(),ktw);
        if(0==comm->getProcessId()){printf("[%i] write index %i saved\n",comm->getProcessId(),ktw);};
      };
      comm->barrier();
      if(0==comm->getProcessId()){printf("[0] run complete\n");};
    };



  private:
    void buildModel(){
      if(0==comm->getProcessId()){printf("[0] build model\n");};
      //auto xds=(std::shared_ptr<PanNDE::MultiVariate>(*)[xducers.size()])xducers.data();
      //std::shared_ptr<PanNDE::MultiVariate>(*xds)[3]
      auto xds=(std::shared_ptr<PanNDE::MultiVariate>(*)[3])xducers.data();

      model=std::make_shared<HostSolver::HostRSGTDUT>(
              HostSolver::HostRSGTDUT(sim_inputs->scalar(timing_names.dt),
                                      domain,solution,xds,xducers.size()/3,comm));//xducer
    };
    void setZeroInitialCondition(){
      if(0==comm->getProcessId()){printf("[0] zeroing initial conditions\n");};
      for(int kd=0;kd<3;kd++){
        solution->emplaceField(state_names.V[kd],
                      field_maker->makeManagedField(domain->mesh(),PanNDE::Field::NODE));
      };
      for(int ki=0;ki<3;ki++){
        for(int kj=ki;kj<3;kj++){
          solution->emplaceField(state_names.S[ki][kj],
                        field_maker->makeManagedField(domain->mesh(),PanNDE::Field::CELL));
        };
      };
    };
    std::shared_ptr<PanNDE::MultiVariate> makeHannStub(double center[3],double radius,
                                                       int Ncycle,double freq,double phase){
      printf("  center: %g %g %g\n",center[0],center[1],center[2]);
      printf("  radius: %g\n",radius);
      printf("  frequency,phase: %g,%g\n",freq,phase);
      printf("  Ncycle: %i\n",Ncycle);
      double pt1[3];double pt2[3];
      int32_t box[8];
      domain->mesh()->cell(0,box);
      domain->mesh()->nodeCoordinate(box[0],pt1);
      domain->mesh()->nodeCoordinate(box[6],pt2);
      double ds=pt2[2]-pt1[2];
      auto vertical_excitation=std::make_shared<Stubs::HannTransducerZ>(
                    Stubs::HannTransducerZ(center,radius,ds));
      vertical_excitation->configure(Ncycle,freq,phase);
      return std::move(vertical_excitation);
    };

    void partitionAndGetLocalMesh(){
      auto gmesh=gateway->getMesh(mesh_maker);
      if(nullptr!=gmesh){
        presentBoundingBox(gmesh);
        partitioner->partitionMesh(comm->getNumberOfProcesses(),gmesh);
      };
      auto lmeshes=partitioner->distributeMeshPartitions(p_ids,mesh_maker);
      domain->mesh()=lmeshes->at(0);
      solution->mesh()=lmeshes->at(0);
    };
    void loadElasticFields(){
      if(0==comm->getProcessId()){printf("[0] load fields\n");};
      auto gfield=gateway->getField(property_names.density,field_maker);
      auto lfields=partitioner->distributeFieldPartitions(p_ids,field_maker,gfield);
      domain->emplaceField(property_names.density,lfields->at(0));
      for(int ki=0;ki<6;ki++){
        for(int kj=ki;kj<6;kj++){
          gfield=gateway->getField(property_names.CIJ[ki][kj],field_maker);
          lfields=partitioner->distributeFieldPartitions(p_ids,field_maker,gfield);
          domain->emplaceField(property_names.CIJ[ki][kj],lfields->at(0));
        };
      };
      //if(0==comm->getProcessId()){
      //  for(int k=0;k<16;k++){presentCIJ(k);};
      //};
    };
    void loadStubTransducerData(){
      if(0==comm->getProcessId()){printf("[0] load transducer\n");};
      PanNDE::TransducerWindowedSineParameterNames xdNames;
      int Nxd=gateway->getValue(xdNames.transducerCount());
      xducers.reserve(3*Nxd);
      double freq;
      double phase;
      int Ncycle;
      double radius;
      double center[3];

      center0[0]=gateway->getValue(xdNames.centerX(0));
      center0[1]=gateway->getValue(xdNames.centerY(0));
      center0[2]=gateway->getValue(xdNames.centerZ(0));
      printf("Loading %i transducers:\n",Nxd);
      for(int kxd=0;kxd<Nxd;kxd++){
        xducers.push_back(nullptr);
        xducers.push_back(nullptr);
        center[0]=gateway->getValue(xdNames.centerX(kxd));
        center[1]=gateway->getValue(xdNames.centerY(kxd));
        center[2]=gateway->getValue(xdNames.centerZ(kxd));
        freq=gateway->getValue(xdNames.frequency(kxd));
        Ncycle=gateway->getValue(xdNames.cycleCount(kxd));
        phase=gateway->getValue(xdNames.phase(kxd));
        radius=gateway->getValue(xdNames.radius(kxd));
        xducers.push_back(makeHannStub(center,radius,Ncycle,freq,phase));
      };
    };
    void loadTimingData(){
      if(0==comm->getProcessId()){printf("[0] load timing data\n");};
      dt=gateway->getValue(timing_names.dt);
      sim_inputs->emplaceScalar(timing_names.dt,dt);
      write_times=gateway->getArray(timing_names.write_times,dbl_maker);
      sim_inputs->emplaceArray(timing_names.write_times,write_times);
      if(0==comm->getProcessId()){printf("[0]   dt=%e\n",dt);};
    };
    template<class T>
    inline std::shared_ptr<T> makeShared(){
      return std::move(std::make_shared<T>(T()));
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
    void presentCIJ(int cidx){
      char buff[100];
      sprintf(buff,"[%i] CIJ at cell %i, z=%g:\n",comm->getProcessId(),cidx,cellDepth(cidx));
      std::string output(buff);
      for(int ki=0;ki<6;ki++){
        for(int kj=0;kj<6;kj++){
          sprintf(buff,"  %e",domain->field(property_names.CIJ[ki][kj])->atCell(cidx));
          output.append(buff);
        };
        output+="\n";
      };
      printf("%s\n",output.c_str());
    };
    double cellDepth(int cidx){
      double pt[3];double z=0.;
      int box[8];
      domain->mesh()->cell(cidx,box);
      for(int kb=0;kb<8;kb++){
        domain->mesh()->nodeCoordinate(box[kb],pt);
        z+=pt[2];
      };
      z=z/8.;
      return z;
    }

    std::string input_file;
    std::string output_base;

    double dt;
    std::shared_ptr<PanNDE::Array<double>> write_times=nullptr;

    std::shared_ptr<PanNDE::Communicator> comm=nullptr;
    std::shared_ptr<Controller::Gateway> gateway=nullptr;
    std::shared_ptr<Controller::Partitioner> partitioner=nullptr;
    std::shared_ptr<PanNDE::Array<int32_t>> p_ids=nullptr;

    std::shared_ptr<PanNDE::FieldBundle> domain=nullptr;
    std::shared_ptr<PanNDE::FieldBundle> solution=nullptr;
    std::shared_ptr<PanNDE::DataBundle<double>> sim_inputs=nullptr;
    //3 multivariates required to define a single transducer
    //std::shared_ptr<PanNDE::MultiVariate> xducer[3]={nullptr,nullptr,nullptr};
    //dual xd
    double center0[3];
    double center1[3];
    //std::shared_ptr<PanNDE::MultiVariate> xducer[2][3]={{nullptr,nullptr,nullptr},
    //                                                    {nullptr,nullptr,nullptr}};
    std::vector<std::shared_ptr<PanNDE::MultiVariate>> xducers;

    std::shared_ptr<PanNDE::Model> model=nullptr;

    std::shared_ptr<PanNDE::MeshFactory> mesh_maker=nullptr;
    std::shared_ptr<PanNDE::FieldFactory> field_maker=nullptr;
    std::shared_ptr<PanNDE::ArrayFactory<double>> dbl_maker=nullptr;
    std::shared_ptr<PanNDE::ArrayFactory<int32_t>> i32_maker=nullptr;
    std::shared_ptr<PanNDE::ArrayFactory<int64_t>> i64_maker=nullptr;

    PanNDE::ElasticMaterialNames property_names;
    PanNDE::TimeDomainMetadataNames timing_names;
    PanNDE::TransducerExcitationNames xducer_names;
    PanNDE::ElasticStateNames state_names;
};

int main(int argc, char** argv){
  Main driver(argc,argv);
  driver.makeModel();
  driver.runModel();
  return 0;
};
