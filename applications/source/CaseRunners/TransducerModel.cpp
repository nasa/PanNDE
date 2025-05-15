//Transducer.cpp

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
#include "HostRSGTDUT.hpp"

#include "Boss.hpp"

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
      xd_manager=Boss::ArbSigTransducerManager::makeShared(file_manager,data_manager,factories,comm);
      transducer_node_ids=factories->makeManagedI64Array();
      node_coords_bundle=factories->makeManagedF64DataBundle();
      node_coords_bundle->emplaceArray("node_x",factories->makeManagedF64Array());
      node_coords_bundle->emplaceArray("node_y",factories->makeManagedF64Array());
      node_coords_bundle->emplaceArray("node_z",factories->makeManagedF64Array());
      velocities_bundle=factories->makeManagedF64DataBundle();
      velocities_bundle->emplaceArray("Vx",factories->makeManagedF64Array());
      velocities_bundle->emplaceArray("Vy",factories->makeManagedF64Array());
      velocities_bundle->emplaceArray("Vz",factories->makeManagedF64Array());
    };
    ~Main(){comm->Finalize();};

    void makeModel(/*std::string input_file*/){
      printf("[%i] commence model creation\n",comm->getProcessId());
      file_manager->openInputFile();
      partition();
      loadTransducerData();
      loadTimingData();
      setZeroInitialCondition();
      buildModel();
      setupTransducerNodes();
      printf("[%i] model creation complete\n",comm->getProcessId());
      printf("[%i] Print Transducer Nodes size %lli \n",comm->getProcessId(), transducer_node_ids->size());
      velocities_bundle->array("Vx")->resize(transducer_node_ids->size());
      velocities_bundle->array("Vy")->resize(transducer_node_ids->size());
      velocities_bundle->array("Vz")->resize(transducer_node_ids->size());
      file_manager->writeTable("node_coords", node_coords_bundle); 
    };

    void runModel(){
      double time=0;
      double xdvals[4]={center->at(0),center->at(1),center->at(2),0};
      for(int ktw=0;ktw<write_times->size();ktw++){
        mainReport("advancing to t=%g\n",write_times->at(ktw));
        while(time<write_times->at(ktw)){
          xdvals[3]=time;
          mainReport("excitation: f(%g)=%e\n",time,xducers.at(2)->evalAt(xdvals));
          time=model->solve();
          comm->broadcastValue(&time);
        };
        comm->barrier();
        updateVelocityBundle(model->getStates());
        mainReport("t=%e saving write index %i\n",time,ktw);
        file_manager->writeSolution(model->getStates(),ktw);
        file_manager->writeTable("velocities", velocities_bundle, ktw);
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
    void writeVelocityTables(int ktw){

    };
    void buildModel(){
      mainReport("build model\n");
      auto xds=(std::shared_ptr<PanNDE::MultiVariate>(*)[3])xducers.data();
      model=std::make_shared<HostSolver::HostRSGTDUT>(
              HostSolver::HostRSGTDUT(dt,data_manager->getDomain(),data_manager->getSolution(),
                                      xds,xducers.size()/3, comm));
    };
    void setupTransducerNodes(){
      printf("[%i] Setting up Transducer Nodes\n",comm->getProcessId());
      int64_t nodeCount = 0;
      
      auto amesh=data_manager->getMesh();
      for(int k=0;k<amesh->nodeCount();k++){
        if(checkXds(amesh,k)){nodeCount++;};
      };
      printf("[%i] number of transducer nodes: %lli\n",comm->getProcessId(),nodeCount);
      transducer_node_ids->resize(nodeCount);
      node_coords_bundle->array("node_x")->resize(nodeCount);
      node_coords_bundle->array("node_y")->resize(nodeCount);
      node_coords_bundle->array("node_z")->resize(nodeCount);
      int64_t nk = 0;
      for(int k=0;k<amesh->nodeCount();k++){
        if (checkXds(amesh,k)){
          transducer_node_ids->at(nk)=k;
          double coord[3];
          amesh->nodeCoordinate(k,coord);
          node_coords_bundle->array("node_x")->at(nk)=coord[0];
          node_coords_bundle->array("node_y")->at(nk)=coord[1];
          node_coords_bundle->array("node_z")->at(nk)=coord[2];
          //printf("[%i] Transducer Node %lli at (%g, %g, %g) \n",comm->getProcessId(), nk, coord[0],coord[1], coord[2]);
          nk++;
        };
      };
    };
    bool checkXds(std::shared_ptr<PanNDE::Mesh> amesh, int n_idx){
      bool check=false;
      double coord[3];
      amesh->nodeCoordinate(n_idx,coord);
      for(int k=0;k<xducers.size();k++){
        check=((nullptr==xducers.at(k))?0.:(xducers.at(k)->checkAt(coord)));
        if(check){return check;};
      };
      return check;
    };
    void updateVelocityBundle(std::shared_ptr<PanNDE::FieldBundle> solution=nullptr){
      for(int k=0;k<transducer_node_ids->size();k++){
          velocities_bundle->array("Vx")->at(k)=solution->field(state_names.V[0])->at(transducer_node_ids->at(k));
          velocities_bundle->array("Vy")->at(k)=solution->field(state_names.V[1])->at(transducer_node_ids->at(k));
          velocities_bundle->array("Vz")->at(k)=solution->field(state_names.V[2])->at(transducer_node_ids->at(k));
        }
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
    void loadTransducerData(){
      mainReport("load transducers\n");
      int Nxd=xd_manager->getTransducerCount();
      xducers.reserve(3*Nxd);
      printf("[%i] Loading %i transducers:\n",comm->getProcessId(),Nxd);
      for(int kxd=0;kxd<Nxd;kxd++){
        printf("[%i] Loading transducer %i of %i:\n",comm->getProcessId(),kxd,Nxd);
        xducers.push_back(nullptr);
        xducers.push_back(nullptr);
        xducers.push_back(xd_manager->extractTransducer(kxd));
      };
      printf("[%i] Transducer array size %i :\n",comm->getProcessId(),3*Nxd);
      center=xd_manager->getTransducerCenter(0);
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
    std::shared_ptr<Boss::ArbSigTransducerManager> xd_manager=nullptr;

    std::shared_ptr<PanNDE::Array<double>> center=nullptr;
    std::vector<std::shared_ptr<PanNDE::MultiVariate>> xducers;

    double dt;
    std::shared_ptr<PanNDE::Array<double>> write_times=nullptr;
    std::shared_ptr<PanNDE::Model> model=nullptr;

    std::shared_ptr<PanNDE::Array<int64_t>> transducer_node_ids=nullptr;
    std::shared_ptr<PanNDE::DataBundle<double>> node_coords_bundle=nullptr;
    std::shared_ptr<PanNDE::DataBundle<double>> velocities_bundle=nullptr;
    PanNDE::ElasticStateNames state_names;

};

int main(int argc, char** argv){
  Main driver(argc,argv);
  driver.makeModel();
  driver.runModel();
  return 0;
};