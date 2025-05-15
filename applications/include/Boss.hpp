/*! \headerfile Boss.hpp "applications/include/Boss.hpp"
* "Boss.hpp" contains class definitions for creating UT simulations
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
#include <getopt.h>
#include <cstdint>
#include <cmath>

#include "StandardNames.hpp"

#include "Array.hpp"
#include "Field.hpp"
#include "Mesh.hpp"
#include "Univariate.hpp"
#include "MultiVariate.hpp"
#include "Communicator.hpp"
#include "Partitioner.hpp"
#include "Gateway.hpp"

#include "HostArray.hpp"
#include "HostField.hpp"
#include "HexMesh.hpp"
#include "HostUnivariate.hpp"
#include "HostNormalXDZ.hpp"
#include "MetisPartitioner.hpp"
#include "VTKGateway.hpp"

namespace Boss {
  template<class T>
  inline std::shared_ptr<T> makeSharedObject(){
    return std::move(std::make_shared<T>(T()));
  };
  /*! \class HostFactoryManager Boss.hpp "applications/include/Boss.hpp"
  * Factories for shared pointers to different arrays and data bundles necessary for TransduverModel.cpp and WaterColumnModel.cpp 
  */
  class HostFactoryManager {
    public:
      /*!
      * Creates factories for shared pointers to different arrays and data bundles. 
      */
      HostFactoryManager(){
        mesh_maker=makeSharedObject<HostData::HexMeshFactory>();
        field_maker=makeSharedObject<HostData::HostFieldFactory>();
        dbl_maker=makeSharedObject<HostData::HostArrayFactory<double>>();
        i32_maker=makeSharedObject<HostData::HostArrayFactory<int32_t>>();
        i64_maker=makeSharedObject<HostData::HostArrayFactory<int64_t>>();
      };
      static
      /*!
      * Creates the data structures factories. 
      */
      std::shared_ptr<HostFactoryManager> makeShared(){
        auto factories=std::make_shared<HostFactoryManager>(HostFactoryManager());
        return std::move(factories);
      };
      /*!
      * Creates a field bundle.
      */
      std::shared_ptr<PanNDE::FieldBundle> makeEmptyManagedFieldBundle(){
        return field_maker->makeEmptyManagedFieldBundle();
      };
      /*!
      * Creates a float64 bundle.
      */
      std::shared_ptr<PanNDE::DataBundle<double>> makeManagedF64DataBundle(){
        return dbl_maker->makeManagedDataBundle();
      };
      /*!
      * Creates an int32 bundle.
      */
      std::shared_ptr<PanNDE::DataBundle<int32_t>> makeManagedI32DataBundle(){
        return i32_maker->makeManagedDataBundle();
      };
      /*!
      * Creates an int64 bundle.
      */
      std::shared_ptr<PanNDE::DataBundle<int64_t>> makeManagedI64DataBundle(){
        return i64_maker->makeManagedDataBundle();
      };
      /*!
      * Creates a float64 array.
      */
      std::shared_ptr<PanNDE::Array<double>> makeManagedF64Array(){
        return dbl_maker->makeManagedArray();
      };
      /*!
      * Creates an int32 array.
      */
      std::shared_ptr<PanNDE::Array<int32_t>> makeManagedI32Array(){
        return i32_maker->makeManagedArray();
      };
      /*!
      * Creates an int64 array.
      */
      std::shared_ptr<PanNDE::Array<int64_t>> makeManagedI64Array(){
        return i64_maker->makeManagedArray();
      };
      /*!
      * Creates a field of cells.
      */
      std::shared_ptr<PanNDE::Field> makeManagedCellField(std::shared_ptr<PanNDE::Mesh> mesh){
        return field_maker->makeManagedField(mesh,PanNDE::Field::CELL);
      };
      /*!
      * Creates a field of nodes.
      */
      std::shared_ptr<PanNDE::Field> makeManagedNodeField(std::shared_ptr<PanNDE::Mesh> mesh){
        return field_maker->makeManagedField(mesh,PanNDE::Field::NODE);
      };
      /*!
      * Creates a mesh factory.
      */
      std::shared_ptr<PanNDE::MeshFactory> getMeshFactory(){return mesh_maker;};
      /*!
      * Creates a field factory.
      */
      std::shared_ptr<PanNDE::FieldFactory> getFieldFactory(){return field_maker;};
      /*!
      * Creates a float64 factory.
      */
      std::shared_ptr<PanNDE::ArrayFactory<double>> getF64ArrayFactory(){return dbl_maker;};
      /*!
      * Creates an int32 factory.
      */
      std::shared_ptr<PanNDE::ArrayFactory<int32_t>> getI32ArrayFactory(){return i32_maker;};
      /*!
      * Creates an int64 factory.
      */
      std::shared_ptr<PanNDE::ArrayFactory<int64_t>> getI64ArrayFactory(){return i64_maker;};

    private:
      std::shared_ptr<PanNDE::MeshFactory> mesh_maker=nullptr; /**< The mesh factory */
      std::shared_ptr<PanNDE::FieldFactory> field_maker=nullptr; /**< The field factory */
      std::shared_ptr<PanNDE::ArrayFactory<double>> dbl_maker=nullptr; /**< The factory for an array of float64 */
      std::shared_ptr<PanNDE::ArrayFactory<int32_t>> i32_maker=nullptr; /**< The factory for an array of int32 */
      std::shared_ptr<PanNDE::ArrayFactory<int64_t>> i64_maker=nullptr; /**< The factory for an array of int64 */
  };
  /*! \class VTKIOManager Boss.hpp "applications/include/Boss.hpp"
  * Managing the VTK data reading and writing necessary for TransduverModel.cpp and WaterColumnModel.cpp 
  */
  class VTKIOManager {
    public:
      /*!
      * Creates the VTKIO Manager.
      */
      VTKIOManager(int& argc, char**& argv,
                   std::shared_ptr<HostFactoryManager> factories,
                   std::shared_ptr<PanNDE::Communicator> comm){
        input_file=ParseInputFile(argc,argv,comm->getProcessId());
        output_base=ParseOutputBase(argc,argv,comm->getProcessId());
        this->comm=comm;
        gateway=VTKIO::VTKGateway::makeShared(comm);
        this->factories=factories;
      };
      /*!
      * VTKIO Manager destructor.
      */
      ~VTKIOManager(){gateway->close();};

      static
      /*!
      * Creates shared pointer to the VTKIO Manager.
      */
      std::shared_ptr<VTKIOManager> makeShared(int& argc, char**& argv,
                                               std::shared_ptr<HostFactoryManager> factories,
                                               std::shared_ptr<PanNDE::Communicator> comm){
        auto file_manager=std::make_shared<VTKIOManager>(VTKIOManager(argc,argv,factories,comm));
        return std::move(file_manager);
      };
      /*!
      * Uses the gateway to open an input file. 
      */
      void openInputFile(){
        comm->barrier();
        gateway->open(input_file);
        return;
      };
      /*!
      * Gets a shared pointer to the mesh. 
      */
      std::shared_ptr<PanNDE::Mesh> readMesh(){return gateway->getMesh(factories->getMeshFactory());};
      /*!
      * Gets a shared pointer to a field.  
      */
      std::shared_ptr<PanNDE::Field> getField(std::string field_name){
        return gateway->getField(field_name,factories->getFieldFactory());
      };
      /*!
      * Gets a value from the gateway.   
      */
      double getValue(std::string value_name){return gateway->getValue(value_name);};
      /*!
      * Gets a shared pointer to an array.  
      */
      std::shared_ptr<PanNDE::Array<double>> getArray(std::string array_name){
        return gateway->getArray(array_name,factories->getF64ArrayFactory());
      };
      /*!
      * Uses the gateway to write the full solution. 
      */
      void writeSolution(std::shared_ptr<PanNDE::FieldBundle> state,int tidx){
        gateway->writeSolution(output_base,state,tidx);
      };
      /*!
      * Uses the gateway to write a table.  
      */
      void writeTable(std::string fname_no_ext,
                      std::shared_ptr<PanNDE::DataBundle<double>> data){
        gateway->writeTable(fname_no_ext,data);
      };
      /*!
      * Uses the gateway to write a table at a specific time index.  
      */
      void writeTable(std::string fname_no_ext,
                      std::shared_ptr<PanNDE::DataBundle<double>> data,
                      int tid){
        gateway->writeTable(fname_no_ext,data,tid);
      };
      /*!
      * Uses the gateway to write a table with meta data.  
      */
      void writeTable(std::string fname_no_ext,
                      std::shared_ptr<PanNDE::DataBundle<double>> data,
                      std::shared_ptr<PanNDE::DataBundle<double>> meta){
        gateway->writeTable(fname_no_ext,data,meta);
      };
      /*!
      * Uses the gateway to write a table with meta data at a time index.  
      */
      void writeTable(std::string fname_no_ext,
                      std::shared_ptr<PanNDE::DataBundle<double>> data,
                      std::shared_ptr<PanNDE::DataBundle<double>> meta,int tid){
        gateway->writeTable(fname_no_ext,data,meta,tid);
      };
    private:
    /*!
    * Parse the input file.   
    */
    std::string ParseInputFile(int& argc, char**& argv,int rank=0){
      int c=getopt(argc, argv, "f:");
      if(-1==c){
        if(0==rank){printf("required -f <filename> argument\n");};
        exit(0);
      };
      return std::string(optarg);
    };
    /*!
    * Parse the outputfiles base name.   
    */
    std::string ParseOutputBase(int& argc, char**& argv,int rank=0){
      int c=getopt(argc, argv, "o:");
      if(-1==c){
        if(0==rank){printf("required -o <output_base> argument\n");};
        exit(0);
      };
      return std::string(optarg);
    };
      std::shared_ptr<PanNDE::Communicator> comm=nullptr; /**< The mpi communicator */
      std::shared_ptr<Controller::Gateway> gateway=nullptr; /**< The IO gateway */
      std::shared_ptr<HostFactoryManager> factories=nullptr; /**< The factories */
      std::string input_file; /**< The input filename */
      std::string output_base; /**< The ouput filenames' base */
  };
  /*! \class METISManager Boss.hpp "applications/include/Boss.hpp"
  * Managing the METIS functions for partitioning necessary for TransduverModel.cpp and WaterColumnModel.cpp 
  */
  class METISManager {
    public:
      /*!
      * Constructor for the METIS manger. 
      */
      METISManager(std::shared_ptr<HostFactoryManager> factories,
                   std::shared_ptr<PanNDE::Communicator> comm){
        partitioner=std::make_shared<NetMPI::MetisPartitioner>(
                      NetMPI::MetisPartitioner(factories->getI32ArrayFactory(),
                                               factories->getI64ArrayFactory(),
                                               factories->getF64ArrayFactory(),
                                               comm));
        this->comm=comm;
        this->factories=factories;
        getPIDasArray();
      };
      static
      /*!
      * Creates ashared pointer for the METIS Manager.
      */ 
      std::shared_ptr<METISManager> makeShared(std::shared_ptr<HostFactoryManager> factories,
                                               std::shared_ptr<PanNDE::Communicator> comm){
        auto partitioner=std::make_shared<METISManager>(METISManager(factories,comm));
        return std::move(partitioner);
      };
      /*!
      * Wrapper for the Internal partitioning and distributing of the mesh.
      */ 
      std::shared_ptr<PanNDE::Mesh> partitionAndDistribute(std::shared_ptr<PanNDE::Mesh> global_mesh){
        if(comm->getNumberOfProcesses()==1){
          return global_mesh;
        };
        auto lmesh=meshPartitionDistribution(global_mesh);
        return lmesh;
      };
      /*!
      * Wrapper for the Internal partitioned fields.
      */ 
      std::shared_ptr<PanNDE::Field> distributeFieldPartitions(std::shared_ptr<PanNDE::Field> global_field){
        if(comm->getNumberOfProcesses()==1){
          return global_field;
        };
        auto lfield=fieldDistribution(global_field);
        return lfield;
      };
    private:
      /*!
      * Wrapper for the partitioner for the fields,
      */
      std::shared_ptr<PanNDE::Field> fieldDistribution(std::shared_ptr<PanNDE::Field> global_field){
        auto lfields=partitioner->distributeFieldPartitions(p_ids,factories->getFieldFactory(),global_field);
        return lfields->at(0);
      };
      /*!
      * Wrapper for the partitioner for the mesh.
      */
      std::shared_ptr<PanNDE::Mesh> meshPartitionDistribution(std::shared_ptr<PanNDE::Mesh> global_mesh){
        if(nullptr!=global_mesh){
          printf("  partitioning ...\n");
          partitioner->partitionMesh(comm->getNumberOfProcesses(),global_mesh);
        };
        printf("  distributing ...\n");
        auto lmeshes=partitioner->distributeMeshPartitions(p_ids,factories->getMeshFactory());
        return lmeshes->at(0);
      };
      /*!
      * Return the process IDs as an array. 
      */
      void getPIDasArray(){
        p_ids=factories->makeManagedI32Array();
        p_ids->resize(1);
        p_ids->at(0)=comm->getProcessId();
      };
      std::shared_ptr<PanNDE::Array<int32_t>> p_ids=nullptr; /**< The array of process IDs */
      std::shared_ptr<Controller::Partitioner> partitioner=nullptr; /**< The partitioner */
      std::shared_ptr<HostFactoryManager> factories=nullptr; /**< The factories */
      std::shared_ptr<PanNDE::Communicator> comm=nullptr; /**< The communicator */
  };
  /*! \class DataManager Boss.hpp "applications/include/Boss.hpp"
  * Performing some data management tasks necessary for TransduverModel.cpp and WaterColumnModel.cpp 
  */
  class DataManager {
    public:
      /*!
      * Constructor for the DataManager. 
      */
      DataManager(std::shared_ptr<HostFactoryManager> factories,
                  std::shared_ptr<PanNDE::Communicator> comm){
        domain=factories->makeEmptyManagedFieldBundle();
        solution=factories->makeEmptyManagedFieldBundle();
        this->factories=factories;
      };
      static
      /*!
      * Creates ashared pointer for the DataManager.
      */ 
      std::shared_ptr<DataManager> makeShared(std::shared_ptr<HostFactoryManager> factories,
                                              std::shared_ptr<PanNDE::Communicator> comm){
        auto dm=std::make_shared<DataManager>(DataManager(factories,comm));
        return std::move(dm);
      };
      /*!
      * Assigns the local mesh to the the pointers for other meshs
      */
      void assignMesh(std::shared_ptr<PanNDE::Mesh> local_mesh){
        mesh=local_mesh;
        domain->mesh()=local_mesh;
        solution->mesh()=local_mesh;
      };
      /*!
      * Returns a shared pointer to the mesh. 
      */
      std::shared_ptr<PanNDE::Mesh> getMesh(){return mesh;};
      /*!
      * Places the fields in the domain
      */
      void emplaceDomainField(std::string field_name,std::shared_ptr<PanNDE::Field> local_field){
        domain->emplaceField(field_name,local_field);
      };
      /*!
      * Places the Cells in the solution
      */
      void emplaceNewCellSolutionField(std::string field_name){
        solution->emplaceField(field_name,factories->makeManagedCellField(mesh));
      };
      /*!
      * Places the Nodes in the solution
      */
      void emplaceNewNodeSolutionField(std::string field_name){
        solution->emplaceField(field_name,factories->makeManagedNodeField(mesh));
      };
      /*!
      * Returns a shared pointer to the domain. 
      */
      std::shared_ptr<PanNDE::FieldBundle> getDomain(){return domain;};
      /*!
      * Returns a shared pointer to the solution. 
      */
      std::shared_ptr<PanNDE::FieldBundle> getSolution(){return solution;};

    private:
      std::shared_ptr<PanNDE::Mesh> mesh=nullptr; /**< The mesh */
      std::shared_ptr<PanNDE::FieldBundle> domain=nullptr; /**< The Field bundle of the domain */
      std::shared_ptr<PanNDE::FieldBundle> solution=nullptr; /**< The Field bundle of the solution */
      std::shared_ptr<HostFactoryManager> factories=nullptr; /**< The factories */
  };
  /*! \class ArbSigTransducerManager Boss.hpp "applications/include/Boss.hpp"
  * Functions to create and manage a the transducers necessary for TransduverModel.cpp.
  */
  class ArbSigTransducerManager {
    public:
      /*!
      * Constructor for the ArbSigTransducerManager. 
      */
      ArbSigTransducerManager(std::shared_ptr<VTKIOManager> file_manager,
                              std::shared_ptr<DataManager> data_manager,
                              std::shared_ptr<HostFactoryManager> factories,
                              std::shared_ptr<PanNDE::Communicator> comm){
        this->file_manager=file_manager;
        this->factories=factories;
        this->data_manager=data_manager;
        this->comm=comm;
      };

      static
      /*!
      * Creates ashared pointer for the ArbSigTransducerManager.
      */ 
      std::shared_ptr<ArbSigTransducerManager> makeShared(std::shared_ptr<VTKIOManager> file_manager,
                                                          std::shared_ptr<DataManager> data_manager,
                                                          std::shared_ptr<HostFactoryManager> factories,
                                                          std::shared_ptr<PanNDE::Communicator> comm){
        auto xd= std::make_shared<ArbSigTransducerManager>(ArbSigTransducerManager(
                                  file_manager,data_manager,factories,comm));
        return std::move(xd);
      };
      /*!
      * Returns the number of transducers.
      */ 
      int getTransducerCount(){
        Nxd=file_manager->getValue(xdNames.transducerCount());
        return Nxd;
      };
      /*!
      * Extracts the information for a transducer from an input file
      */ 
      std::shared_ptr<PanNDE::MultiVariate> extractTransducer(int xd_idx){
        dz=getCellSize(data_manager->getMesh());
        auto center=getTransducerCenter(xd_idx);
        double radius=file_manager->getValue(xdNames.radius(xd_idx));
        auto sig_times=file_manager->getArray(xdNames.signalTimes(xd_idx));
        auto sig_values=file_manager->getArray(xdNames.signalValues(xd_idx));
        std::shared_ptr<PanNDE::Univariate> signal=
                    HostData::HostUnivariate::makeShared(sig_times,sig_values);
        std::shared_ptr<PanNDE::MultiVariate> xd=
                    Stubs::HostNormalXDZ::makeShared(center->data(),radius,dz,signal);
        return std::move(xd);
      };
      /*!
      * Returns the center for a transducer
      */ 
      std::shared_ptr<PanNDE::Array<double>> getTransducerCenter(int xd_idx){
        auto center=factories->makeManagedF64Array();
        center->resize(3);
        center->at(0)=file_manager->getValue(xdNames.centerX(xd_idx));
        center->at(1)=file_manager->getValue(xdNames.centerY(xd_idx));
        center->at(2)=file_manager->getValue(xdNames.centerZ(xd_idx));
        return std::move(center);
      };

    private:
      /*!
      * gets the size of cell 
      */ 
      double getCellSize(std::shared_ptr<PanNDE::Mesh> mesh,int cell_id=0){
        auto cell=mesh->cell(cell_id);
        if(cell->size()!=8){throw std::logic_error("invalid mesh for this transducer manager");};
        auto pt1=mesh->nodeCoordinate(cell->at(0));
        auto pt2=mesh->nodeCoordinate(cell->at(6));
        return std::abs(pt1->at(2)-pt2->at(2));
      };
      int Nxd=0; /**< number of transducers */
      double dz=0.; /**< cell size */
      PanNDE::TransducerWithArbitrarySignalNames xdNames; /**< transducer names */
      std::shared_ptr<VTKIOManager> file_manager=nullptr; /**< io manager */
      std::shared_ptr<DataManager> data_manager=nullptr; /**< data manager */
      std::shared_ptr<HostFactoryManager> factories=nullptr; /**< factories */
      std::shared_ptr<PanNDE::Communicator> comm=nullptr; /**< communicator */
  };
  /*! \class RoundWaterExcitationManager Boss.hpp "applications/include/Boss.hpp"
  * Functions to create and manage the water column excitations for WaterColumnModel.cpp 
  *
  */  
  class RoundWaterExcitationManager {
    public:
      /*!
      * Constructor for the Water Column. 
      */
      RoundWaterExcitationManager(std::shared_ptr<VTKIOManager> file_manager,
                              std::shared_ptr<DataManager> data_manager,
                              std::shared_ptr<HostFactoryManager> factories,
                              std::shared_ptr<PanNDE::Communicator> comm){
        this->file_manager=file_manager;
        this->factories=factories;
        this->data_manager=data_manager;
        this->comm=comm;
      };

      static
      /*!
      * creates a shared pointer for  Water Column. 
      */
      std::shared_ptr<RoundWaterExcitationManager> makeShared(std::shared_ptr<VTKIOManager> file_manager,
                                                          std::shared_ptr<DataManager> data_manager,
                                                          std::shared_ptr<HostFactoryManager> factories,
                                                          std::shared_ptr<PanNDE::Communicator> comm){
        auto xd= std::make_shared<RoundWaterExcitationManager>(RoundWaterExcitationManager(
                                  file_manager,data_manager,factories,comm));
        return std::move(xd);
      };
      /*!
      * returns the number of water columns
      */
      int getExcitationCount(){
        Nxd=file_manager->getValue(xdNames.excitationCount());
        return Nxd;
      };
      /*!
      * Extracts the information for a water column X excitation from an input file
      */ 
      std::shared_ptr<PanNDE::MultiVariate> extractXXExcitation(int xd_idx){
        dz=getCellSize(data_manager->getMesh());
        auto center=getExcitationCenter(xd_idx);
        double radius=file_manager->getValue(xdNames.radius(xd_idx));
        auto sig_times=file_manager->getArray(xdNames.signalTimes(xd_idx));
        auto sigxx_values=file_manager->getArray(xdNames.signalXXValues(xd_idx));
        std::shared_ptr<PanNDE::Univariate> signal=
                    HostData::HostUnivariate::makeShared(sig_times,sigxx_values);
        std::shared_ptr<PanNDE::MultiVariate> xd=
                    Stubs::HostNormalXDZ::makeShared(center->data(),radius,dz,signal);
        return std::move(xd);
      };
      /*!
      * Extracts the information for a water column Y excitation from an input file
      */ 
      std::shared_ptr<PanNDE::MultiVariate> extractYYExcitation(int xd_idx){
        dz=getCellSize(data_manager->getMesh());
        auto center=getExcitationCenter(xd_idx);
        double radius=file_manager->getValue(xdNames.radius(xd_idx));
        auto sig_times=file_manager->getArray(xdNames.signalTimes(xd_idx));
        auto sigyy_values=file_manager->getArray(xdNames.signalYYValues(xd_idx));
        std::shared_ptr<PanNDE::Univariate> signal=
                    HostData::HostUnivariate::makeShared(sig_times,sigyy_values);
        std::shared_ptr<PanNDE::MultiVariate> xd=
                    Stubs::HostNormalXDZ::makeShared(center->data(),radius,dz,signal);
        return std::move(xd);
      };
      /*!
      * Extracts the information for a water column Z excitation from an input file
      */
      std::shared_ptr<PanNDE::MultiVariate> extractZZExcitation(int xd_idx){
        dz=getCellSize(data_manager->getMesh());
        auto center=getExcitationCenter(xd_idx);
        double radius=file_manager->getValue(xdNames.radius(xd_idx));
        auto sig_times=file_manager->getArray(xdNames.signalTimes(xd_idx));
        auto sigzz_values=file_manager->getArray(xdNames.signalZZValues(xd_idx));
        std::shared_ptr<PanNDE::Univariate> signal=
                    HostData::HostUnivariate::makeShared(sig_times,sigzz_values);
        std::shared_ptr<PanNDE::MultiVariate> xd=
                    Stubs::HostNormalXDZ::makeShared(center->data(),radius,dz,signal);
        return std::move(xd);
      };
      /*!
      * Returns the center for a watercolumn
      */ 
      std::shared_ptr<PanNDE::Array<double>> getExcitationCenter(int xd_idx){
        auto center=factories->makeManagedF64Array();
        center->resize(3);
        center->at(0)=file_manager->getValue(xdNames.centerX(xd_idx));
        center->at(1)=file_manager->getValue(xdNames.centerY(xd_idx));
        center->at(2)=file_manager->getValue(xdNames.centerZ(xd_idx));
        return std::move(center);
      };

    private:
      /*!
      * gets the size of cell 
      */ 
      double getCellSize(std::shared_ptr<PanNDE::Mesh> mesh,int cell_id=0){
        auto cell=mesh->cell(cell_id);
        if(cell->size()!=8){throw std::logic_error("invalid mesh for this transducer manager");};
        auto pt1=mesh->nodeCoordinate(cell->at(0));
        auto pt2=mesh->nodeCoordinate(cell->at(6));
        return std::abs(pt1->at(2)-pt2->at(2));
      };
      int Nxd=0; /**< number of water columns */
      double dz=0.; /**< cell size */
      PanNDE::RoundWaterColumnExcitationNames xdNames; /**< transducer names */
      std::shared_ptr<VTKIOManager> file_manager=nullptr; /**< io manager */
      std::shared_ptr<DataManager> data_manager=nullptr; /**< data manager */
      std::shared_ptr<HostFactoryManager> factories=nullptr; /**< factories */
      std::shared_ptr<PanNDE::Communicator> comm=nullptr; /**< communicator */
  };
};
