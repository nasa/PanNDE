/*! \headerfile MPICommunicator.hpp "modules/NetMPI/include/MPICommunicator.hpp"
* "MPICommunicator.hpp" contains the class implementing the data communication among processes. 
* This wraps MPI
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

#include <cstdio>

#include <memory>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <mpi.h>
#include <tuple>

#include "Mesh.hpp"
#include "Field.hpp"
#include "Array.hpp"
#include "Communicator.hpp"

#include "internal/MPIType.hpp"
#include "internal/MPITransceiver.hpp"
#include "internal/MPIHaloExchanger.hpp"

namespace NetMPI {
  /*! 
  * \class MPICommunicator MPICommunicator.hpp "modules/NetMPI/include/MPICommunicator.hpp"
  * "MPICommunicator.hpp" implements the data communication among processes. 
  * This wraps MPI
  *
  */
  class MPICommunicator : public PanNDE::Communicator {
    public:
      /*!
      * constructor
      * \param comm MPI communicator to setup
      */
      MPICommunicator(MPI_Comm comm=MPI_COMM_WORLD){
        if(isInitialized()){
          setupMPI(comm);
        };
      };
      /*!
      * constructor
      * \param mesh mesh for comms patterns for halo exchange
      * \param comm MPI communicator to setup
      */
      MPICommunicator(std::shared_ptr<PanNDE::Mesh> mesh,MPI_Comm comm=MPI_COMM_WORLD){
        if(isInitialized()){
          setupMPI(comm);
          determineExchangePattern(mesh);
        };
      };
      /*!
      * destructor
      */
      ~MPICommunicator(){/*Finalize();*/};

      /*! 
      * Initialize a communicator object. As many initializers (e.g. MPI, VTK, etc.) 
      * require the arguments to `int main()', they are passed through here.
      * \param argc number of arguments in argv
      * \param argv the command line arguments provided to int main()
      */
      void Init(int& argc,char**& argv)override{
        if(!isInitialized()){
          MPI_Init(&argc,&argv);
          setupMPI(MPI_COMM_WORLD);
        };
      };
      /*! 
      * Initialize a communicator object. As many initializers (e.g. MPI, VTK, etc.) 
      * require the arguments to `int main()', they are passed through here. 
      * The data exchange pattern can be determined for final configuration if the mesh is passed in
      * \param argc number of arguments in argv
      * \param argv the command line arguments provided to int main()
      * \param mesh the local mesh for simulation
      */
      void Init(int& argc,char**& argv,std::shared_ptr<PanNDE::Mesh> mesh)override{
        if(!isInitialized()){
          MPI_Init(&argc,&argv);
          setupMPI(MPI_COMM_WORLD);
        };
        determineExchangePattern(mesh);
      };
      /*!
      * get the id number of the process in the local communicator
      */
      int getProcessId()override{ifNotInitializedError();return rank;};
      /*!
      * get the id number of the process in the global communicator
      */
      int getGlobalId()override{
        ifNotInitializedError();
        int grank;
        MPI_Comm_rank(MPI_COMM_WORLD,&grank);
        return grank;
      };
      /*!
      * get the number of processes on the local communicator
      */
      int getNumberOfProcesses()override{ifNotInitializedError();return comm_size;};
      /*!
      * get the number of processes on the global communicator
      */
      int getNumberOfAllProcesses()override{
        ifNotInitializedError();
        int gcomm_size;
        MPI_Comm_size(netMPI_comm,&gcomm_size);
        return gcomm_size;
      };
      /*!
      * release commmunicator assets
      */
      void Finalize()override{
        int flag;MPI_Finalized(&flag);
        if(!flag){MPI_Finalize();};
      };

      /*!
      * use the local mesh to determine the halo exchange pattern. Not recommended approach, but
      * retained for backwards compatibility
      * \param mesh the local mesh from which the communication pattern can be determined
      */
      void determineExchangePattern(std::shared_ptr<PanNDE::Mesh> mesh)override{
        exchanger->determineExchangePattern(mesh);
      };
      /*!
      * after having determined the exchange pattern, use that pattern to configure the halo exchange for
      * the specific field. Not recommended approach, but retained for backwards compatibility
      * \param keyname field name
      * \param field the field data to be registered for halo exchange
      */
      void setupDataLinks(std::string keyname,std::shared_ptr<PanNDE::Field> field)override{
        exchanger->setupDataLinks(keyname,field);
      };
      /*!
      * Use the field bundle to determine the exchange pattern and register all fields contained for data
      * exchange.
      * \param fields the field bundle to set up for halo exchange
      */
      void setupDataLinks(std::shared_ptr<PanNDE::FieldBundle> fields)override{
        exchanger->setupDataLinks(fields);
      };
      
      /*!
      * Start asynchronous halo exchange on field by name
      * \param keyname field name
      */
      void startHaloExchange(std::string keyname)override{
        exchanger->startHaloExchange(keyname);
      };
      /*!
      * Wait for completion of halo exchange on field by name
      * \param keyname field name
      */
      void waitUntilDone(std::string keyname)override{
        exchanger->waitUntilDone(keyname);
      };

      //Send/receives: add types as appropriate. templates would be _nice_ but virtual
      //  templates are disallowed
      /*!
      * Send array to process
      * \param array the array to send
      * \param receiver process to receive the array
      */
      void sendArray(std::shared_ptr<PanNDE::Array<double>> array,int receiver)override{
        sendArray(array->data(),array->size(),receiver);
      };
      /*!
      * Send array to process
      * \param array the array to send
      * \param N the length of the array
      * \param receiver process to receive the array
      */
      void sendArray(double* array,int N,int receiver)override{
        transceiver->isendArray<double>(array,N,receiver);
      };
      /*!
      * Receive array from process
      * \param maker the array factory method required to build the array
      * \param sender process which sent the array
      */
      std::shared_ptr<PanNDE::Array<double>> recvArray(
                      std::shared_ptr<PanNDE::ArrayFactory<double>> maker,int sender)override{
        return std::move(transceiver->irecvArray<double>(maker,sender));
      };
      /*!
      * broadcast array by process
      * \param array the array to be either populated or sent
      * \param sender process which broadcasts the array
      */
      void broadcastArray(std::shared_ptr<PanNDE::Array<double>>& array,int sender=0)override{
        transceiver->bcastArray<double>(array,sender);
      };

      /*!
      * Send array to process
      * \param array the array to send
      * \param receiver process to receive the array
      */
      void sendArray(std::shared_ptr<PanNDE::Array<int32_t>> array,int receiver)override{
        sendArray(array->data(),array->size(),receiver);
      };
      /*!
      * Send array to process
      * \param array the array to send
      * \param N the length of the array
      * \param receiver process to receive the array
      */
      void sendArray(int32_t* array,int N,int receiver)override{
        transceiver->isendArray<int32_t>(array,N,receiver);
      };
      /*!
      * Receive array from process
      * \param maker the array factory method required to build the array
      * \param sender process which sent the array
      */
      std::shared_ptr<PanNDE::Array<int32_t>> recvArray(
                      std::shared_ptr<PanNDE::ArrayFactory<int32_t>> maker,int sender)override{
        return std::move(transceiver->irecvArray<int32_t>(maker,sender));};
      /*!
      * broadcast array by process
      * \param array the array to be either populated or sent
      * \param sender process which broadcasts the array
      */
      void broadcastArray(std::shared_ptr<PanNDE::Array<int32_t>>& array,int sender=0)override{
        transceiver->bcastArray<int32_t>(array,sender);
      };

      /*!
      * Send array to process
      * \param array the array to send
      * \param receiver process to receive the array
      */
      void sendArray(std::shared_ptr<PanNDE::Array<int64_t>> array,int receiver)override{
        sendArray(array->data(),array->size(),receiver);
      };
      /*!
      * Send array to process
      * \param array the array to send
      * \param N the length of the array
      * \param receiver process to receive the array
      */
      void sendArray(int64_t* array,int N,int receiver)override{
        transceiver->isendArray<int64_t>(array,N,receiver);
      };
      /*!
      * Receive array from process
      * \param maker the array factory method required to build the array
      * \param sender process which sent the array
      */
      std::shared_ptr<PanNDE::Array<int64_t>> recvArray(
                      std::shared_ptr<PanNDE::ArrayFactory<int64_t>> maker,int sender)override{
        return std::move(transceiver->irecvArray<int64_t>(maker,sender));};
      /*!
      * broadcast array by process
      * \param array the array to be either populated or sent
      * \param sender process which broadcasts the array
      */
      void broadcastArray(std::shared_ptr<PanNDE::Array<int64_t>>& array,int sender=0)override{
        transceiver->bcastArray<int64_t>(array,sender);
      };

      /*!
      * send value to process
      * \param value the value to send
      * \param receiver process to send value
      */
      void sendValue(double value,int receiver)override{
        transceiver->isendValue<double>(value,receiver);
      };
      /*!
      * receive value from process
      * \param value the address to which to write
      * \param sender process which sent the value
      */
      void recvValue(double* value,int sender)override{
        transceiver->irecvValue<double>(value,sender);
      };
      /*!
      * broadcast value to all processes
      * \param value the address to which to write (if receiving) or send (if broadcasting)
      * \param sender process which sent the value
      */
      void broadcastValue(double* value,int sender=0)override{
        transceiver->broadcastValue<double>(value,sender);
      };
      /*!
      * all broadcast value for placement in array by process index
      * \param value the value to send
      * \param maker the array factory to build the synthesized array
      */
      std::shared_ptr<PanNDE::Array<double>> allGatherValue(double value,
                      std::shared_ptr<PanNDE::ArrayFactory<double>> maker)override{
        return std::move(transceiver->allGatherValue<double>(value,maker));
      };

      /*!
      * send value to process
      * \param value the value to send
      * \param receiver process to send value
      */
      void sendValue(int32_t value,int receiver)override{
        transceiver->isendValue<int32_t>(value,receiver);
      };
      /*!
      * receive value from process
      * \param value the address to which to write
      * \param sender process which sent the value
      */
      void recvValue(int32_t* value,int sender)override{
        transceiver->irecvValue<int32_t>(value,sender);
      };
      /*!
      * broadcast value to all processes
      * \param value the address to which to write (if receiving) or send (if broadcasting)
      * \param sender process which sent the value
      */
      void broadcastValue(int32_t* value,int sender=0)override{
        transceiver->broadcastValue<int32_t>(value,sender);
      };
      /*!
      * all broadcast value for placement in array by process index
      * \param value the value to send
      * \param maker the array factory to build the synthesized array
      */
      std::shared_ptr<PanNDE::Array<int32_t>> allGatherValue(int32_t value,
                      std::shared_ptr<PanNDE::ArrayFactory<int32_t>> maker)override{
        return std::move(transceiver->allGatherValue<int32_t>(value,maker));
      };

      /*!
      * send value to process
      * \param value the value to send
      * \param receiver process to send value
      */
      void sendValue(int64_t value,int receiver)override{
        transceiver->isendValue<int64_t>(value,receiver);
      };
      /*!
      * receive value from process
      * \param value the address to which to write
      * \param sender process which sent the value
      */
      void recvValue(int64_t* value,int sender)override{
        transceiver->irecvValue<int64_t>(value,sender);
      };
      /*!
      * broadcast value to all processes
      * \param value the address to which to write (if receiving) or send (if broadcasting)
      * \param sender process which sent the value
      */
      void broadcastValue(int64_t* value,int sender=0)override{
        transceiver->broadcastValue<int64_t>(value,sender);
      };
      /*!
      * all broadcast value for placement in array by process index
      * \param value the value to send
      * \param maker the array factory to build the synthesized array
      */
      std::shared_ptr<PanNDE::Array<int64_t>> allGatherValue(int64_t value,
                      std::shared_ptr<PanNDE::ArrayFactory<int64_t>> maker)override{
        return std::move(transceiver->allGatherValue<int64_t>(value,maker));
      };

      /*!
      * wait for all transmissions/receptions to complete
      */
      void waitall()override{transceiver->waitall();};

      /*!
      * wait for all processes to arrive
      */
      void barrier()override{MPI_Barrier(netMPI_comm);};

    private:
      void ifNotInitializedError(){if(!isInitialized()){throw std::runtime_error("MPI Not initialized");};};
      int isInitialized(){int flag;MPI_Initialized(&flag);return flag;};
      void setupMPI(MPI_Comm comm){
        MPI_Comm_dup(comm,&netMPI_comm);
        MPI_Comm_rank(netMPI_comm,&rank);
        MPI_Comm_size(netMPI_comm,&comm_size);
        transceiver=std::make_shared<NetMPI::MPITransceiver>(NetMPI::MPITransceiver(comm));
        exchanger=std::make_shared<NetMPI::MPIHaloExchanger>(NetMPI::MPIHaloExchanger(comm));
      //  setupCounters();
      };
      //void setupCounters(){
      //  send_comm_counter.resize(0);send_comm_counter.resize(comm_size,0);
      //  recv_comm_counter.resize(0);recv_comm_counter.resize(comm_size,0);
      //};

      MPI_Comm netMPI_comm;
      int rank,comm_size;

      //std::vector<int> recv_comm_counter;
      //std::vector<int> send_comm_counter;

      std::shared_ptr<NetMPI::MPITransceiver> transceiver=nullptr;
      std::shared_ptr<NetMPI::MPIHaloExchanger> exchanger=nullptr;
  };
};
