/*! \headerfile MPICommunicator.hpp "modules/NetMPI/include/MPICommunicator.hpp"
* "MPICommunicator.hpp" implements the PanNDE::Communicator interface using the Message Passing
* Interface (MPI) standard for parallel computing. It encapsulates all MPI-specific operations,
* isolating the rest of the framework from direct MPI dependencies.
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
  /*! \class MPICommunicator MPICommunicator.hpp "modules/NetMPI/include/MPICommunicator.hpp"
  *
  * Implements the PanNDE::Communicator interface using MPI for parallel communication.
  * 
  * This class provides a complete MPI-based implementation of the PanNDE communication
  * abstraction layer, handling both point-to-point and collective operations as well as
  * specialized mesh-based halo exchanges. It delegates specialized operations to two
  * helper classes:
  * 
  * - MPITransceiver: Handles basic data transfer operations (send/receive, broadcast, etc.)
  * - MPIHaloExchanger: Manages ghost region updates based on mesh partitioning
  * 
  * The class carefully manages MPI initialization state and provides a clean interface
  * that shields the rest of the framework from direct MPI dependencies and complexities.
  * 
  */
  class MPICommunicator : public PanNDE::Communicator {
    public:
      /*!
      * Constructs a communicator with an optional MPI communicator.
      * 
      * Will only perform initialization if MPI is already initialized.
      * 
      * \param comm MPI_Comm The MPI communicator to use (default: MPI_COMM_WORLD)
      */
      MPICommunicator(MPI_Comm comm=MPI_COMM_WORLD){
        if(isInitialized()){setupMPI(comm);};
      };

      /*!
      * Constructs a communicator with a mesh and optional MPI communicator.
      * 
      * This variant immediately configures halo exchange patterns based on the
      * provided mesh, if MPI is already initialized.
      * 
      * \param mesh std::shared_ptr<PanNDE::Mesh> Mesh for determining communication patterns
      * \param comm MPI_Comm The MPI communicator to use (default: MPI_COMM_WORLD)
      */
      MPICommunicator(std::shared_ptr<PanNDE::Mesh> mesh,MPI_Comm comm=MPI_COMM_WORLD){
        if(isInitialized()){
          setupMPI(comm);
          determineExchangePattern(mesh);
        };
      };

      /*!
      * Destructor.
      */
      ~MPICommunicator(){};

      /*!
      * Creates a shared pointer to an MPICommunicator.
      * 
      * Initializes MPI if not already initialized, then creates and configures a communicator.
      * 
      * \param comm MPI_Comm The MPI communicator to use (default: MPI_COMM_WORLD)
      * \return std::shared_ptr<NetMPI::MPICommunicator> The created communicator
      */
      static inline
      std::shared_ptr<NetMPI::MPICommunicator> makeShared(MPI_Comm comm=MPI_COMM_WORLD){
        if(!isInitialized()){MPI_Init(NULL,NULL);};
        auto communicator=std::make_shared<NetMPI::MPICommunicator>(NetMPI::MPICommunicator(comm));
        return std::move(communicator);
      };

      /*!
      * Creates a shared pointer to an MPICommunicator using command-line arguments.
      * 
      * Initializes MPI with the provided arguments if not already initialized,
      * then creates and configures a communicator.
      * 
      * \param argc int& Reference to argument count
      * \param argv char**& Reference to argument array
      * \param comm MPI_Comm The MPI communicator to use (default: MPI_COMM_WORLD)
      * \return std::shared_ptr<NetMPI::MPICommunicator> The created communicator
      */
      static inline
      std::shared_ptr<NetMPI::MPICommunicator> makeShared(int& argc, char**& argv, 
                                                          MPI_Comm comm=MPI_COMM_WORLD){
        if(!isInitialized()){MPI_Init(&argc,&argv);};
        auto communicator=std::make_shared<NetMPI::MPICommunicator>(NetMPI::MPICommunicator(comm));
        return std::move(communicator);
      };

      /*!
      * Initializes the MPI communication system with command-line arguments.
      * 
      * If MPI is not already initialized, this method initializes it with the
      * provided arguments and configures the communicator.
      * 
      * \param argc int& Reference to argument count
      * \param argv char**& Reference to argument array
      */
      void Init(int& argc,char**& argv) override {
        if(!isInitialized()){
          MPI_Init(&argc,&argv);
          setupMPI(MPI_COMM_WORLD);
        };
      };

      /*!
      * Initializes the MPI communication system with arguments and mesh information.
      * 
      * This variant initializes MPI with the provided arguments if necessary,
      * then configures the communicator and sets up mesh-based exchange patterns.
      * 
      * \param argc int& Reference to argument count
      * \param argv char**& Reference to argument array
      * \param mesh std::shared_ptr<PanNDE::Mesh> Mesh for determining communication patterns
      */
      void Init(int& argc,char**& argv,std::shared_ptr<PanNDE::Mesh> mesh) override {
        Init(argc,argv);
        determineExchangePattern(mesh);
      };

      /*!
      * Gets the process rank in the local communicator.
      * 
      * \return int The process rank (0-based)
      * \throws std::runtime_error If MPI is not initialized
      */
      int getProcessId() override {ifNotInitializedError();return rank;};

      /*!
      * Gets the process rank in the global (world) communicator.
      * 
      * \return int The global process rank (0-based)
      * \throws std::runtime_error If MPI is not initialized
      */
      int getGlobalId() override {
        ifNotInitializedError();
        int grank;
        MPI_Comm_rank(MPI_COMM_WORLD,&grank);
        return grank;
      };

      /*!
      * Gets the number of processes in the local communicator.
      * 
      * \return int The number of processes
      * \throws std::runtime_error If MPI is not initialized
      */
      int getNumberOfProcesses() override {ifNotInitializedError();return comm_size;};

      /*!
      * Gets the number of processes in the global communicator.
      * 
      * \return int The total number of processes
      * \throws std::runtime_error If MPI is not initialized
      */
      int getNumberOfAllProcesses() override {
        ifNotInitializedError();
        int gcomm_size;
        MPI_Comm_size(netMPI_comm,&gcomm_size);
        return gcomm_size;
      };

      /*!
      * Finalizes the MPI communication system.
      * 
      * This should be called before program termination to ensure
      * proper cleanup of MPI resources, but only if MPI hasn't already
      * been finalized.
      */
      void Finalize() override {
        int flag;MPI_Finalized(&flag);
        if(!flag){MPI_Finalize();};
      };

      /*!
      * Determines the communication pattern for halo exchanges based on mesh topology.
      * 
      * Analyzes the given mesh to establish which processes need to exchange data
      * and which mesh elements are involved in those exchanges.
      * 
      * \param mesh std::shared_ptr<PanNDE::Mesh> The mesh to analyze
      */
      void determineExchangePattern(std::shared_ptr<PanNDE::Mesh> mesh) override {
        exchanger->determineExchangePattern(mesh);
      };

      /*!
      * Registers a field for halo exchange using the previously determined pattern.
      * 
      * \param keyname std::string Name identifier for the field
      * \param field std::shared_ptr<PanNDE::Field> The field to register
      */
      void setupDataLinks(std::string keyname,std::shared_ptr<PanNDE::Field> field) override {
        exchanger->setupDataLinks(keyname,field);
      };

      /*!
      * Sets up halo exchange for all fields in a field bundle.
      * 
      * This is the preferred method for configuring halo exchanges as it
      * handles all fields in the bundle at once.
      * 
      * \param fields std::shared_ptr<PanNDE::FieldBundle> Bundle containing fields to register
      */
      void setupDataLinks(std::shared_ptr<PanNDE::FieldBundle> fields) override {
        exchanger->setupDataLinks(fields);
      };
      
      /*!
      * Initiates an asynchronous halo exchange operation for a field.
      * 
      * \param keyname std::string Name of the field to update
      */
      void startHaloExchange(std::string keyname) override {
        exchanger->startHaloExchange(keyname);
      };

      /*!
      * Waits for completion of a previously initiated halo exchange.
      * 
      * \param keyname std::string Name of the field being updated
      */
      void waitUntilDone(std::string keyname) override {
        exchanger->waitUntilDone(keyname);
      };

      //----------------------------------------------------------------------
      // Double Array Communication
      //----------------------------------------------------------------------

      /*!
      * Sends a double array to another process.
      * 
      * \param array std::shared_ptr<PanNDE::Array<double>> The array to send
      * \param receiver int Rank of the receiving process
      */
      void sendArray(std::shared_ptr<PanNDE::Array<double>> array,int receiver) override {
        sendArray(array->data(),array->size(),receiver);
      };

      /*!
      * Sends a raw double array to another process.
      * 
      * \param array double* Pointer to the array data
      * \param N int Number of elements in the array
      * \param receiver int Rank of the receiving process
      */
      void sendArray(double* array,int N,int receiver) override {
        transceiver->isendArray<double>(array,N,receiver);
      };

      /*!
      * Receives a double array from another process.
      * 
      * \param maker std::shared_ptr<PanNDE::ArrayFactory<double>> Factory to create the array
      * \param sender int Rank of the sending process
      * \return std::shared_ptr<PanNDE::Array<double>> The received array
      */
      std::shared_ptr<PanNDE::Array<double>> recvArray(
                      std::shared_ptr<PanNDE::ArrayFactory<double>> maker,int sender) override {
        return std::move(transceiver->irecvArray<double>(maker,sender));
      };

      /*!
      * Broadcasts a double array from one process to all others.
      * 
      * \param array std::shared_ptr<PanNDE::Array<double>>& Reference to the array
      * \param sender int Rank of the broadcasting process (default: 0)
      */
      void broadcastArray(std::shared_ptr<PanNDE::Array<double>>& array,int sender=0) override {
        transceiver->bcastArray<double>(array,sender);
      };

      //----------------------------------------------------------------------
      // Int32 Array Communication
      //----------------------------------------------------------------------

      /*!
      * Sends an int32 array to another process.
      * 
      * \param array std::shared_ptr<PanNDE::Array<int32_t>> The array to send
      * \param receiver int Rank of the receiving process
      */
      void sendArray(std::shared_ptr<PanNDE::Array<int32_t>> array,int receiver) override {
        sendArray(array->data(),array->size(),receiver);
      };

      /*!
      * Sends a raw int32 array to another process.
      * 
      * \param array int32_t* Pointer to the array data
      * \param N int Number of elements in the array
      * \param receiver int Rank of the receiving process
      */
      void sendArray(int32_t* array,int N,int receiver) override {
        transceiver->isendArray<int32_t>(array,N,receiver);
      };

      /*!
      * Receives an int32 array from another process.
      * 
      * \param maker std::shared_ptr<PanNDE::ArrayFactory<int32_t>> Factory to create the array
      * \param sender int Rank of the sending process
      * \return std::shared_ptr<PanNDE::Array<int32_t>> The received array
      */
      std::shared_ptr<PanNDE::Array<int32_t>> recvArray(
                      std::shared_ptr<PanNDE::ArrayFactory<int32_t>> maker,int sender) override {
        return std::move(transceiver->irecvArray<int32_t>(maker,sender));
      };

      /*!
      * Broadcasts an int32 array from one process to all others.
      * 
      * \param array std::shared_ptr<PanNDE::Array<int32_t>>& Reference to the array
      * \param sender int Rank of the broadcasting process (default: 0)
      */
      void broadcastArray(std::shared_ptr<PanNDE::Array<int32_t>>& array,int sender=0) override {
        transceiver->bcastArray<int32_t>(array,sender);
      };

      //----------------------------------------------------------------------
      // Int64 Array Communication
      //----------------------------------------------------------------------

      /*!
      * Sends an int64 array to another process.
      * 
      * \param array std::shared_ptr<PanNDE::Array<int64_t>> The array to send
      * \param receiver int Rank of the receiving process
      */
      void sendArray(std::shared_ptr<PanNDE::Array<int64_t>> array,int receiver) override {
        sendArray(array->data(),array->size(),receiver);
      };

      /*!
      * Sends a raw int64 array to another process.
      * 
      * \param array int64_t* Pointer to the array data
      * \param N int Number of elements in the array
      * \param receiver int Rank of the receiving process
      */
      void sendArray(int64_t* array,int N,int receiver) override {
        transceiver->isendArray<int64_t>(array,N,receiver);
      };

      /*!
      * Receives an int64 array from another process.
      * 
      * \param maker std::shared_ptr<PanNDE::ArrayFactory<int64_t>> Factory to create the array
      * \param sender int Rank of the sending process
      * \return std::shared_ptr<PanNDE::Array<int64_t>> The received array
      */
      std::shared_ptr<PanNDE::Array<int64_t>> recvArray(
                      std::shared_ptr<PanNDE::ArrayFactory<int64_t>> maker,int sender) override {
        return std::move(transceiver->irecvArray<int64_t>(maker,sender));
      };

      /*!
      * Broadcasts an int64 array from one process to all others.
      * 
      * \param array std::shared_ptr<PanNDE::Array<int64_t>>& Reference to the array
      * \param sender int Rank of the broadcasting process (default: 0)
      */
      void broadcastArray(std::shared_ptr<PanNDE::Array<int64_t>>& array,int sender=0) override {
        transceiver->bcastArray<int64_t>(array,sender);
      };

      //----------------------------------------------------------------------
      // Double Scalar Communication
      //----------------------------------------------------------------------

      /*!
      * Sends a double value to another process.
      * 
      * \param value double The value to send
      * \param receiver int Rank of the receiving process
      */
      void sendValue(double value,int receiver) override {
        transceiver->isendValue<double>(value,receiver);
      };

      /*!
      * Receives a double value from another process.
      * 
      * \param value double* Pointer to store the received value
      * \param sender int Rank of the sending process
      */
      void recvValue(double* value,int sender) override {
        transceiver->irecvValue<double>(value,sender);
      };

      /*!
      * Broadcasts a double value from one process to all others.
      * 
      * \param value double* Pointer to the value (source or destination)
      * \param sender int Rank of the broadcasting process (default: 0)
      */
      void broadcastValue(double* value,int sender=0) override {
        transceiver->broadcastValue<double>(value,sender);
      };

      /*!
      * Gathers a double value from all processes into an array.
      * 
      * \param value double The local value to contribute
      * \param maker std::shared_ptr<PanNDE::ArrayFactory<double>> Factory to create the result array
      * \return std::shared_ptr<PanNDE::Array<double>> Array of gathered values
      */
      std::shared_ptr<PanNDE::Array<double>> allGatherValue(double value,
                      std::shared_ptr<PanNDE::ArrayFactory<double>> maker) override {
        return std::move(transceiver->allGatherValue<double>(value,maker));
      };

      //----------------------------------------------------------------------
      // Int32 Scalar Communication
      //----------------------------------------------------------------------

      /*!
      * Sends an int32 value to another process.
      * 
      * \param value int32_t The value to send
      * \param receiver int Rank of the receiving process
      */
      void sendValue(int32_t value,int receiver) override {
        transceiver->isendValue<int32_t>(value,receiver);
      };

      /*!
      * Receives an int32 value from another process.
      * 
      * \param value int32_t* Pointer to store the received value
      * \param sender int Rank of the sending process
      */
      void recvValue(int32_t* value,int sender) override {
        transceiver->irecvValue<int32_t>(value,sender);
      };

      /*!
      * Broadcasts an int32 value from one process to all others.
      * 
      * \param value int32_t* Pointer to the value (source or destination)
      * \param sender int Rank of the broadcasting process (default: 0)
      */
      void broadcastValue(int32_t* value,int sender=0) override {
        transceiver->broadcastValue<int32_t>(value,sender);
      };

      /*!
      * Gathers an int32 value from all processes into an array.
      * 
      * \param value int32_t The local value to contribute
      * \param maker std::shared_ptr<PanNDE::ArrayFactory<int32_t>> Factory to create the result array
      * \return std::shared_ptr<PanNDE::Array<int32_t>> Array of gathered values
      */
      std::shared_ptr<PanNDE::Array<int32_t>> allGatherValue(int32_t value,
                      std::shared_ptr<PanNDE::ArrayFactory<int32_t>> maker) override {
        return std::move(transceiver->allGatherValue<int32_t>(value,maker));
      };

      //----------------------------------------------------------------------
      // Int64 Scalar Communication
      //----------------------------------------------------------------------

      /*!
      * Sends an int64 value to another process.
      * 
      * \param value int64_t The value to send
      * \param receiver int Rank of the receiving process
      */
      void sendValue(int64_t value,int receiver) override {
        transceiver->isendValue<int64_t>(value,receiver);
      };

      /*!
      * Receives an int64 value from another process.
      * 
      * \param value int64_t* Pointer to store the received value
      * \param sender int Rank of the sending process
      */
      void recvValue(int64_t* value,int sender) override {
        transceiver->irecvValue<int64_t>(value,sender);
      };

      /*!
      * Broadcasts an int64 value from one process to all others.
      * 
      * \param value int64_t* Pointer to the value (source or destination)
      * \param sender int Rank of the broadcasting process (default: 0)
      */
      void broadcastValue(int64_t* value,int sender=0) override {
        transceiver->broadcastValue<int64_t>(value,sender);
      };

      /*!
      * Gathers an int64 value from all processes into an array.
      * 
      * \param value int64_t The local value to contribute
      * \param maker std::shared_ptr<PanNDE::ArrayFactory<int64_t>> Factory to create the result array
      * \return std::shared_ptr<PanNDE::Array<int64_t>> Array of gathered values
      */
      std::shared_ptr<PanNDE::Array<int64_t>> allGatherValue(int64_t value,
                      std::shared_ptr<PanNDE::ArrayFactory<int64_t>> maker) override {
        return std::move(transceiver->allGatherValue<int64_t>(value,maker));
      };

      //----------------------------------------------------------------------
      // Synchronization Operations
      //----------------------------------------------------------------------

      /*!
      * Waits for all pending communication operations to complete.
      */
      void waitall() override {transceiver->waitall();};

      /*!
      * Creates a synchronization point for all processes.
      */
      void barrier() override {MPI_Barrier(netMPI_comm);};

    private:
      /*!
      * Checks if MPI is initialized and throws an error if not.
      * \throws std::runtime_error If MPI is not initialized
      */
      static
      void ifNotInitializedError(){if(!isInitialized()){throw std::runtime_error("MPI Not initialized");};};
      
      /*!
      * Checks if MPI is initialized.
      * \return int Non-zero if MPI is initialized, zero otherwise
      */
      static
      int isInitialized(){int flag;int err=MPI_Initialized(&flag);return flag;};
      
      /*!
      * Sets up the MPI environment and helper objects.
      * \param comm MPI_Comm The MPI communicator to use
      */
      void setupMPI(MPI_Comm comm){
        MPI_Comm_dup(comm,&netMPI_comm);
        MPI_Comm_rank(netMPI_comm,&rank);
        MPI_Comm_size(netMPI_comm,&comm_size);
        transceiver=std::make_shared<NetMPI::MPITransceiver>(NetMPI::MPITransceiver(comm));
        exchanger=std::make_shared<NetMPI::MPIHaloExchanger>(NetMPI::MPIHaloExchanger(comm));
      };

      //! MPI communicator used by this instance
      MPI_Comm netMPI_comm;
      
      //! Process rank in the local communicator
      int rank;
      
      //! Number of processes in the local communicator
      int comm_size;

      //! Helper object for point-to-point and collective communications
      std::shared_ptr<NetMPI::MPITransceiver> transceiver=nullptr;
      
      //! Helper object for mesh-based halo exchanges
      std::shared_ptr<NetMPI::MPIHaloExchanger> exchanger=nullptr;
  };
};