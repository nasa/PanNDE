/*! \headerfile Communicator.hpp "modules/PanNDE/include/Communicator.hpp"
* "Communicator.hpp" defines the abstract interface for parallel communication in PanNDE.
* It encapsulates all inter-process data exchange operations, isolating the rest of the
* framework from the specific communication library (e.g., MPI) implementation details.
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

#include "Mesh.hpp"
#include "Field.hpp"
#include "Array.hpp"

namespace PanNDE {
  /*! \class Communicator Communicator.hpp "modules/PanNDE/include/Communicator.hpp"
  *
  * Defines the interface for parallel communication operations in distributed simulations.
  * 
  * The Communicator provides an abstraction layer over specific parallel communication
  * libraries (such as MPI), allowing the rest of the PanNDE framework to remain
  * independent of the underlying communication implementation. This supports:
  *
  * 1. Process management and identification
  * 2. Point-to-point communication (send/receive)
  * 3. Collective operations (broadcast, gather)
  * 4. Mesh-based halo exchange patterns for ghost region updates
  * 5. Synchronization primitives (barriers)
  *
  * Implementations of this interface should handle all necessary setup, data
  * transfer, and cleanup operations for the target parallel environment.
  *
  */
  class Communicator {
    public:
      //----------------------------------------------------------------------
      // Initialization and Process Information
      //----------------------------------------------------------------------
      
      /*! 
      * Initializes the communication system.
      * 
      * Many communication libraries (e.g., MPI) require access to command-line arguments
      * for initialization. This method handles that initialization process.
      * 
      * \param argc Reference to the count of command-line arguments
      * \param argv Reference to the command-line argument array
      */
      virtual void Init(int& argc,char**& argv) =0;
      
      /*! 
      * Initializes the communication system with mesh information.
      * 
      * This variant accepts a mesh to determine communication patterns during
      * initialization, optimizing the setup for domain decomposition.
      * 
      * \param argc Reference to the count of command-line arguments
      * \param argv Reference to the command-line argument array
      * \param mesh The local partition mesh used to determine communication patterns
      */
      virtual void Init(int& argc,char**& argv,std::shared_ptr<PanNDE::Mesh> mesh) =0;
      
      /*!
      * Gets the rank (ID) of this process within its local communicator.
      * 
      * This ID uniquely identifies the process within its local communication group.
      * 
      * \return int Process rank in the local communicator (0-based)
      */
      virtual int getProcessId() =0;
      
      /*!
      * Gets the global rank (ID) of this process.
      * 
      * For hierarchical communicators, this returns the process ID in the
      * global/world communicator context.
      * 
      * \return int Global process rank (0-based)
      */
      virtual int getGlobalId() =0;
      
      /*!
      * Gets the total number of processes in the local communicator.
      * 
      * \return int Number of processes in the local communicator
      */
      virtual int getNumberOfProcesses() =0;
      
      /*!
      * Gets the total number of processes across all communicators.
      * 
      * For hierarchical communicators, this returns the total process count
      * in the global/world context.
      * 
      * \return int Total number of processes globally
      */
      virtual int getNumberOfAllProcesses() =0;
      
      /*!
      * Finalizes the communication system and releases resources.
      * 
      * This method should be called before program termination to ensure
      * proper cleanup of communication resources.
      */
      virtual void Finalize() =0;

      //----------------------------------------------------------------------
      // Mesh-Based Communication Pattern Setup
      //----------------------------------------------------------------------
      
      /*!
      * Analyzes a mesh to determine the communication pattern for halo exchanges.
      * 
      * This method examines mesh connectivity across partition boundaries to
      * establish which processes need to exchange data and which mesh elements
      * are involved in those exchanges.
      * 
      * Note: This is retained for backward compatibility. Using setupDataLinks
      * with a FieldBundle is recommended for new code.
      * 
      * \param mesh The partitioned mesh to analyze
      */
      virtual void determineExchangePattern(std::shared_ptr<PanNDE::Mesh> mesh) =0;
      
      /*!
      * Registers a field for halo exchange using a previously determined pattern.
      * 
      * This method prepares a specific field for halo exchange operations based
      * on the pattern determined by determineExchangePattern().
      * 
      * Note: This is retained for backward compatibility. Using setupDataLinks
      * with a FieldBundle is recommended for new code.
      * 
      * \param keyname Name identifier for the field
      * \param field The field to register for halo exchange
      */
      virtual void setupDataLinks(std::string keyname,std::shared_ptr<PanNDE::Field> field) =0;
      
      /*!
      * Sets up halo exchange for all fields in a field bundle.
      * 
      * This convenience method determines the exchange pattern from the bundle's
      * mesh and registers all contained fields for halo exchange in a single operation.
      * This is the preferred approach for configuring parallel data exchange.
      * 
      * \param fields The field bundle containing all fields to register
      */
      virtual void setupDataLinks(std::shared_ptr<PanNDE::FieldBundle> fields) =0;
      
      //----------------------------------------------------------------------
      // Asynchronous Halo Exchange Operations
      //----------------------------------------------------------------------
      
      /*!
      * Initiates an asynchronous halo exchange operation for a specific field.
      * 
      * This starts the communication to update ghost regions for the named field
      * but does not wait for completion.
      * 
      * \param keyname Name identifier of the field to update
      */
      virtual void startHaloExchange(std::string keyname) =0;
      
      /*!
      * Waits for completion of a previously initiated halo exchange.
      * 
      * This blocks until the asynchronous halo exchange for the named field
      * has completed, ensuring that ghost regions are fully updated.
      * 
      * \param keyname Name identifier of the field being updated
      */
      virtual void waitUntilDone(std::string keyname) =0;

      //----------------------------------------------------------------------
      // Double Array Communication
      //----------------------------------------------------------------------
      
      /*!
      * Sends a double array to another process.
      * 
      * \param array The array to send
      * \param receiver Rank of the receiving process
      */
      virtual void sendArray(std::shared_ptr<PanNDE::Array<double>> array,int receiver) =0;
      
      /*!
      * Sends a raw double array to another process.
      * 
      * \param array Pointer to the array data
      * \param N Number of elements in the array
      * \param receiver Rank of the receiving process
      */
      virtual void sendArray(double* array,int N,int receiver) =0;
      
      /*!
      * Receives a double array from another process.
      * 
      * \param maker Array factory to create the received array
      * \param sender Rank of the sending process
      * \return std::shared_ptr<PanNDE::Array<double>> The received array
      */
      virtual std::shared_ptr<PanNDE::Array<double>> recvArray(
                              std::shared_ptr<PanNDE::ArrayFactory<double>> maker,int sender) =0;
      
      /*!
      * Broadcasts a double array from one process to all others.
      * 
      * For the sender, the array contains the data to broadcast.
      * For receivers, the array will be resized and populated with received data.
      * 
      * \param array Reference to the array (source on sender, destination on receivers)
      * \param sender Rank of the broadcasting process (default: 0)
      */
      virtual void broadcastArray(std::shared_ptr<PanNDE::Array<double>>& array,int sender=0) =0;

      //----------------------------------------------------------------------
      // Int32 Array Communication
      //----------------------------------------------------------------------
      
      /*!
      * Sends an int32 array to another process.
      * 
      * \param array The array to send
      * \param receiver Rank of the receiving process
      */
      virtual void sendArray(std::shared_ptr<PanNDE::Array<int32_t>> array,int receiver) =0;
      
      /*!
      * Sends a raw int32 array to another process.
      * 
      * \param array Pointer to the array data
      * \param N Number of elements in the array
      * \param receiver Rank of the receiving process
      */
      virtual void sendArray(int32_t* array,int N,int receiver) =0;
      
      /*!
      * Receives an int32 array from another process.
      * 
      * \param maker Array factory to create the received array
      * \param sender Rank of the sending process
      * \return std::shared_ptr<PanNDE::Array<int32_t>> The received array
      */
      virtual std::shared_ptr<PanNDE::Array<int32_t>> recvArray(
                              std::shared_ptr<PanNDE::ArrayFactory<int32_t>> maker,int sender) =0;
      
      /*!
      * Broadcasts an int32 array from one process to all others.
      * 
      * For the sender, the array contains the data to broadcast.
      * For receivers, the array will be resized and populated with received data.
      * 
      * \param array Reference to the array (source on sender, destination on receivers)
      * \param sender Rank of the broadcasting process (default: 0)
      */
      virtual void broadcastArray(std::shared_ptr<PanNDE::Array<int32_t>>& array,int sender=0) =0;

      //----------------------------------------------------------------------
      // Int64 Array Communication
      //----------------------------------------------------------------------
      
      /*!
      * Sends an int64 array to another process.
      * 
      * \param array The array to send
      * \param receiver Rank of the receiving process
      */
      virtual void sendArray(std::shared_ptr<PanNDE::Array<int64_t>> array,int receiver) =0;
      
      /*!
      * Sends a raw int64 array to another process.
      * 
      * \param array Pointer to the array data
      * \param N Number of elements in the array
      * \param receiver Rank of the receiving process
      */
      virtual void sendArray(int64_t* array,int N,int receiver) =0;
      
      /*!
      * Receives an int64 array from another process.
      * 
      * \param maker Array factory to create the received array
      * \param sender Rank of the sending process
      * \return std::shared_ptr<PanNDE::Array<int64_t>> The received array
      */
      virtual std::shared_ptr<PanNDE::Array<int64_t>> recvArray(
                              std::shared_ptr<PanNDE::ArrayFactory<int64_t>> maker,int sender) =0;
      
      /*!
      * Broadcasts an int64 array from one process to all others.
      * 
      * For the sender, the array contains the data to broadcast.
      * For receivers, the array will be resized and populated with received data.
      * 
      * \param array Reference to the array (source on sender, destination on receivers)
      * \param sender Rank of the broadcasting process (default: 0)
      */
      virtual void broadcastArray(std::shared_ptr<PanNDE::Array<int64_t>>& array,int sender=0) =0;

      //----------------------------------------------------------------------
      // Double Scalar Communication
      //----------------------------------------------------------------------
      
      /*!
      * Sends a double value to another process.
      * 
      * \param value The value to send
      * \param receiver Rank of the receiving process
      */
      virtual void sendValue(double value,int receiver) =0;
      
      /*!
      * Receives a double value from another process.
      * 
      * \param value Pointer to the location where the received value will be stored
      * \param sender Rank of the sending process
      */
      virtual void recvValue(double* value,int sender) =0;
      
      /*!
      * Broadcasts a double value from one process to all others.
      * 
      * \param value Pointer to the value (source on sender, destination on receivers)
      * \param sender Rank of the broadcasting process (default: 0)
      */
      virtual void broadcastValue(double* value,int sender=0) =0;
      
      /*!
      * Gathers a double value from all processes into an array.
      * 
      * Each process contributes a value, and all processes receive an array
      * containing all contributed values ordered by process rank.
      * 
      * \param value The local value to contribute
      * \param maker Array factory to create the result array
      * \return std::shared_ptr<PanNDE::Array<double>> Array of gathered values
      */
      virtual std::shared_ptr<PanNDE::Array<double>> allGatherValue(double value,
                              std::shared_ptr<PanNDE::ArrayFactory<double>> maker) =0;

      //----------------------------------------------------------------------
      // Int32 Scalar Communication
      //----------------------------------------------------------------------
      
      /*!
      * Sends an int32 value to another process.
      * 
      * \param value The value to send
      * \param receiver Rank of the receiving process
      */
      virtual void sendValue(int32_t value,int receiver) =0;
      
      /*!
      * Receives an int32 value from another process.
      * 
      * \param value Pointer to the location where the received value will be stored
      * \param sender Rank of the sending process
      */
      virtual void recvValue(int32_t* value,int sender) =0;
      
      /*!
      * Broadcasts an int32 value from one process to all others.
      * 
      * \param value Pointer to the value (source on sender, destination on receivers)
      * \param sender Rank of the broadcasting process (default: 0)
      */
      virtual void broadcastValue(int32_t* value,int sender=0) =0;
      
      /*!
      * Gathers an int32 value from all processes into an array.
      * 
      * Each process contributes a value, and all processes receive an array
      * containing all contributed values ordered by process rank.
      * 
      * \param value The local value to contribute
      * \param maker Array factory to create the result array
      * \return std::shared_ptr<PanNDE::Array<int32_t>> Array of gathered values
      */
      virtual std::shared_ptr<PanNDE::Array<int32_t>> allGatherValue(int32_t value,
                              std::shared_ptr<PanNDE::ArrayFactory<int32_t>> maker) =0;

      //----------------------------------------------------------------------
      // Int64 Scalar Communication
      //----------------------------------------------------------------------
      
      /*!
      * Sends an int64 value to another process.
      * 
      * \param value The value to send
      * \param receiver Rank of the receiving process
      */
      virtual void sendValue(int64_t value,int receiver) =0;
      
      /*!
      * Receives an int64 value from another process.
      * 
      * \param value Pointer to the location where the received value will be stored
      * \param sender Rank of the sending process
      */
      virtual void recvValue(int64_t* value,int sender) =0;
      
      /*!
      * Broadcasts an int64 value from one process to all others.
      * 
      * \param value Pointer to the value (source on sender, destination on receivers)
      * \param sender Rank of the broadcasting process (default: 0)
      */
      virtual void broadcastValue(int64_t* value,int sender=0) =0;
      
      /*!
      * Gathers an int64 value from all processes into an array.
      * 
      * Each process contributes a value, and all processes receive an array
      * containing all contributed values ordered by process rank.
      * 
      * \param value The local value to contribute
      * \param maker Array factory to create the result array
      * \return std::shared_ptr<PanNDE::Array<int64_t>> Array of gathered values
      */
      virtual std::shared_ptr<PanNDE::Array<int64_t>> allGatherValue(int64_t value,
                              std::shared_ptr<PanNDE::ArrayFactory<int64_t>> maker) =0;

      //----------------------------------------------------------------------
      // Synchronization Operations
      //----------------------------------------------------------------------
      
      /*!
      * Waits for all pending communication operations to complete.
      * 
      * This ensures that all previously initiated send and receive operations
      * have finished before proceeding.
      */
      virtual void waitall() =0;
      
      /*!
      * Creates a synchronization point for all processes.
      * 
      * This blocks until all processes have reached this barrier, ensuring
      * global synchronization across the communicator.
      */
      virtual void barrier() =0;
  };
};