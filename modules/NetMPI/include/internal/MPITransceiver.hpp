/*! \headerfile MPITransceiver.hpp "modules/NetMPI/include/internal/MPITransceiver.hpp"
* "MPITransceiver.hpp" implements basic point-to-point and collective MPI communication
* operations with type safety and asynchronous operation handling.
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

#include <map>
#include <vector>
#include <mpi.h>

#include "Array.hpp"

#include "MPIType.hpp"

namespace NetMPI {
  /*! \class MPITransceiver MPITransceiver.hpp "modules/NetMPI/include/internal/MPITransceiver.hpp"
  *
  * Implements low-level MPI communication operations with type safety and request management.
  * 
  * MPITransceiver provides templated methods for sending and receiving data between
  * processes using MPI, handling both point-to-point and collective operations.
  * It automatically manages communication tags and request objects for asynchronous
  * operations, and uses the MPIType template to ensure type-safe communications.
  * 
  * This class is primarily used by the MPICommunicator to implement its communication
  * operations while isolating the MPI-specific details.
  *
  */
  class MPITransceiver {
    public:
      /*!
      * Constructs a transceiver with a specified MPI communicator.
      * 
      * \param comm MPI_Comm The MPI communicator to use (default: MPI_COMM_WORLD)
      */
      MPITransceiver(MPI_Comm comm=MPI_COMM_WORLD){
        setupMPI(comm);
      };

      /*!
      * Asynchronously sends an array to a specific process.
      * 
      * This is a two-phase operation. First, the array length is sent,
      * then the array data itself is sent in a non-blocking fashion.
      * 
      * \tparam T The data type of array elements
      * \param array T* Pointer to the array data
      * \param length int Number of elements in the array
      * \param receiver int Rank of the receiving process
      */
      template<typename T>
      void isendArray(T* array,int length,int receiver){
        waitForSendTag(isendValue(length,receiver));
        int tag=createSendTag(receiver);
        MPI_Isend(array,length,NetMPI::MPIType<T>(),receiver,tag,
                  xceiver_comm,&send_request_queue.at(tag));
      };

      /*!
      * Asynchronously receives an array from a specific process.
      * 
      * This is a two-phase operation. First, the array length is received,
      * then a properly sized array is created and populated with the data.
      * 
      * \tparam T The data type of array elements
      * \param maker std::shared_ptr<PanNDE::ArrayFactory<T>> Factory for creating the array
      * \param sender int Rank of the sending process
      * \return std::shared_ptr<PanNDE::Array<T>> The received array
      */
      template<typename T>
      std::shared_ptr<PanNDE::Array<T>> irecvArray(
                      std::shared_ptr<PanNDE::ArrayFactory<T>> maker,int sender){
        int length;
        waitForRecvTag(irecvValue(&length,sender));
        int tag=createRecvTag(sender);
        auto array=maker->makeManagedArray();array->resize(length);
        MPI_Irecv(array->data(),length,NetMPI::MPIType<T>(),sender,tag,
                  xceiver_comm,&recv_request_queue.at(tag));
        return std::move(array);
      };

      /*!
      * Broadcasts an array from one process to all others.
      * 
      * This is a collective operation that ensures all processes have
      * the same array data. The array size is broadcast first, then the data.
      * 
      * \tparam T The data type of array elements
      * \param array std::shared_ptr<PanNDE::Array<T>>& Reference to the array
      * \param sender int Rank of the broadcasting process (default: 0)
      */
      template<typename T>
      void bcastArray(std::shared_ptr<PanNDE::Array<T>>& array,int sender){
        int N=array->size();
        MPI_Bcast(&N,1,NetMPI::MPIType<T>(),sender,xceiver_comm);
        if(rank!=sender){array->resize(N);};
        MPI_Bcast(array->data(),N,NetMPI::MPIType<T>(),sender,xceiver_comm);
      };

      /*!
      * Asynchronously sends a single value to a specific process.
      * 
      * \tparam T The data type of the value
      * \param value T The value to send
      * \param receiver int Rank of the receiving process
      * \return int The communication tag associated with this operation
      */
      template<typename T>
      int isendValue(T value,int receiver){
        int tag=createSendTag(receiver);
        MPI_Isend(&value,1,NetMPI::MPIType<T>(),receiver,tag,
                  xceiver_comm,&send_request_queue.at(tag));
        return tag;
      };

      /*!
      * Asynchronously receives a single value from a specific process.
      * 
      * \tparam T The data type of the value
      * \param value T* Pointer to where the received value should be stored
      * \param sender int Rank of the sending process
      * \return int The communication tag associated with this operation
      */
      template<typename T>
      int irecvValue(T* value,int sender){
        int tag=createRecvTag(sender);
        MPI_Irecv(value,1,NetMPI::MPIType<T>(),sender,tag,
                  xceiver_comm,&recv_request_queue.at(tag));
        return tag;
      };

      /*!
      * Broadcasts a single value from one process to all others.
      * 
      * This is a collective operation that ensures all processes have
      * the same value.
      * 
      * \tparam T The data type of the value
      * \param value T* Pointer to the value (source on sender, destination on receivers)
      * \param sender int Rank of the broadcasting process (default: 0)
      */
      template<typename T>
      void broadcastValue(T* value,int sender=0){
        MPI_Bcast(value,1,NetMPI::MPIType<T>(),sender,xceiver_comm);
      };

      /*!
      * Gathers a value from all processes into an array.
      * 
      * Each process contributes a value, and all processes receive
      * the complete set of values.
      * 
      * \tparam T The data type of the value
      * \param value T The local value to contribute
      * \param maker std::shared_ptr<PanNDE::ArrayFactory<T>> Factory for creating the result array
      * \return std::shared_ptr<PanNDE::Array<T>> Array containing gathered values
      */
      template<typename T>
      std::shared_ptr<PanNDE::Array<T>> allGatherValue(T value,
                      std::shared_ptr<PanNDE::ArrayFactory<T>> maker){
        auto array=maker->makeManagedArray();array->resize(comm_size);
        MPI_Allgather(&value,1,NetMPI::MPIType<T>(),
                      array->data(),1,NetMPI::MPIType<T>(),xceiver_comm);
        return std::move(array);
      };

      /*!
      * Waits for all pending send and receive operations to complete.
      * 
      * This method ensures that all asynchronous operations initiated by
      * this transceiver have completed before returning, and clears all
      * request queues.
      */
      void waitall(){
        for(auto it=send_request_queue.begin();it!=send_request_queue.end();it++){
          MPI_Wait(&(it->second),MPI_STATUS_IGNORE);
        };
        send_request_queue.clear();
        for(auto it=recv_request_queue.begin();it!=recv_request_queue.end();it++){
          MPI_Wait(&(it->second),MPI_STATUS_IGNORE);
        };
        recv_request_queue.clear();
      };

    private:
      /*!
      * Sets up the MPI environment and request tracking structures.
      * \param comm MPI_Comm The MPI communicator to use
      */
      void setupMPI(MPI_Comm comm){
        MPI_Comm_dup(comm,&xceiver_comm);
        MPI_Comm_rank(xceiver_comm,&rank);
        MPI_Comm_size(xceiver_comm,&comm_size);
        setupCounters();
      };

      /*!
      * Initializes communication counters for all processes.
      */
      void setupCounters(){
        send_comm_counter.resize(0);send_comm_counter.resize(comm_size,0);
        recv_comm_counter.resize(0);recv_comm_counter.resize(comm_size,0);
      };

      /*!
      * Creates a unique tag for a send operation and initializes its request object.
      * \param remote_rank int The rank of the target process
      * \return int The created tag value
      */
      int createSendTag(int remote_rank){
        int tag=send_comm_counter.at(remote_rank)*comm_size+remote_rank;
        send_comm_counter.at(remote_rank)++;
        send_request_queue.emplace(tag,MPI_Request());
        return tag;
      };

      /*!
      * Creates a unique tag for a receive operation and initializes its request object.
      * \param remote_rank int The rank of the source process
      * \return int The created tag value
      */
      int createRecvTag(int remote_rank){
        int tag=recv_comm_counter.at(remote_rank)*comm_size+rank;
        recv_comm_counter.at(remote_rank)++;
        recv_request_queue.emplace(tag,MPI_Request());
        return tag;
      };

      /*!
      * Waits for a specific send operation to complete.
      * \param tag int The tag identifying the operation
      */
      void waitForSendTag(int tag){
        MPI_Wait(&send_request_queue.at(tag),MPI_STATUS_IGNORE);
        send_request_queue.erase(tag);
      };

      /*!
      * Waits for a specific receive operation to complete.
      * \param tag int The tag identifying the operation
      */
      void waitForRecvTag(int tag){
        MPI_Wait(&recv_request_queue.at(tag),MPI_STATUS_IGNORE);
        recv_request_queue.erase(tag);
      };

      //! The MPI communicator used for all operations
      MPI_Comm xceiver_comm;
      //! Process rank within the communicator
      int rank;
      //! Total number of processes in the communicator
      int comm_size;
      //! Counters for generating unique receive tags for each remote process
      std::vector<int> recv_comm_counter;
      //! Counters for generating unique send tags for each remote process
      std::vector<int> send_comm_counter;
      //! Map of active receive requests indexed by tag
      std::map<int,MPI_Request> recv_request_queue;
      //! Map of active send requests indexed by tag
      std::map<int,MPI_Request> send_request_queue;
  };
};