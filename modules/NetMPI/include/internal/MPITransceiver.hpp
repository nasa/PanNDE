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

#include <map>
#include <vector>
#include <mpi.h>

#include "Array.hpp"

#include "MPIType.hpp"

namespace NetMPI {
  class MPITransceiver {
    public:
      MPITransceiver(MPI_Comm comm=MPI_COMM_WORLD){
        setupMPI(comm);
      };

      template<typename T>
      void isendArray(T* array,int length,int receiver){
        waitForSendTag(isendValue(length,receiver));
        int tag=createSendTag(receiver);
        MPI_Isend(array,length,NetMPI::MPIType<T>(),receiver,tag,
                  xceiver_comm,&send_request_queue.at(tag));
      };
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
      template<typename T>
      void bcastArray(std::shared_ptr<PanNDE::Array<T>>& array,int sender){
        int N=array->size();
        MPI_Bcast(&N,1,NetMPI::MPIType<T>(),sender,xceiver_comm);
        if(rank!=sender){array->resize(N);};
        MPI_Bcast(array->data(),N,NetMPI::MPIType<T>(),sender,xceiver_comm);
      };
      template<typename T>
      int isendValue(T value,int receiver){
        int tag=createSendTag(receiver);
        MPI_Isend(&value,1,NetMPI::MPIType<T>(),receiver,tag,
                  xceiver_comm,&send_request_queue.at(tag));
        return tag;
      };
      template<typename T>
      int irecvValue(T* value,int sender){
        int tag=createRecvTag(sender);
        MPI_Irecv(value,1,NetMPI::MPIType<T>(),sender,tag,
                  xceiver_comm,&recv_request_queue.at(tag));
        return tag;
      };
      template<typename T>
      void broadcastValue(T* value,int sender=0){
        MPI_Bcast(value,1,NetMPI::MPIType<T>(),sender,xceiver_comm);
      };
      template<typename T>
      std::shared_ptr<PanNDE::Array<T>> allGatherValue(T value,
                      std::shared_ptr<PanNDE::ArrayFactory<T>> maker){
        auto array=maker->makeManagedArray();array->resize(comm_size);
        MPI_Allgather(&value,1,NetMPI::MPIType<T>(),
                      array->data(),1,NetMPI::MPIType<T>(),xceiver_comm);
        return std::move(array);
      };

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
      void setupMPI(MPI_Comm comm){
        MPI_Comm_dup(comm,&xceiver_comm);
        MPI_Comm_rank(xceiver_comm,&rank);
        MPI_Comm_size(xceiver_comm,&comm_size);
        setupCounters();
      };
      void setupCounters(){
        send_comm_counter.resize(0);send_comm_counter.resize(comm_size,0);
        recv_comm_counter.resize(0);recv_comm_counter.resize(comm_size,0);
      };
      int createSendTag(int remote_rank){
        int tag=send_comm_counter.at(remote_rank)*comm_size+remote_rank;
        send_comm_counter.at(remote_rank)++;
        send_request_queue.emplace(tag,MPI_Request());
        return tag;
      };
      int createRecvTag(int remote_rank){
        int tag=recv_comm_counter.at(remote_rank)*comm_size+rank;
        recv_comm_counter.at(remote_rank)++;
        recv_request_queue.emplace(tag,MPI_Request());
        return tag;
      };
      void waitForSendTag(int tag){
        MPI_Wait(&send_request_queue.at(tag),MPI_STATUS_IGNORE);
        send_request_queue.erase(tag);
      };
      void waitForRecvTag(int tag){
        MPI_Wait(&recv_request_queue.at(tag),MPI_STATUS_IGNORE);
        recv_request_queue.erase(tag);
      };
      MPI_Comm xceiver_comm;
      int rank,comm_size;
      std::vector<int> recv_comm_counter;
      std::vector<int> send_comm_counter;
      std::map<int,MPI_Request> recv_request_queue;
      std::map<int,MPI_Request> send_request_queue;
  };
};