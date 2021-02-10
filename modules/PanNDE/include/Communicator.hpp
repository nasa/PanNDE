/*! \headerfile Communicator.hpp "modules/PanNDE/include/Communicator.hpp"
* "Communicator.hpp" contains the class definition encapsulating 
* the data communication among processes. This interface is intended to prevent leakage of 
* any chosen inter-process communications library (e.g. MPI) into other parts of the software stack
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
  * Defines the data communcations methods for inter-process data exchange.
  *
  */
  class Communicator {
    public:
      /*! 
      * Initialize a communicator object. As many initializers (e.g. MPI, VTK, etc.) 
      * require the arguments to `int main()', they are passed through here.
      * \param argc number of arguments in argv
      * \param argv the command line arguments provided to int main()
      */
      virtual void Init(int& argc,char**& argv) =0;
      /*! 
      * Initialize a communicator object. As many initializers (e.g. MPI, VTK, etc.) 
      * require the arguments to `int main()', they are passed through here. 
      * The data exchange pattern can be determined for final configuration if the mesh is passed in
      * \param argc number of arguments in argv
      * \param argv the command line arguments provided to int main()
      * \param mesh the local mesh for simulation
      */
      virtual void Init(int& argc,char**& argv,std::shared_ptr<PanNDE::Mesh> mesh) =0;
      /*!
      * get the id number of the process in the local communicator
      */
      virtual int getProcessId() =0;
      /*!
      * get the id number of the process in the global communicator
      */
      virtual int getGlobalId() =0;
      /*!
      * get the number of processes on the local communicator
      */
      virtual int getNumberOfProcesses() =0;
      /*!
      * get the number of processes on the global communicator
      */
      virtual int getNumberOfAllProcesses() =0;
      /*!
      * release commmunicator assets
      */
      virtual void Finalize() =0;

      /*!
      * use the local mesh to determine the halo exchange pattern. Not recommended approach, but
      * retained for backwards compatibility
      * \param mesh the local mesh from which the communication pattern can be determined
      */
      virtual void determineExchangePattern(std::shared_ptr<PanNDE::Mesh> mesh) =0;
      /*!
      * after having determined the exchange pattern, use that pattern to configure the halo exchange for
      * the specific field. Not recommended approach, but retained for backwards compatibility
      * \param keyname field name
      * \param field the field data to be registered for halo exchange
      */
      virtual void setupDataLinks(std::string keyname,std::shared_ptr<PanNDE::Field> field) =0;
      /*!
      * Use the field bundle to determine the exchange pattern and register all fields contained for data
      * exchange.
      * \param fields the field bundle to set up for halo exchange
      */
      virtual void setupDataLinks(std::shared_ptr<PanNDE::FieldBundle> fields) =0;
      
      /*!
      * Start asynchronous halo exchange on field by name
      * \param keyname field name
      */
      virtual void startHaloExchange(std::string keyname) =0;
      /*!
      * Wait for completion of halo exchange on field by name
      * \param keyname field name
      */
      virtual void waitUntilDone(std::string keyname) =0;

      //Send/receives: add types as appropriate. templates would be _nice_ but virtual
      //  templates are disallowed
      /*!
      * Send array to process
      * \param array the array to send
      * \param receiver process to receive the array
      */
      virtual void sendArray(std::shared_ptr<PanNDE::Array<double>> array,int receiver) =0;
      /*!
      * Send array to process
      * \param array the array to send
      * \param N the length of the array
      * \param receiver process to receive the array
      */
      virtual void sendArray(double* array,int N,int receiver) =0;
      /*!
      * Receive array from process
      * \param maker the array factory method required to build the array
      * \param sender process which sent the array
      */
      virtual std::shared_ptr<PanNDE::Array<double>> recvArray(
                              std::shared_ptr<PanNDE::ArrayFactory<double>> maker,int sender) =0;
      /*!
      * broadcast array by process
      * \param array the array to be either populated or sent
      * \param sender process which broadcasts the array
      */
      virtual void broadcastArray(std::shared_ptr<PanNDE::Array<double>>& array,int sender=0) =0;

      /*!
      * Send array to process
      * \param array the array to send
      * \param receiver process to receive the array
      */
      virtual void sendArray(std::shared_ptr<PanNDE::Array<int32_t>> array,int receiver) =0;
      /*!
      * Send array to process
      * \param array the array to send
      * \param N the length of the array
      * \param receiver process to receive the array
      */
      virtual void sendArray(int32_t* array,int N,int receiver) =0;
      /*!
      * Receive array from process
      * \param maker the array factory method required to build the array
      * \param sender process which sent the array
      */
      virtual std::shared_ptr<PanNDE::Array<int32_t>> recvArray(
                              std::shared_ptr<PanNDE::ArrayFactory<int32_t>> maker,int sender) =0;
      /*!
      * broadcast array by process
      * \param array the array to be either populated or sent
      * \param sender process which broadcasts the array
      */
      virtual void broadcastArray(std::shared_ptr<PanNDE::Array<int32_t>>& array,int sender=0) =0;

      /*!
      * Send array to process
      * \param array the array to send
      * \param receiver process to receive the array
      */
      virtual void sendArray(std::shared_ptr<PanNDE::Array<int64_t>> array,int receiver) =0;
      /*!
      * Send array to process
      * \param array the array to send
      * \param N the length of the array
      * \param receiver process to receive the array
      */
      virtual void sendArray(int64_t* array,int N,int receiver) =0;
      /*!
      * Receive array from process
      * \param maker the array factory method required to build the array
      * \param sender process which sent the array
      */
      virtual std::shared_ptr<PanNDE::Array<int64_t>> recvArray(
                              std::shared_ptr<PanNDE::ArrayFactory<int64_t>> maker,int sender) =0;
      /*!
      * broadcast array by process
      * \param array the array to be either populated or sent
      * \param sender process which broadcasts the array
      */
      virtual void broadcastArray(std::shared_ptr<PanNDE::Array<int64_t>>& array,int sender=0) =0;

      /*!
      * send value to process
      * \param value the value to send
      * \param receiver process to send value
      */
      virtual void sendValue(double value,int receiver) =0;
      /*!
      * receive value from process
      * \param value the address to which to write
      * \param sender process which sent the value
      */
      virtual void recvValue(double* value,int sender) =0;
      /*!
      * broadcast value to all processes
      * \param value the address to which to write (if receiving) or send (if broadcasting)
      * \param sender process which sent the value
      */
      virtual void broadcastValue(double* value,int sender=0) =0;
      /*!
      * all broadcast value for placement in array by process index
      * \param value the value to send
      * \param maker the array factory to build the synthesized array
      */
      virtual std::shared_ptr<PanNDE::Array<double>> allGatherValue(double value,
                              std::shared_ptr<PanNDE::ArrayFactory<double>> maker) =0;

      /*!
      * send value to process
      * \param value the value to send
      * \param receiver process to send value
      */
      virtual void sendValue(int32_t value,int receiver) =0;
      /*!
      * receive value from process
      * \param value the address to which to write
      * \param sender process which sent the value
      */
      virtual void recvValue(int32_t* value,int sender) =0;
      /*!
      * broadcast value to all processes
      * \param value the address to which to write (if receiving) or send (if broadcasting)
      * \param sender process which sent the value
      */
      virtual void broadcastValue(int32_t* value,int sender=0) =0;
      /*!
      * all broadcast value for placement in array by process index
      * \param value the value to send
      * \param maker the array factory to build the synthesized array
      */
      virtual std::shared_ptr<PanNDE::Array<int32_t>> allGatherValue(int32_t value,
                              std::shared_ptr<PanNDE::ArrayFactory<int32_t>> maker) =0;

      /*!
      * send value to process
      * \param value the value to send
      * \param receiver process to send value
      */
      virtual void sendValue(int64_t value,int receiver) =0;
      /*!
      * receive value from process
      * \param value the address to which to write
      * \param sender process which sent the value
      */
      virtual void recvValue(int64_t* value,int sender) =0;
      /*!
      * broadcast value to all processes
      * \param value the address to which to write (if receiving) or send (if broadcasting)
      * \param sender process which sent the value
      */
      virtual void broadcastValue(int64_t* value,int sender=0) =0;
      /*!
      * all broadcast value for placement in array by process index
      * \param value the value to send
      * \param maker the array factory to build the synthesized array
      */
      virtual std::shared_ptr<PanNDE::Array<int64_t>> allGatherValue(int64_t value,
                              std::shared_ptr<PanNDE::ArrayFactory<int64_t>> maker) =0;

      /*!
      * wait for all transmissions/receptions to complete
      */
      virtual void waitall() =0;
      /*!
      * wait for all processes to arrive
      */
      virtual void barrier() =0;
  };
};