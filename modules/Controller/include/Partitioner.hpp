/*! \headerfile Partitioner.hpp "modules/Controller/include/Partitioner.hpp"
* "Partitioner.hpp" contains the abstract interface for domain decomposition and
* distribution of meshes and fields across multiple processes for parallel computing.
* It defines methods for partitioning global meshes and distributing the resulting
* partitions and their associated field data.
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

#include "Communicator.hpp"

namespace Controller {
  /*! \class Partitioner Partitioner.hpp "modules/Controller/include/Partitioner.hpp"
  * 
  * Defines an abstract interface for domain decomposition in parallel computations.
  * 
  * The Partitioner interface provides methods for dividing a global mesh into
  * smaller partitions for distributed parallel computing. It supports both the
  * partitioning process (dividing the computational domain) and the distribution
  * process (sending mesh and field partitions to appropriate processes).
  * 
  * Concrete implementations of this interface handle specific partitioning algorithms
  * and communication protocols to efficiently distribute the workload across multiple
  * processes while minimizing communication overhead.
  *
  */
  class Partitioner {
    public:
      /*!
      * Partitions a global mesh into multiple parts.
      * This method divides a monolithic mesh into the specified number of partitions,
      * typically using graph partitioning algorithms to minimize inter-partition
      * communication while balancing computational load.
      * 
      * \param Nparts int Number of partitions to create
      * \param global_mesh std::shared_ptr<PanNDE::Mesh> The monolithic mesh to partition
      */
      virtual void partitionMesh(int Nparts,std::shared_ptr<PanNDE::Mesh> global_mesh) =0;
      
      /*!
      * Distributes partitioned mesh objects to requesting processes.
      * This method sends the requested mesh partitions from the source process
      * to the calling process. Each process can request multiple partitions.
      * 
      * \param partitionIds std::shared_ptr<PanNDE::Array<int>> IDs of partitions requested by this process
      * \param maker std::shared_ptr<PanNDE::MeshFactory> Factory for creating mesh objects on the receiving end
      * \param mesh_source_rank int Rank of the process holding the partitioned global mesh (default: 0)
      * \return std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> Array of distributed mesh partitions
      */
      virtual std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> distributeMeshPartitions(
                                  std::shared_ptr<PanNDE::Array<int>> partitionIds,
                                  std::shared_ptr<PanNDE::MeshFactory> maker,
                                  int mesh_source_rank=0) =0;

      /*!
      * Distributes field data corresponding to mesh partitions.
      * This method sends the field data associated with requested mesh partitions
      * from the source process to the calling process, ensuring field values are
      * properly mapped to the partitioned meshes.
      * 
      * \param partitionIds std::shared_ptr<PanNDE::Array<int>> IDs of partitions requested by this process
      * \param maker std::shared_ptr<PanNDE::FieldFactory> Factory for creating field objects on the receiving end
      * \param global_field std::shared_ptr<PanNDE::Field> The global field to partition (if nullptr, assumed to exist on source_rank)
      * \param field_source_rank int Rank of the process holding the global field (default: 0)
      * \return std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Field>>> Array of distributed field partitions
      */
      virtual std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Field>>> distributeFieldPartitions(
                                  std::shared_ptr<PanNDE::Array<int>> partitionIds,
                                  std::shared_ptr<PanNDE::FieldFactory> maker,
                                  std::shared_ptr<PanNDE::Field> global_field=nullptr,
                                  int field_source_rank=0) =0;
  };
};