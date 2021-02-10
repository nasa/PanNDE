/*! \headerfile Partitioner.hpp "modules/Controller/include/Partitioner.hpp"
* "Partitioner.hpp" contains the class definition for partitioning and distributing a mesh 
* and it's associated fields
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
  * Defines the required externally faced operations for partitioning and distributing a mesh 
  * and it's associated fields
  *
  */
  class Partitioner {
    public:
      /*!
      * partition a mesh
      * \param Nparts number of parts to decompose the mesh into
      * \param global_mesh the monolithic mesh to partition
      */
      virtual void partitionMesh(int Nparts,std::shared_ptr<PanNDE::Mesh> global_mesh) =0;
      
      /*!
      * distribute the mesh
      * \param partitionIds which partitions the local process is requesting
      * \param maker the mesh factory
      * \param mesh_source_rank process which has the global mesh and is distributing the mesh partitions
      */
      virtual std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> distributeMeshPartitions(
                                  std::shared_ptr<PanNDE::Array<int>> partitionIds,
                                  std::shared_ptr<PanNDE::MeshFactory> maker,
                                  int mesh_source_rank=0) =0;

      /*!
      * distribute a field
      * \param partitionIds which partitions the local process is requesting
      * \param maker the field factory
      * \param global_field the field that must be partitioned that is currently on the field_source_rank process
      * \param field_source_rank process which has the global field and is distributing the field partitions
      */
      virtual std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Field>>> distributeFieldPartitions(
                                  std::shared_ptr<PanNDE::Array<int>> partitionIds,
                                  std::shared_ptr<PanNDE::FieldFactory> maker,
                                  std::shared_ptr<PanNDE::Field> global_field=nullptr,
                                  int field_source_rank=0) =0;
  };
};
