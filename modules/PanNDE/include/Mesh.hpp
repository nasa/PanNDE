/*! \headerfile Mesh.hpp "modules/PanNDE/include/Mesh.hpp"
* "Mesh.hpp" contains class definitions for mesh and mesh factory objects. 
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

#include <cstdint>
#include <memory>

#include "Array.hpp"

namespace PanNDE{
  /*! \class Mesh Mesh.hpp "modules/PanNDE/include/Mesh.hpp"
  *
  * Defines the methods required to interact with a Mesh 
  *
  */
  class Mesh{
    public:
      //mesh component definitions
      struct Node{
        std::shared_ptr<PanNDE::Array<double>> coord;
        int64_t global_id;
        int32_t owner_partition;
      };
      struct Cell{
        std::shared_ptr<PanNDE::Array<int32_t>> nodes;
        int64_t global_id;
        int32_t owner_partition;
      };
      
      /*!
      * get the number of nodes in the mesh
      */
      virtual int32_t nodeCount() =0;
      /*!
      * get the number of cells in the mesh
      */
      virtual int32_t cellCount() =0;
      /*!
      * get the partition identifier for this mesh
      */
      virtual int32_t partitionId() =0;
      
      /*!
      * get the nodal coordinates
      */
      virtual void nodeCoordinate(int32_t node_id,double coord[3]) =0;
      /*!
      * get the nodal coordinates
      */
      virtual std::shared_ptr<PanNDE::Array<double>> nodeCoordinate(int32_t node_id) =0;

      /*!
      * get the nodes which define the corners of the cell. A well constructed mesh implementation 
      * should return the nodes in a CGNS compliant winding scheme
      */
      virtual void cell(int32_t cell_id,int32_t* nodes) =0;
      /*!
      * get the nodes which define the corners of the cell. A well constructed mesh implementation 
      * should return the node indices in a CGNS compliant winding scheme
      */
      virtual std::shared_ptr<PanNDE::Array<int32_t>> cell(int32_t cell_id) =0;
      
      /*!
      * get the cells to which the node serves as corners. A well constructed mesh implementation 
      * should return the cell indices in a CGNS compliant winding scheme
      */
      virtual void connectedCells(int32_t node_id,int32_t* cells)=0;
      /*!
      * get the cells to which the node serves as corners. A well constructed mesh implementation 
      * should return the cell indices in a CGNS compliant winding scheme
      */
      virtual std::shared_ptr<PanNDE::Array<int32_t>> connectedCells(int32_t node_id) =0;

      /*!
      * get the global node index for a given local index (for a partitioned mesh)
      */
      virtual int64_t globalNodeId(int32_t node_id) =0;
      /*!
      * get the partition which owns the node (i.e., updates values at this node)
      */
      virtual int32_t nodeHomePartition(int32_t node_id) =0;

      /*!
      * get the global cell index for a given local index (for a partitioned mesh)
      */
      virtual int64_t globalCellId(int32_t cell_id) =0;
      /*!
      * get the partition which owns the cell (i.e., updates values at this cell)
      */
      virtual int32_t cellHomePartition(int32_t cell_id) =0;

      /*!
      * create a copy of another mesh
      * \param other the source mesh for the copy
      */
      virtual void copy(PanNDE::Mesh* other) =0;
      /*!
      * create a copy of another mesh
      * \param other the source mesh for the copy
      */
      virtual void copy(std::shared_ptr<PanNDE::Mesh> other) =0;
  };

  /*! \class MeshFactory Mesh.hpp "modules/PanNDE/include/Mesh.hpp"
  *
  * Defines a factory class to create a Mesh object
  *
  */
  class MeshFactory{
    public:
      /*!
      * make a node
      * \param coord[3] nodal cartesian coordinates
      * \param global_id the global identifier of the node
      * \param owner which mesh partition the node belongs to
      */
      virtual PanNDE::Mesh::Node makeNode(double coord[3],
                                          int64_t global_id,
                                          int32_t owner=0) =0;
      /*!
      * make a cell
      * \param nodes nodal indicies which form the corners of the cell
      * \param Nnodes number of nodes which define the cell (e.g. 8 for hexes, 4 for tets)
      * \param global_id the global identifier of the cell
      * \param owner which mesh partition the cell belongs to
      */
      virtual PanNDE::Mesh::Cell makeCell(int32_t* nodes,int32_t Nnodes,
                                          int64_t global_id,
                                          int32_t owner=0) =0;

      
      /*! 
      * create a shared mesh
      * \param nodes a C-style array (not preferred) of Node structs
      * \param Nnodes the length of the nodes array
      * \param cells a C-style array (not preferred) of cell structs
      * \param Ncells the length of the cells array
      * \param partition the partition identifier of the mesh to be created
      */
      virtual std::shared_ptr<PanNDE::Mesh> makeManagedMesh(
                      PanNDE::Mesh::Node* nodes,int64_t Nnodes,
                      PanNDE::Mesh::Cell* cells,int64_t Ncells,int32_t partition=0) =0;
      /*! 
      * create a shared mesh
      * \param nodes a PanNDE array of Node structs
      * \param cells a PanNDE array of cell structs
      * \param partition the partition identifier of the mesh to be created
      */
      virtual std::shared_ptr<PanNDE::Mesh> makeManagedMesh(
                      std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Node>> nodes,
                      std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Cell>> cells,int32_t partition=0) =0;
      /*! 
      * create a mesh on a raw pointer. not recommended, call deleteMesh() if used.
      * \param nodes a C-style array (not preferred) of Node structs
      * \param Nnodes the length of the nodes array
      * \param cells a C-style array (not preferred) of cell structs
      * \param Ncells the length of the cells array
      * \param partition the partition identifier of the mesh to be created
      */
      virtual PanNDE::Mesh* newMesh(
                      PanNDE::Mesh::Node* nodes,int64_t Nnodes,
                      PanNDE::Mesh::Cell* cells,int64_t Ncells,int32_t partition=0) =0;
      /*! 
      * create a mesh on a raw pointer. not recommended, call deleteMesh() if used.
      * \param nodes a PanNDE array of Node structs
      * \param cells a PanNDE array of cell structs
      * \param partition the partition identifier of the mesh to be created
      */
      virtual PanNDE::Mesh* newMesh(
                      std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Node>> nodes,
                      std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Cell>> cells,int32_t partition=0) =0;
      /*!
      * delete meshes created using newMesh()
      * \param mesh mesh to be deleted
      */
      virtual void deleteMesh(PanNDE::Mesh* mesh) =0;

      /*! 
      * create a mesh array
      */
      virtual std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> makeManagedMeshArray() =0;
  };
};