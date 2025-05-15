/*! \headerfile Mesh.hpp "modules/PanNDE/include/Mesh.hpp"
* "Mesh.hpp" defines the fundamental spatial discretization interfaces for the PanNDE framework.
* It provides abstractions for computational meshes, supporting both serial and parallel
* simulations across various element types and partitioning strategies.
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
  * Defines the core interface for computational meshes in the PanNDE framework.
  *
  * A mesh represents the spatial discretization of a domain into cells for numerical 
  * simulations. This interface supports structured and unstructured meshes with
  * various element types (primarily hexahedral elements), and provides methods
  * for accessing geometry, topology, and partitioning information.
  *
  * Meshes form the foundation for fields (solution variables) and are essential
  * components for all physics models in the framework. Implementations must
  * carefully maintain connectivity relationships and support parallel computing
  * through mesh partitioning.
  *
  */
  class Mesh{
    public:
      //! Node structure representing mesh vertices
      struct Node{
        std::shared_ptr<PanNDE::Array<double>> coord;  //!< 3D coordinates [x,y,z]
        int64_t global_id;                             //!< Global identifier in partitioned meshes
        int32_t owner_partition;                       //!< Partition that owns this node
      };

      //! Cell structure representing mesh elements
      struct Cell{
        std::shared_ptr<PanNDE::Array<int32_t>> nodes; //!< Node indices defining the cell
        int64_t global_id;                             //!< Global identifier in partitioned meshes
        int32_t owner_partition;                       //!< Partition that owns this cell
      };
      
      /*!
      * Gets the total number of nodes in this mesh partition.
      *
      * This includes both locally owned nodes and ghost nodes from
      * neighboring partitions.
      *
      * \return int32_t Number of nodes
      */
      virtual int32_t nodeCount() =0;

      /*!
      * Gets the total number of cells in this mesh partition.
      *
      * This includes both locally owned cells and ghost cells from
      * neighboring partitions.
      *
      * \return int32_t Number of cells
      */
      virtual int32_t cellCount() =0;

      /*!
      * Gets the partition identifier for this mesh.
      *
      * For single-partition (serial) meshes, this is typically 0.
      * For parallel meshes, this is a unique identifier for the partition.
      *
      * \return int32_t Partition ID
      */
      virtual int32_t partitionId() =0;
      
      /*!
      * Gets the coordinates of a node using a pre-allocated array.
      *
      * \param node_id int32_t Local node index
      * \param coord double[3] Output array to receive the x,y,z coordinates
      */
      virtual void nodeCoordinate(int32_t node_id,double coord[3]) =0;

      /*!
      * Gets the coordinates of a node as a shared array.
      *
      * \param node_id int32_t Local node index
      * \return std::shared_ptr<PanNDE::Array<double>> Array containing x,y,z coordinates
      */
      virtual std::shared_ptr<PanNDE::Array<double>> nodeCoordinate(int32_t node_id) =0;

      /*!
      * Gets the node indices of a cell using a pre-allocated array.
      *
      * The nodes should be returned in a CGNS-compliant ordering. For hexahedral cells,
      * this means vertices are numbered in counter-clockwise order for each face,
      * starting with the bottom face and then the top face.
      *
      * \param cell_id int32_t Local cell index
      * \param nodes int32_t* Pre-allocated array to receive node indices
      */
      virtual void cell(int32_t cell_id,int32_t* nodes) =0;

      /*!
      * Gets the node indices of a cell as a shared array.
      *
      * The nodes should be returned in a CGNS-compliant ordering. For hexahedral cells,
      * this means vertices are numbered in counter-clockwise order for each face,
      * starting with the bottom face and then the top face.
      *
      * \param cell_id int32_t Local cell index
      * \return std::shared_ptr<PanNDE::Array<int32_t>> Array of node indices
      */
      virtual std::shared_ptr<PanNDE::Array<int32_t>> cell(int32_t cell_id) =0;
      
      /*!
      * Gets the cell indices connected to a node using a pre-allocated array.
      *
      * This provides the dual connectivity (inverse of the cell-to-nodes mapping).
      * For well-structured meshes, the ordering should be consistent with CGNS
      * winding direction from the node's perspective.
      *
      * \param node_id int32_t Local node index
      * \param cells int32_t* Pre-allocated array to receive cell indices
      */
      virtual void connectedCells(int32_t node_id,int32_t* cells)=0;

      /*!
      * Gets the cell indices connected to a node as a shared array.
      *
      * This provides the dual connectivity (inverse of the cell-to-nodes mapping).
      * For well-structured meshes, the ordering should be consistent with CGNS
      * winding direction from the node's perspective.
      *
      * \param node_id int32_t Local node index
      * \return std::shared_ptr<PanNDE::Array<int32_t>> Array of cell indices
      */
      virtual std::shared_ptr<PanNDE::Array<int32_t>> connectedCells(int32_t node_id) =0;

      /*!
      * Maps a local node index to its global node index.
      *
      * For partitioned meshes, this allows cross-referencing nodes between partitions.
      * For single-partition meshes, this typically returns the same value as the input.
      *
      * \param node_id int32_t Local node index
      * \return int64_t Global node index
      */
      virtual int64_t globalNodeId(int32_t node_id) =0;

      /*!
      * Gets the partition ID that owns a node.
      *
      * Ghost nodes (those replicated from neighboring partitions) will report
      * a different owner than the current partition.
      *
      * \param node_id int32_t Local node index
      * \return int32_t Owner partition ID
      */
      virtual int32_t nodeHomePartition(int32_t node_id) =0;

      /*!
      * Maps a local cell index to its global cell index.
      *
      * For partitioned meshes, this allows cross-referencing cells between partitions.
      * For single-partition meshes, this typically returns the same value as the input.
      *
      * \param cell_id int32_t Local cell index
      * \return int64_t Global cell index
      */
      virtual int64_t globalCellId(int32_t cell_id) =0;

      /*!
      * Gets the partition ID that owns a cell.
      *
      * Ghost cells (those replicated from neighboring partitions) will report
      * a different owner than the current partition.
      *
      * \param cell_id int32_t Local cell index
      * \return int32_t Owner partition ID
      */
      virtual int32_t cellHomePartition(int32_t cell_id) =0;

      /*!
      * Creates a copy of another mesh using a raw pointer.
      *
      * The implementation should perform a deep copy of all mesh data.
      *
      * \param other PanNDE::Mesh* The source mesh to copy
      */
      virtual void copy(PanNDE::Mesh* other) =0;

      /*!
      * Creates a copy of another mesh using a shared pointer.
      *
      * The implementation should perform a deep copy of all mesh data.
      *
      * \param other std::shared_ptr<PanNDE::Mesh> The source mesh to copy
      */
      virtual void copy(std::shared_ptr<PanNDE::Mesh> other) =0;

      //! Virtual destructor to ensure proper cleanup in derived classes
      virtual ~Mesh(){};
  };

  /*! \class MeshFactory Mesh.hpp "modules/PanNDE/include/Mesh.hpp"
  *
  * Factory interface for creating mesh objects and components.
  *
  * This interface provides methods to create nodes, cells, and complete meshes
  * with different memory management strategies. Implementations should focus on
  * efficient creation of mesh components for their target architecture and
  * maintain proper connectivity relationships.
  *
  */
  class MeshFactory{
    public:
      /*!
      * Creates a node with specified coordinates and attributes.
      *
      * \param coord double[3] Node coordinates in 3D space (x,y,z)
      * \param global_id int64_t Global identifier for the node
      * \param owner int32_t ID of the partition that owns this node (default: 0)
      * \return PanNDE::Mesh::Node A node structure with the specified properties
      */
      virtual PanNDE::Mesh::Node makeNode(double coord[3],
                                          int64_t global_id,
                                          int32_t owner=0) =0;

      /*!
      * Creates a cell with specified connectivity and attributes.
      *
      * \param nodes int32_t* Array of node indices that form the cell
      * \param Nnodes int32_t Number of nodes in the cell (e.g., 8 for hexahedra, 4 for tetrahedra)
      * \param global_id int64_t Global identifier for the cell
      * \param owner int32_t ID of the partition that owns this cell (default: 0)
      * \return PanNDE::Mesh::Cell A cell structure with the specified properties
      */
      virtual PanNDE::Mesh::Cell makeCell(int32_t* nodes,int32_t Nnodes,
                                          int64_t global_id,
                                          int32_t owner=0) =0;

      
      /*!
      * Creates a mesh from arrays of nodes and cells with shared pointer management.
      *
      * This is the preferred method for creating meshes as it provides automatic
      * memory management.
      *
      * \param nodes PanNDE::Mesh::Node* Array of node structures
      * \param Nnodes int64_t Number of nodes
      * \param cells PanNDE::Mesh::Cell* Array of cell structures
      * \param Ncells int64_t Number of cells
      * \param partition int32_t Partition ID for the mesh (default: 0)
      * \return std::shared_ptr<PanNDE::Mesh> A shared pointer to the newly created mesh
      */
      virtual std::shared_ptr<PanNDE::Mesh> makeManagedMesh(
                      PanNDE::Mesh::Node* nodes,int64_t Nnodes,
                      PanNDE::Mesh::Cell* cells,int64_t Ncells,int32_t partition=0) =0;

      /*!
      * Creates a mesh from PanNDE arrays of nodes and cells with shared pointer management.
      *
      * This version uses PanNDE arrays instead of raw pointers, which can provide
      * better safety and flexibility.
      *
      * \param nodes std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Node>> Array of nodes
      * \param cells std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Cell>> Array of cells
      * \param partition int32_t Partition ID for the mesh (default: 0)
      * \return std::shared_ptr<PanNDE::Mesh> A shared pointer to the newly created mesh
      */
      virtual std::shared_ptr<PanNDE::Mesh> makeManagedMesh(
                      std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Node>> nodes,
                      std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Cell>> cells,int32_t partition=0) =0;

      /*!
      * Creates a mesh from arrays of nodes and cells with manual memory management.
      *
      * This method is not recommended for general use and requires manual deletion
      * using deleteMesh(). It's provided primarily for compatibility with C-style interfaces.
      *
      * \param nodes PanNDE::Mesh::Node* Array of node structures
      * \param Nnodes int64_t Number of nodes
      * \param cells PanNDE::Mesh::Cell* Array of cell structures
      * \param Ncells int64_t Number of cells
      * \param partition int32_t Partition ID for the mesh (default: 0)
      * \return PanNDE::Mesh* A raw pointer to the newly created mesh
      */
      virtual PanNDE::Mesh* newMesh(
                      PanNDE::Mesh::Node* nodes,int64_t Nnodes,
                      PanNDE::Mesh::Cell* cells,int64_t Ncells,int32_t partition=0) =0;

      /*!
      * Creates a mesh from PanNDE arrays of nodes and cells with manual memory management.
      *
      * This method is not recommended for general use and requires manual deletion
      * using deleteMesh(). It's provided primarily for compatibility with C-style interfaces.
      *
      * \param nodes std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Node>> Array of nodes
      * \param cells std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Cell>> Array of cells
      * \param partition int32_t Partition ID for the mesh (default: 0)
      * \return PanNDE::Mesh* A raw pointer to the newly created mesh
      */
      virtual PanNDE::Mesh* newMesh(
                      std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Node>> nodes,
                      std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Cell>> cells,int32_t partition=0) =0;

      /*!
      * Deletes a mesh that was created using newMesh().
      *
      * This method must be called for any mesh created with newMesh() to prevent memory leaks.
      *
      * \param mesh PanNDE::Mesh* Pointer to the mesh to delete
      */
      virtual void deleteMesh(PanNDE::Mesh* mesh) =0;

      /*!
      * Creates an array to store multiple meshes.
      *
      * This is useful for multi-domain simulations or when working with
      * multiple related meshes (e.g., a sequence of evolving meshes).
      *
      * \return std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> An array of mesh pointers
      */
      virtual std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> makeManagedMeshArray() =0;
  };
};