/*! \headerfile HexMesh.hpp "modules/HostData/include/HexMesh.hpp"
* "HexMesh.hpp" contains class implementations for hexahedral (voxel) mesh representation
* and mesh creation. It provides host-based implementations of the PanNDE mesh interfaces 
* for creating, storing, and accessing structured hexahedral meshes.
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
#include <array>

//required interfaces
#include "Array.hpp"
#include "Mesh.hpp"

//same module inclusions
#include "HostArray.hpp"

namespace HostData {
  /*! \class HexMesh HexMesh.hpp "modules/HostData/include/HexMesh.hpp"
  * 
  * Implements a host-based hexahedral (voxel) mesh data structure.
  * This class provides a concrete implementation of the PanNDE::Mesh interface 
  * for structured hexahedral meshes, supporting partitioned meshes for parallel 
  * computation. Mesh data includes node coordinates, cell connectivity, partition
  * ownership information, and global/local indexing mappings.
  *
  */
  class HexMesh : public PanNDE::Mesh {
    public:
      /*!
      * Gets the total number of nodes in this mesh partition.
      * \return int32_t The number of nodes
      */
      int32_t nodeCount()override{return node_owner.size();};
      
      /*!
      * Gets the total number of cells in this mesh partition.
      * \return int32_t The number of cells
      */
      int32_t cellCount()override{return hex_owner.size();};
      
      /*!
      * Gets the partition identifier for this mesh.
      * \return int32_t The partition ID
      */
      int32_t partitionId()override{return mypartition;};
      
      /*!
      * Gets the coordinates of a specific node.
      * \param node_id int32_t Index of the node
      * \param coord double[3] Output array to hold X, Y, Z coordinates
      */
      void nodeCoordinate(int32_t node_id,double coord[3])override{
        coord[0]=X.at(node_id);coord[1]=Y.at(node_id);coord[2]=Z.at(node_id);
      };
      
      /*!
      * Gets the coordinates of a specific node as an array.
      * \param node_id int32_t Index of the node
      * \return std::shared_ptr<PanNDE::Array<double>> Array containing X, Y, Z coordinates
      */
      std::shared_ptr<PanNDE::Array<double>> nodeCoordinate(int32_t node_id)override{
        auto coord=darraymfg->makeManagedArray();
        coord->resize(0);coord->resize(3);
        coord->at(0)=X.at(node_id);coord->at(1)=Y.at(node_id);coord->at(2)=Z.at(node_id);
        return std::move(coord);
      };

      /*!
      * Gets the node indices defining a specific hexahedral cell.
      * \param cell_id int32_t Index of the cell
      * \param nodes int32_t* Output array to hold the 8 node indices in CGNS order
      */
      void cell(int32_t cell_id,int32_t* nodes)override{
        for(int kb=0;kb<8;kb++){
          nodes[kb]=hexes.at(cell_id).at(kb);
        };
      };
      
      /*!
      * Gets the node indices defining a specific hexahedral cell as an array.
      * \param cell_id int32_t Index of the cell
      * \return std::shared_ptr<PanNDE::Array<int32_t>> Array containing 8 node indices in CGNS order
      */
      std::shared_ptr<PanNDE::Array<int32_t>> cell(int32_t cell_id)override{
        auto nodes=iarraymfg->makeManagedArray();
        nodes->resize(0);nodes->resize(8);
        for(int kb=0;kb<8;kb++){nodes->at(kb)=hexes.at(cell_id).at(kb);};
        return std::move(nodes);
      };

      /*!
      * Gets the cell indices connected to a specific node.
      * \param node_id int32_t Index of the node
      * \param cells int32_t* Output array to hold up to 8 cell indices; -1 indicates no cell
      */
      void connectedCells(int32_t node_id,int32_t* cells)override{
        for(int kb=0;kb<8;kb++){cells[kb]=node_connected_cells.at(node_id).at(kb);};
      };
      
      /*!
      * Gets the cell indices connected to a specific node as an array.
      * \param node_id int32_t Index of the node
      * \return std::shared_ptr<PanNDE::Array<int32_t>> Array containing up to 8 cell indices; -1 indicates no cell
      */
      std::shared_ptr<PanNDE::Array<int32_t>> connectedCells(int32_t node_id)override{
        auto box=iarraymfg->makeManagedArray();
        box->resize(8);
        for(int kb=0;kb<8;kb++){box->at(kb)=node_connected_cells.at(node_id).at(kb);};
        return std::move(box);
      };

      /*!
      * Gets the global node index for a local node ID.
      * \param node_id int32_t Local index of the node
      * \return int64_t Global (domain-wide) node index
      */
      int64_t globalNodeId(int32_t node_id)override{return node_global_id.at(node_id);};
      
      /*!
      * Gets the home partition ID for a node.
      * \param node_id int32_t Local index of the node
      * \return int32_t ID of the partition that owns this node
      */
      int32_t nodeHomePartition(int32_t node_id)override{return node_owner.at(node_id);};

      /*!
      * Gets the global cell index for a local cell ID.
      * \param cell_id int32_t Local index of the cell
      * \return int64_t Global (domain-wide) cell index
      */
      int64_t globalCellId(int32_t cell_id)override{return hex_global_id.at(cell_id);};
      
      /*!
      * Gets the home partition ID for a cell.
      * \param cell_id int32_t Local index of the cell
      * \return int32_t ID of the partition that owns this cell
      */
      int32_t cellHomePartition(int32_t cell_id)override{return hex_owner.at(cell_id);};

      /*!
      * Creates a copy of another mesh.
      * \param other PanNDE::Mesh* Pointer to the source mesh
      * \throw std::out_of_range If node indices in cells are invalid
      */
      void copy(PanNDE::Mesh* other)override{
        mypartition=other->partitionId();
        resizeVectors(other->nodeCount(),other->cellCount());
        for(int k=0;k<other->nodeCount();k++){
          setNode(k,other->nodeCoordinate(k),
                  other->nodeHomePartition(k),
                  other->globalNodeId(k));
        };
        for(int k=0;k<other->cellCount();k++){
          setCell(k,other->cell(k),other->cellHomePartition(k),other->globalCellId(k));
        };
        buildLinks();
      };
      
      /*!
      * Creates a copy of another mesh.
      * \param other std::shared_ptr<PanNDE::Mesh> Shared pointer to the source mesh
      * \throw std::out_of_range If node indices in cells are invalid
      */
      void copy(std::shared_ptr<PanNDE::Mesh> other)override{copy(other.get());};

      /*!
      * Default constructor creating an empty mesh.
      */
      HexMesh(){};
      
      /*!
      * Constructs a mesh from arrays of nodes and cells.
      * \param nodes std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Node>> Array of node definitions
      * \param cells std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Cell>> Array of cell definitions
      * \param partition int32_t Partition ID for this mesh (default: 0)
      * \throw std::out_of_range If node indices in cells are invalid
      */
      HexMesh(std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Node>> nodes,
              std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Cell>> cells,int32_t partition=0){
        mypartition=partition;
        resizeVectors(nodes->size(),cells->size());
        double coord[3];
        for(int k=0;k<nodes->size();k++){
          setNode(k,nodes->at(k).coord,nodes->at(k).owner_partition,nodes->at(k).global_id);
        };
        for(int k=0;k<cells->size();k++){
          setCell(k,cells->at(k).nodes,cells->at(k).owner_partition,cells->at(k).global_id);
        };
        buildLinks();
      };
      
    private:
      /*!
      * Sets the node indices for a specific cell.
      * \param cell_id int Index of the cell to modify
      * \param box std::shared_ptr<PanNDE::Array<int32_t>> Array of 8 node indices
      * \param owner int32_t Partition ID that owns this cell
      * \param g_id int32_t Global ID for this cell
      * \throw std::out_of_range If any node index is invalid
      */
      void setCell(int cell_id,std::shared_ptr<PanNDE::Array<int32_t>> box,
                   int32_t owner,int32_t g_id){
        for(int kb=0;kb<8;kb++){
          if(node_owner.size()<=box->at(kb)){
            throw std::out_of_range(std::string("invalid node id: Node count: ")+
                                    std::to_string(node_owner.size())+
                                    std::string(", provided index: ")+
                                    std::to_string(box->at(kb)));
          };
          hexes.at(cell_id).at(kb)=box->at(kb);
        };
        hex_owner.at(cell_id)=owner;
        hex_global_id.at(cell_id)=g_id;
      };
      
      /*!
      * Sets the coordinates and metadata for a specific node.
      * \param node_id int Index of the node to modify
      * \param coord std::shared_ptr<PanNDE::Array<double>> Array of 3 coordinate values
      * \param owner int32_t Partition ID that owns this node
      * \param g_id int32_t Global ID for this node
      */
      void setNode(int node_id,std::shared_ptr<PanNDE::Array<double>> coord,
                   int32_t owner,int32_t g_id){
        X.at(node_id)=coord->at(0);
        Y.at(node_id)=coord->at(1);
        Z.at(node_id)=coord->at(2);
        node_owner.at(node_id)=owner;
        node_global_id.at(node_id)=g_id;
      };
      
      /*!
      * Resizes all internal vectors to accommodate the specified number of nodes and cells.
      * \param Nnodes int Number of nodes to allocate space for
      * \param Ncells int Number of cells to allocate space for
      */
      void resizeVectors(int Nnodes,int Ncells){
        X.resize(0);X.resize(Nnodes);
        Y.resize(0);Y.resize(Nnodes);
        Z.resize(0);Z.resize(Nnodes);
        node_owner.resize(0);node_owner.resize(Nnodes);
        node_global_id.resize(0);node_global_id.resize(Nnodes);
        hexes.resize(0);hexes.resize(Ncells);
        hex_owner.resize(0);hex_owner.resize(Ncells);
        hex_global_id.resize(0);hex_global_id.resize(Ncells);
      };
      
      /*!
      * Builds the connectivity mapping between nodes and cells.
      * For each node, this determines which cells include that node as a vertex.
      */
      void buildLinks(){
        node_connected_cells.resize(0);node_connected_cells.resize(node_owner.size());
        for(int32_t k=0;k<node_owner.size();k++){node_connected_cells.at(k).fill(-1);};
        for(int32_t k=0;k<hexes.size();k++){
          node_connected_cells.at(hexes.at(k).at(0)).at(6)=k;
          node_connected_cells.at(hexes.at(k).at(1)).at(7)=k;
          node_connected_cells.at(hexes.at(k).at(2)).at(4)=k;
          node_connected_cells.at(hexes.at(k).at(3)).at(5)=k;
          node_connected_cells.at(hexes.at(k).at(4)).at(2)=k;
          node_connected_cells.at(hexes.at(k).at(5)).at(3)=k;
          node_connected_cells.at(hexes.at(k).at(6)).at(0)=k;
          node_connected_cells.at(hexes.at(k).at(7)).at(1)=k;
        };
      };

      //! Factory for creating arrays of double values
      std::shared_ptr<HostData::HostArrayFactory<double>> darraymfg=
              std::make_shared<HostData::HostArrayFactory<double>>(HostData::HostArrayFactory<double>());
              
      //! Factory for creating arrays of integer values
      std::shared_ptr<HostData::HostArrayFactory<int32_t>> iarraymfg=
              std::make_shared<HostData::HostArrayFactory<int32_t>>(HostData::HostArrayFactory<int32_t>());

      //! X coordinates for all nodes
      std::vector<double> X;
      
      //! Y coordinates for all nodes
      std::vector<double> Y;
      
      //! Z coordinates for all nodes
      std::vector<double> Z;
      
      //! Partition ownership information for each node
      std::vector<int32_t> node_owner;
      
      //! Global ID mapping for each node
      std::vector<int64_t> node_global_id;
      
      //! For each node, stores up to 8 connected cells (with -1 indicating no connection)
      std::vector<std::array<int32_t,8>> node_connected_cells;

      //! For each cell, stores the 8 nodes that form its vertices in CGNS order
      std::vector<std::array<int32_t,8>> hexes;
      
      //! Partition ownership information for each cell
      std::vector<int32_t> hex_owner;
      
      //! Global ID mapping for each cell
      std::vector<int64_t> hex_global_id;

      //! Partition ID of this mesh instance
      int32_t mypartition;
  };

  /*! \class HexMeshFactory HexMesh.hpp "modules/HostData/include/HexMesh.hpp"
  * 
  * Factory class for creating hexahedral mesh objects and components.
  * Implements the PanNDE::MeshFactory interface to provide methods for creating
  * nodes, cells, and complete mesh objects for hexahedral meshes.
  *
  */
  class HexMeshFactory : public PanNDE::MeshFactory {
    public:
      /*!
      * Creates a shared pointer to a new HexMeshFactory instance.
      * \return std::shared_ptr<HostData::HexMeshFactory> Shared pointer to the new factory
      */
      static
      std::shared_ptr<HostData::HexMeshFactory> makeShared(){
        auto maker=std::make_shared<HostData::HexMeshFactory>(HostData::HexMeshFactory());
        return std::move(maker);
      };
      
      /*!
      * Creates a mesh node structure with the specified properties.
      * \param coord[3] double Array of X, Y, Z coordinates for the node
      * \param global_id int64_t Global ID for the node
      * \param owner int32_t Partition ID that owns this node (default: 0)
      * \return PanNDE::Mesh::Node Node structure with the specified properties
      */
      PanNDE::Mesh::Node makeNode(double coord[3],
                                  int64_t global_id,
                                  int32_t owner=0)override{
        PanNDE::Mesh::Node node;
        auto darraymfg=std::make_shared<HostData::HostArrayFactory<double>>(HostData::HostArrayFactory<double>());
        node.coord=darraymfg->makeManagedArray();
        node.coord->resize(3);
        for(int k=0;k<3;k++){node.coord->at(k)=coord[k];};
        node.global_id=global_id;
        node.owner_partition=owner;
        return node;
      };
      
      /*!
      * Creates a mesh cell structure with the specified properties.
      * \param nodes int32_t* Array of node indices forming the cell vertices
      * \param Nnodes int32_t Number of nodes in the cell (typically 8 for hexahedral cells)
      * \param global_id int64_t Global ID for the cell
      * \param owner int32_t Partition ID that owns this cell (default: 0)
      * \return PanNDE::Mesh::Cell Cell structure with the specified properties
      */
      PanNDE::Mesh::Cell makeCell(int32_t* nodes,int32_t Nnodes,
                                  int64_t global_id,
                                  int32_t owner=0)override{
        PanNDE::Mesh::Cell cell;
        auto iarraymfg=std::make_shared<HostData::HostArrayFactory<int32_t>>(HostData::HostArrayFactory<int32_t>());
        cell.nodes=iarraymfg->makeManagedArray();
        cell.nodes->resize(Nnodes);
        for(int k=0;k<Nnodes;k++){cell.nodes->at(k)=nodes[k];};
        cell.global_id=global_id;
        cell.owner_partition=owner;
        return cell;
      };

      /*!
      * Creates a mesh from raw node and cell arrays.
      * \param nodes PanNDE::Mesh::Node* Array of node structures
      * \param Nnodes int64_t Number of nodes in the array
      * \param cells PanNDE::Mesh::Cell* Array of cell structures
      * \param Ncells int64_t Number of cells in the array
      * \param partition int32_t Partition ID for this mesh (default: 0)
      * \return std::shared_ptr<PanNDE::Mesh> Shared pointer to the created mesh
      */
      std::shared_ptr<PanNDE::Mesh> makeManagedMesh(
                      PanNDE::Mesh::Node* nodes,int64_t Nnodes,
                      PanNDE::Mesh::Cell* cells,int64_t Ncells,int32_t partition=0)override{
        auto narray=makeArray<PanNDE::Mesh::Node>(nodes,Nnodes);
        auto carray=makeArray<PanNDE::Mesh::Cell>(cells,Ncells);
        return std::make_shared<HostData::HexMesh>(HostData::HexMesh(narray,carray,partition));
      };
      
      /*!
      * Creates a mesh from node and cell arrays.
      * \param nodes std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Node>> Array of node structures
      * \param cells std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Cell>> Array of cell structures
      * \param partition int32_t Partition ID for this mesh (default: 0)
      * \return std::shared_ptr<PanNDE::Mesh> Shared pointer to the created mesh
      */
      std::shared_ptr<PanNDE::Mesh> makeManagedMesh(
                      std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Node>> nodes,
                      std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Cell>> cells,int32_t partition=0)override{
        return std::make_shared<HostData::HexMesh>(HostData::HexMesh(nodes,cells,partition));
      };
      
      /*!
      * Creates a raw pointer to a mesh object (caller takes ownership).
      * \param nodes PanNDE::Mesh::Node* Array of node structures
      * \param Nnodes int64_t Number of nodes in the array
      * \param cells PanNDE::Mesh::Cell* Array of cell structures
      * \param Ncells int64_t Number of cells in the array
      * \param partition int32_t Partition ID for this mesh (default: 0)
      * \return PanNDE::Mesh* Pointer to the created mesh
      * \warning Use with caution - caller is responsible for memory management
      */
      PanNDE::Mesh* newMesh(
                      PanNDE::Mesh::Node* nodes,int64_t Nnodes,
                      PanNDE::Mesh::Cell* cells,int64_t Ncells,int32_t partition=0)override{
        auto narray=makeArray<PanNDE::Mesh::Node>(nodes,Nnodes);
        auto carray=makeArray<PanNDE::Mesh::Cell>(cells,Ncells);
        PanNDE::Mesh* mesh=new HostData::HexMesh(narray,carray,partition);
        return mesh;
      };
      
      /*!
      * Creates a raw pointer to a mesh object (caller takes ownership).
      * \param nodes std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Node>> Array of node structures
      * \param cells std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Cell>> Array of cell structures
      * \param partition int32_t Partition ID for this mesh (default: 0)
      * \return PanNDE::Mesh* Pointer to the created mesh
      * \warning Use with caution - caller is responsible for memory management
      */
      PanNDE::Mesh* newMesh(
                      std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Node>> nodes,
                      std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Cell>> cells,int32_t partition=0)override{
        PanNDE::Mesh* mesh=new HostData::HexMesh(nodes,cells,partition);
        return mesh;
      };

      /*!
      * Deletes a mesh created using newMesh().
      * \param mesh PanNDE::Mesh* Pointer to the mesh to delete
      */
      void deleteMesh(PanNDE::Mesh* mesh)override{delete mesh;};

      /*!
      * Creates an empty array to store mesh objects.
      * \return std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> Shared pointer to the new array
      */
      std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> makeManagedMeshArray()override{
        return std::move(makeArray<std::shared_ptr<PanNDE::Mesh>>());
      };
      
    private:
      /*!
      * Helper method to create an array from raw data.
      * \param raw T* Pointer to raw array of elements
      * \param N int64_t Number of elements
      * \return std::shared_ptr<PanNDE::Array<T>> New array containing copies of the elements
      */
      template<typename T>
      std::shared_ptr<PanNDE::Array<T>> makeArray(T* raw,int64_t N){
        auto array=makeArray<T>();
        array->resize(N);
        for(int64_t k=0;k<N;k++){array->at(k)=raw[k];};
        return std::move(array);
      };
      
      /*!
      * Helper method to create an empty array.
      * \return std::shared_ptr<PanNDE::Array<T>> New empty array
      */
      template<typename T>
      std::shared_ptr<PanNDE::Array<T>> makeArray(){
        auto arraymfg=std::make_shared<HostData::HostArrayFactory<T>>(HostData::HostArrayFactory<T>());
        auto array=arraymfg->makeManagedArray();
        return std::move(array);
      };
  };
};