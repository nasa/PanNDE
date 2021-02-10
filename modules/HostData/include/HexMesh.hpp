/*! \headerfile HexMesh.hpp "modules/HostData/include/HexMesh.hpp"
* "HexMesh.hpp" contains class definitions for voxel mesh and voxel mesh factory objects. 
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

//required interfaces
#include "Array.hpp"
#include "Mesh.hpp"

//same module inclusions
#include "HostArray.hpp"

namespace HostData {
  /*! \class HexMesh HexMesh.hpp "modules/HostData/include/HexMesh.hpp"
  * Implements PanNDE::Mesh 
  *
  */
  class HexMesh : public PanNDE::Mesh {
    public:
      /*!
      * get the number of nodes in the mesh
      */
      int32_t nodeCount()override{return node_owner.size();};
      /*!
      * get the number of cells in the mesh
      */
      int32_t cellCount()override{return hex_owner.size();};
      /*!
      * get the partition identifier for this mesh
      */
      int32_t partitionId()override{return mypartition;};
      
      /*!
      * get the nodal coordinates
      */
      void nodeCoordinate(int32_t node_id,double coord[3])override{
        coord[0]=X.at(node_id);coord[1]=Y.at(node_id);coord[2]=Z.at(node_id);
      };
      /*!
      * get the nodal coordinates
      */
      std::shared_ptr<PanNDE::Array<double>> nodeCoordinate(int32_t node_id)override{
        auto coord=darraymfg->makeManagedArray();
        coord->resize(0);coord->resize(3);
        coord->at(0)=X.at(node_id);coord->at(1)=Y.at(node_id);coord->at(2)=Z.at(node_id);
        return std::move(coord);
      };

      /*!
      * get the nodes which define the corners of the cell. A well constructed mesh implementation 
      * should return the nodes in a CGNS compliant winding scheme
      */
      void cell(int32_t cell_id,int32_t* nodes)override{
        for(int kb=0;kb<8;kb++){
          nodes[kb]=hexes.at(cell_id).at(kb);
        };
      };
      /*!
      * get the nodes which define the corners of the cell. Returns the node indices in a 
      * CGNS compliant winding scheme
      */
      std::shared_ptr<PanNDE::Array<int32_t>> cell(int32_t cell_id)override{
        auto nodes=iarraymfg->makeManagedArray();
        nodes->resize(0);nodes->resize(8);
        for(int kb=0;kb<8;kb++){nodes->at(kb)=hexes.at(cell_id).at(kb);};
        return std::move(nodes);
      };

      /*!
      * get the cells to which the node serves as corners. Returns the cell indices in a CGNS 
      * compliant winding scheme
      */
      void connectedCells(int32_t node_id,int32_t* cells)override{
        for(int kb=0;kb<8;kb++){cells[kb]=node_connected_cells.at(node_id).at(kb);};
      };
      /*!
      * get the cells to which the node serves as corners. Returns the cell indices in a CGNS 
      * compliant winding scheme
      */
      std::shared_ptr<PanNDE::Array<int32_t>> connectedCells(int32_t node_id)override{
        auto box=iarraymfg->makeManagedArray();
        box->resize(8);
        for(int kb=0;kb<8;kb++){box->at(kb)=node_connected_cells.at(node_id).at(kb);};
        return std::move(box);
      };

      /*!
      * get the global node index for a given local index (for a partitioned mesh)
      */
      int64_t globalNodeId(int32_t node_id)override{return node_global_id.at(node_id);};
      /*!
      * get the partition which owns the node (i.e., updates values at this node)
      */
      int32_t nodeHomePartition(int32_t node_id)override{return node_owner.at(node_id);};

      /*!
      * get the global cell index for a given local index (for a partitioned mesh)
      */
      int64_t globalCellId(int32_t cell_id)override{return hex_global_id.at(cell_id);};
      /*!
      * get the partition which owns the cell (i.e., updates values at this cell)
      */
      int32_t cellHomePartition(int32_t cell_id)override{return hex_owner.at(cell_id);};

      /*!
      * create a copy of another mesh
      * \param other the source mesh for the copy
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
      * create a copy of another mesh
      * \param other the source mesh for the copy
      */
      void copy(std::shared_ptr<PanNDE::Mesh> other)override{copy(other.get());};

      /*!
      * blank constructor
      */
      HexMesh(){};
      /*!
      * mesh constructor from arrays of nodes, cells and partion identifier
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
      void setCell(int cell_id,std::shared_ptr<PanNDE::Array<int32_t>> box,
                   int32_t owner,int32_t g_id){
        for(int kb=0;kb<8;kb++){
          if(node_owner.size()<=box->at(kb)){throw std::out_of_range("invalid node id");};
          hexes.at(cell_id).at(kb)=box->at(kb);
        };
        hex_owner.at(cell_id)=owner;
        hex_global_id.at(cell_id)=g_id;
      };
      void setNode(int node_id,std::shared_ptr<PanNDE::Array<double>> coord,
                   int32_t owner,int32_t g_id){
        X.at(node_id)=coord->at(0);
        Y.at(node_id)=coord->at(1);
        Z.at(node_id)=coord->at(2);
        node_owner.at(node_id)=owner;
        node_global_id.at(node_id)=g_id;
      };
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

      std::shared_ptr<HostData::HostArrayFactory<double>> darraymfg=
              std::make_shared<HostData::HostArrayFactory<double>>(HostData::HostArrayFactory<double>());
      std::shared_ptr<HostData::HostArrayFactory<int32_t>> iarraymfg=
              std::make_shared<HostData::HostArrayFactory<int32_t>>(HostData::HostArrayFactory<int32_t>());

      std::vector<double> X;
      std::vector<double> Y;
      std::vector<double> Z;
      std::vector<int32_t> node_owner;
      std::vector<int64_t> node_global_id;
      std::vector<std::array<int32_t,8>> node_connected_cells;

      std::vector<std::array<int32_t,8>> hexes;
      std::vector<int32_t> hex_owner;
      std::vector<int64_t> hex_global_id;

      int32_t mypartition;


  };

  /*! \class HexMeshFactory HexMesh.hpp "modules/HostData/include/HexMesh.hpp"
  * Implements a factory class to create a HexMesh object
  *
  */
  class HexMeshFactory : public PanNDE::MeshFactory {
    public:
      /*!
      * make a node
      * \param coord[3] nodal cartesian coordinates
      * \param global_id the global identifier of the node
      * \param owner which mesh partition the node belongs to
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
      * make a cell
      * \param nodes nodal indicies which form the corners of the cell
      * \param Nnodes number of nodes which define the cell (e.g. 8 for hexes, 4 for tets)
      * \param global_id the global identifier of the cell
      * \param owner which mesh partition the cell belongs to
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
      * create a shared mesh
      * \param nodes a C-style array (not preferred) of Node structs
      * \param Nnodes the length of the nodes array
      * \param cells a C-style array (not preferred) of cell structs
      * \param Ncells the length of the cells array
      * \param partition the partition identifier of the mesh to be created
      */
      std::shared_ptr<PanNDE::Mesh> makeManagedMesh(
                      PanNDE::Mesh::Node* nodes,int64_t Nnodes,
                      PanNDE::Mesh::Cell* cells,int64_t Ncells,int32_t partition=0)override{
        auto narray=makeArray<PanNDE::Mesh::Node>(nodes,Nnodes);
        auto carray=makeArray<PanNDE::Mesh::Cell>(cells,Ncells);
        return std::make_shared<HostData::HexMesh>(HostData::HexMesh(narray,carray,partition));
      };
      /*! 
      * create a shared mesh
      * \param nodes a PanNDE array of Node structs
      * \param cells a PanNDE array of cell structs
      * \param partition the partition identifier of the mesh to be created
      */
      std::shared_ptr<PanNDE::Mesh> makeManagedMesh(
                      std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Node>> nodes,
                      std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Cell>> cells,int32_t partition=0)override{
        return std::make_shared<HostData::HexMesh>(HostData::HexMesh(nodes,cells,partition));
      };
      /*! 
      * create a mesh on a raw pointer. not recommended, call deleteMesh() if used.
      * \param nodes a C-style array (not preferred) of Node structs
      * \param Nnodes the length of the nodes array
      * \param cells a C-style array (not preferred) of cell structs
      * \param Ncells the length of the cells array
      * \param partition the partition identifier of the mesh to be created
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
      * create a mesh on a raw pointer. not recommended, call deleteMesh() if used.
      * \param nodes a PanNDE array of Node structs
      * \param cells a PanNDE array of cell structs
      * \param partition the partition identifier of the mesh to be created
      */
      PanNDE::Mesh* newMesh(
                      std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Node>> nodes,
                      std::shared_ptr<PanNDE::Array<PanNDE::Mesh::Cell>> cells,int32_t partition=0)override{
        PanNDE::Mesh* mesh=new HostData::HexMesh(nodes,cells,partition);
        return mesh;
      };

      /*!
      * delete meshes created using newMesh()
      * \param mesh mesh to be deleted
      */
      void deleteMesh(PanNDE::Mesh* mesh)override{delete mesh;};

      /*! 
      * create a mesh array
      */
      std::shared_ptr<PanNDE::Array<std::shared_ptr<PanNDE::Mesh>>> makeManagedMeshArray()override{
        return std::move(makeArray<std::shared_ptr<PanNDE::Mesh>>());
      };
    private:
      template<typename T>
      std::shared_ptr<PanNDE::Array<T>> makeArray(T* raw,int64_t N){
        auto array=makeArray<T>();
        array->resize(N);
        for(int64_t k=0;k<N;k++){array->at(k)=raw[k];};
        return std::move(array);
      };
      template<typename T>
      std::shared_ptr<PanNDE::Array<T>> makeArray(){
        auto arraymfg=std::make_shared<HostData::HostArrayFactory<T>>(HostData::HostArrayFactory<T>());
        auto array=arraymfg->makeManagedArray();
        return std::move(array);
      };
  };
};
