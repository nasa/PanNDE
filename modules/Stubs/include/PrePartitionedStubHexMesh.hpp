/*! \headerfile PrePartitionedStubHexMesh.hpp "modules/Stubs/include/PrePartitionedStubHexMesh.hpp"
* "PrePartitionedStubHexMesh.hpp" contains a simple pre-partitioned hexahedral mesh implementation
* for testing parallel algorithms and domain decomposition approaches.
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

#include <vector>
#include <stdexcept>
#include <cstdint>

#include "Mesh.hpp"
#include "Array.hpp"

namespace Stubs {
  /*! \class PrePartitionedStubHexMesh PrePartitionedStubHexMesh.hpp "modules/Stubs/include/PrePartitionedStubHexMesh.hpp"
  *
  * Implements a pre-partitioned structured hexahedral mesh for testing parallel algorithms.
  * This class directly implements the PanNDE::Mesh interface to provide a simple regular grid
  * that is already split into partitions for use in parallel processing tests.
  * The partitioning is done along the x-direction with proper assignment of ownership
  * and global identifiers to support parallel operations.
  *
  */
  class PrePartitionedStubHexMesh : public PanNDE::Mesh {
    public:
      /*!
      * Constructor that creates a pre-partitioned hexahedral mesh.
      * \param partId int The ID of this partition (0-based index)
      * \param Nparts int The total number of partitions
      * 
      * Creates a structured grid with dimensions specified by the Dims array.
      * The mesh is partitioned along the x-direction with proper ownership assignment
      * at partition boundaries.
      */
      PrePartitionedStubHexMesh(int partId,int Nparts){
        (this->partId)=partId;
        (this->Nparts)=Nparts;
        Nn=1;Nc=1;
        for(int kd=0;kd<3;kd++){Nn*=Dims[kd];};
        for(int kd=0;kd<3;kd++){Nc*=(Dims[kd]-1);};
        makeCells();
        makeNodes();
      };

      /*!
      * Gets the number of nodes in this partition.
      * \return int32_t The number of nodes
      */
      int32_t nodeCount()override{return Nn;};
      
      /*!
      * Gets the number of cells in this partition.
      * \return int32_t The number of cells
      */
      int32_t cellCount()override{return Nc;};
      
      /*!
      * Gets the partition ID of this mesh partition.
      * \return int32_t The partition ID
      */
      int32_t partitionId()override{return partId;};

      /*!
      * Gets the global ID of a node with local ID.
      * \param node_id int32_t The local node ID
      * \return int64_t The global node ID
      */
      int64_t globalNodeId(int32_t node_id)override{return nodes.at(node_id).globalId;};
      
      /*!
      * Gets the home partition of a node.
      * \param node_id int32_t The local node ID
      * \return int32_t The home partition ID
      */
      int32_t nodeHomePartition(int32_t node_id)override{return nodes.at(node_id).owner;};

      /*!
      * Gets the global ID of a cell with local ID.
      * \param cell_id int32_t The local cell ID
      * \return int64_t The global cell ID
      */
      int64_t globalCellId(int32_t cell_id)override{return cells.at(cell_id).globalId;};
      
      /*!
      * Gets the home partition of a cell.
      * \param cell_id int32_t The local cell ID
      * \return int32_t The home partition ID
      */
      int32_t cellHomePartition(int32_t cell_id)override{return cells.at(cell_id).owner;};
      
      /*!
      * Gets the coordinates of a node.
      * \param node_id int32_t The node ID
      * \param coord double[3] The output coordinates (x,y,z)
      */
      void nodeCoordinate(int32_t node_id,double coord[3])override{
        for(int kd=0;kd<3;kd++){coord[kd]=nodes.at(node_id).coord[kd];};
      };
      
      /*!
      * Gets the node IDs that make up a cell.
      * \param cell_id int32_t The cell ID
      * \param nodes int32_t* Array to receive the 8 node IDs
      */
      void cell(int32_t cell_id,int32_t* nodes)override{
        for(int kb=0;kb<8;kb++){nodes[kb]=cells.at(cell_id).box[kb];};
      };
      
      /*!
      * Gets the cell IDs connected to a node.
      * \param node_id int32_t The node ID
      * \param cells int32_t* Array to receive the connected cell IDs
      */
      void connectedCells(int32_t node_id,int32_t* cells)override{
        for(int kb=0;kb<8;kb++){cells[kb]=nodes.at(node_id).box[kb];};
      };

      /*!
      * Not implemented - array-based node coordinates.
      * \throws std::runtime_error Always throws as this is not implemented
      */
      std::shared_ptr<PanNDE::Array<double>> nodeCoordinate(int32_t node_id)override{throwIfArray();return nullptr;};
      
      /*!
      * Not implemented - array-based cell nodes.
      * \throws std::runtime_error Always throws as this is not implemented
      */
      std::shared_ptr<PanNDE::Array<int32_t>> cell(int32_t cell_id)override{throwIfArray();return nullptr;};
      
      /*!
      * Not implemented - array-based connected cells.
      * \throws std::runtime_error Always throws as this is not implemented
      */
      std::shared_ptr<PanNDE::Array<int32_t>> connectedCells(int32_t node_id)override{throwIfArray();return nullptr;};

      /*!
      * Not implemented - copy from pointer.
      * \throws std::runtime_error Always throws as this is not implemented
      */
      void copy(PanNDE::Mesh* other)override{throwIfCopy();};
      
      /*!
      * Not implemented - copy from shared pointer.
      * \throws std::runtime_error Always throws as this is not implemented
      */
      void copy(std::shared_ptr<PanNDE::Mesh> other)override{throwIfCopy();};;

    private:
      /*!
      * Throws error for unimplemented array methods.
      * \throws std::runtime_error With message about stub implementation
      */
      void throwIfArray(){throw std::runtime_error("No Array Factory Available. This is a stub.");};
      
      /*!
      * Throws error for unimplemented copy methods.
      * \throws std::runtime_error With message about stub implementation
      */
      void throwIfCopy(){throw std::runtime_error("This is a stub. No copy method available.");};
      
      //! Internal node structure with coordinates, connections, and partition information
      struct stubnode{
        double coord[3];         //!< 3D coordinates
        int32_t box[8];          //!< Connected cells (up to 8)
        int32_t owner;           //!< Owner partition
        int64_t globalId;        //!< Global node ID
      };
      
      //! Internal cell structure with connectivity and partition information
      struct stubcell{
        int32_t box[8];          //!< Node IDs forming the cell
        int32_t owner;           //!< Owner partition
        int64_t globalId;        //!< Global cell ID
      };
      
      /*!
      * Creates all hexahedral cells in the mesh.
      */
      void makeCells(){
        cells.resize(Nc);
        int idx;
        for(int kx=0;kx<(Dims[0]-1);kx++){
          for(int ky=0;ky<(Dims[1]-1);ky++){
            for(int kz=0;kz<(Dims[2]-1);kz++){
              idx=kz+(Dims[2]-1)*(ky+(Dims[1]-1)*kx);

              cells.at(idx).box[0]=(kz)+(Dims[2]*((ky)+Dims[1]*(kx)));
              cells.at(idx).box[1]=(kz)+(Dims[2]*((ky)+Dims[1]*(kx+1)));
              cells.at(idx).box[2]=(kz)+(Dims[2]*((ky+1)+Dims[1]*(kx+1)));
              cells.at(idx).box[3]=(kz)+(Dims[2]*((ky+1)+Dims[1]*(kx)));
              cells.at(idx).box[4]=(kz+1)+(Dims[2]*((ky)+Dims[1]*(kx)));
              cells.at(idx).box[5]=(kz+1)+(Dims[2]*((ky)+Dims[1]*(kx+1)));
              cells.at(idx).box[6]=(kz+1)+(Dims[2]*((ky+1)+Dims[1]*(kx+1)));
              cells.at(idx).box[7]=(kz+1)+(Dims[2]*((ky+1)+Dims[1]*(kx)));

              cells.at(idx).owner=partId;
              cells.at(idx).globalId=kz+Dims[2]*(ky+Dims[1]*(kx+partId*(Dims[0]-2)));
            };
          };
        };
      };
      
      /*!
      * Creates all nodes in the mesh and establishes ownership across partitions.
      */
      void makeNodes(){
        nodes.resize(Nn);
        int idx;
        for(int kx=0;kx<Dims[0];kx++){
          int kx_global=partId*(Dims[0]-1)+kx;
          int owner=partId;
          if(partId>0 && kx==0){owner=partId-1;};
          if(partId<(Nparts-1) && kx==(Dims[0]-1)){owner=partId+1;};
          for(int ky=0;ky<Dims[1];ky++){
            for(int kz=0;kz<Dims[2];kz++){
              idx=kz+Dims[2]*(ky+Dims[1]*kx);
              nodes.at(idx).coord[0]=kx*ds[0];
              nodes.at(idx).coord[1]=ky*ds[1];
              nodes.at(idx).coord[2]=kz*ds[2];
              nodes.at(idx).owner=owner;
              nodes.at(idx).globalId=kz+Dims[2]*(ky+Dims[1]*kx_global);
              for(int kb=0;kb<8;kb++){nodes.at(idx).box[kb]=-1;};
            };
          };
        };
        for(int kc=0;kc<Nc;kc++){
          nodes.at(cells.at(kc).box[0]).box[6]=kc;
          nodes.at(cells.at(kc).box[1]).box[7]=kc;
          nodes.at(cells.at(kc).box[2]).box[4]=kc;
          nodes.at(cells.at(kc).box[3]).box[5]=kc;
          nodes.at(cells.at(kc).box[4]).box[2]=kc;
          nodes.at(cells.at(kc).box[5]).box[3]=kc;
          nodes.at(cells.at(kc).box[6]).box[0]=kc;
          nodes.at(cells.at(kc).box[7]).box[1]=kc;
        };
      };
      
      //! Grid spacing in each dimension
      double ds[3]={0.1,0.1,0.01};
      //! Number of nodes
      int Nn;
      //! Number of cells
      int Nc;
      //! Number of nodes in each dimension [x,y,z]
      int Dims[3]={11,11,5};
      //! ID of this partition
      int partId;
      //! Total number of partitions
      int Nparts;

      //! Collection of all nodes in the partition
      std::vector<stubnode> nodes;
      //! Collection of all cells in the partition
      std::vector<stubcell> cells;
  };
};