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

#include <vector>
#include <stdexcept>
#include <cstdint>

#include "Mesh.hpp"
#include "Array.hpp"

namespace Stubs {
  class PrePartitionedStubHexMesh : public PanNDE::Mesh {
    public:
      PrePartitionedStubHexMesh(int partId,int Nparts){
        (this->partId)=partId;
        (this->Nparts)=Nparts;
        Nn=1;Nc=1;
        for(int kd=0;kd<3;kd++){Nn*=Dims[kd];};
        for(int kd=0;kd<3;kd++){Nc*=(Dims[kd]-1);};
        makeCells();
        makeNodes();
      };

      int32_t nodeCount()override{return Nn;};
      int32_t cellCount()override{return Nc;};
      int32_t partitionId()override{return partId;};

      int64_t globalNodeId(int32_t node_id)override{return nodes.at(node_id).globalId;};
      int32_t nodeHomePartition(int32_t node_id)override{return nodes.at(node_id).owner;};

      int64_t globalCellId(int32_t cell_id)override{return cells.at(cell_id).globalId;};
      int32_t cellHomePartition(int32_t cell_id)override{return cells.at(cell_id).owner;};
      
      void nodeCoordinate(int32_t node_id,double coord[3])override{
        for(int kd=0;kd<3;kd++){coord[kd]=nodes.at(node_id).coord[kd];};
      };
      void cell(int32_t cell_id,int32_t* nodes)override{
        for(int kb=0;kb<8;kb++){nodes[kb]=cells.at(cell_id).box[kb];};
      };
      void connectedCells(int32_t node_id,int32_t* cells)override{
        for(int kb=0;kb<8;kb++){cells[kb]=nodes.at(node_id).box[kb];};
      };

      std::shared_ptr<PanNDE::Array<double>> nodeCoordinate(int32_t node_id)override{throwIfArray();return nullptr;};
      std::shared_ptr<PanNDE::Array<int32_t>> cell(int32_t cell_id)override{throwIfArray();return nullptr;};
      std::shared_ptr<PanNDE::Array<int32_t>> connectedCells(int32_t node_id)override{throwIfArray();return nullptr;};

      void copy(PanNDE::Mesh* other)override{throwIfCopy();};
      void copy(std::shared_ptr<PanNDE::Mesh> other)override{throwIfCopy();};;

    private:
      void throwIfArray(){throw std::runtime_error("No Array Factory Available. This is a stub.");};
      void throwIfCopy(){throw std::runtime_error("This is a stub. No copy method available.");};
      struct stubnode{
        double coord[3];
        int32_t box[8];
        int32_t owner;
        int64_t globalId;
      };
      struct stubcell{
        int32_t box[8];
        int32_t owner;
        int64_t globalId;
      };
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
      double ds[3]={0.1,0.1,0.01};
      int Nn;int Nc;
      int Dims[3]={11,11,5};
      int partId;int Nparts;

      std::vector<stubnode> nodes;
      std::vector<stubcell> cells;
  };
};