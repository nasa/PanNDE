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

#include <memory>
#include <cstdint>
#include <vector>
#include <stdexcept>

#include "Field.hpp"
#include "Mesh.hpp"

namespace HostSolver {
  class RightHex /*: public PanNDE::Element*/ {
    public:
      RightHex(std::shared_ptr<PanNDE::Mesh> mesh){
        this->mesh=mesh;
        dx_c.resize(mesh->cellCount());
        dy_c.resize(mesh->cellCount());
        dz_c.resize(mesh->cellCount());
        dx_n.resize(mesh->nodeCount());
        dy_n.resize(mesh->nodeCount());
        dz_n.resize(mesh->nodeCount());
        nodePrep();
        cellPrep();
      }
      double DDx_cc(int idx,std::shared_ptr<PanNDE::Field> field)/*override*/{
        int32_t box[8];
        mesh->cell(idx,box);
        int32_t fwd[4]={box[1],box[2],box[5],box[6]};
        int32_t bwd[4]={box[0],box[3],box[4],box[7]};
        double value=ccDifferencer(fwd,bwd,field);
        return (value*0.25/dx_c.at(idx));
      };
      double DDy_cc(int idx,std::shared_ptr<PanNDE::Field> field)/*override*/{
        int32_t box[8];
        mesh->cell(idx,box);
        int32_t fwd[4]={box[2],box[3],box[6],box[7]};
        int32_t bwd[4]={box[0],box[1],box[4],box[5]};
        double value=ccDifferencer(fwd,bwd,field);
        return (value*0.25/dy_c.at(idx));
      };
      double DDz_cc(int idx,std::shared_ptr<PanNDE::Field> field)/*override*/{
        int32_t box[8];
        mesh->cell(idx,box);
        int32_t fwd[4]={box[4],box[5],box[6],box[7]};
        int32_t bwd[4]={box[0],box[1],box[2],box[3]};
        double value=ccDifferencer(fwd,bwd,field);
        return (value*0.25/dz_c.at(idx));
      };
      double DDx_nc(int idx,std::shared_ptr<PanNDE::Field> field)/*override*/{
        int32_t box[8];
        mesh->connectedCells(idx,box);
        int32_t fwd[4]={box[1],box[2],box[5],box[6]};
        int32_t bwd[4]={box[0],box[3],box[4],box[7]};
        double value=ncDifferencer(fwd,bwd,field);
        return (value*0.25/dx_n.at(idx));
      };
      double DDy_nc(int idx,std::shared_ptr<PanNDE::Field> field)/*override*/{
        int32_t box[8];
        mesh->connectedCells(idx,box);
        int32_t fwd[4]={box[2],box[3],box[6],box[7]};
        int32_t bwd[4]={box[0],box[1],box[4],box[5]};
        double value=ncDifferencer(fwd,bwd,field);
        return (value*0.25/dy_n.at(idx));
      };
      double DDz_nc(int idx,std::shared_ptr<PanNDE::Field> field)/*override*/{
        int32_t box[8];
        mesh->connectedCells(idx,box);
        int32_t fwd[4]={box[4],box[5],box[6],box[7]};
        int32_t bwd[4]={box[0],box[1],box[2],box[3]};
        double value=ncDifferencer(fwd,bwd,field);
        return (value*0.25/(dz_n.at(idx)));
      };

      double getNodeDx(int nidx){return dx_n.at(nidx);};
      double getNodeDy(int nidx){return dy_n.at(nidx);};
      double getNodeDz(int nidx){return dz_n.at(nidx);};
    
    private:
      double ncDifferencer(int fwd[4],int bwd[4],std::shared_ptr<PanNDE::Field> field){
        double value=0.0;
        for(int k=0;k<4;k++){
          value+=(((-1==fwd[k])?0.:field->atCell(fwd[k]))-((-1==bwd[k])?0.:field->atCell(bwd[k])));
        };
        return value;
      };
      double ccDifferencer(int fwd[4],int bwd[4],std::shared_ptr<PanNDE::Field> field){
        double value=0.0;
        for(int k=0;k<4;k++){
          value+=(field->atNode(fwd[k])-field->atNode(bwd[k]));
        };
        return value;
      };
      void cellPrep(){
        int32_t box[8];
        double pt_lrb[3];double pt_ulf[3];
        for(int kc=0;kc<mesh->cellCount();kc++){
          mesh->cell(kc,box);
          mesh->nodeCoordinate(box[0],pt_lrb);
          mesh->nodeCoordinate(box[6],pt_ulf);
          dx_c.at(kc)=pt_ulf[0]-pt_lrb[0];
          dy_c.at(kc)=pt_ulf[1]-pt_lrb[1];
          dz_c.at(kc)=pt_ulf[2]-pt_lrb[2];
        };
      };
      void nodePrep(){
        int32_t box[8];
        double pt_node[3];
        double pt_lrb[3];double pt_ulf[3];
        double pt[3];
        int32_t cbox[8];
        for(int kn=0;kn<mesh->nodeCount();kn++){
          mesh->nodeCoordinate(kn,pt_node);
          mesh->nodeCoordinate(kn,pt_lrb);
          mesh->nodeCoordinate(kn,pt_ulf);
          mesh->connectedCells(kn,box);
          for(int kb=0;kb<8;kb++){
            if(-1!=box[kb]){
              mesh->cell(box[kb],cbox);
              for(int kbb=0;kbb<8;kbb++){
                mesh->nodeCoordinate(cbox[kbb],pt);
                for(int kd=0;kd<3;kd++){
                  pt_lrb[kd]=std::min(pt_lrb[kd],pt[kd]);
                  pt_ulf[kd]=std::max(pt_ulf[kd],pt[kd]);
                };
              };
            };
          };
          for(int kd=0;kd<3;kd++){
            if(pt_lrb[kd]==pt_node[kd]){
              pt_lrb[kd]=2.*pt_node[kd]-pt_ulf[kd];
              //if(0>=pt_lrb[kd]){throw std::logic_error("LRB");};
            };
            if(pt_ulf[kd]==pt_node[kd]){
              pt_ulf[kd]=2.*pt_node[kd]-pt_lrb[kd];
              //if(0>=pt_ulf[kd]){throw std::logic_error("ULF");};
            };
          };
          dx_n.at(kn)=0.5*(pt_ulf[0]-pt_lrb[0]);if(0>=dx_n.at(kn)){throw std::logic_error("DX");};
          dy_n.at(kn)=0.5*(pt_ulf[1]-pt_lrb[1]);if(0>=dy_n.at(kn)){throw std::logic_error("DY");};
          dz_n.at(kn)=0.5*(pt_ulf[2]-pt_lrb[2]);if(0>=dz_n.at(kn)){throw std::logic_error("DZ");};
        };
      };

      std::vector<double> dx_c;
      std::vector<double> dy_c;
      std::vector<double> dz_c;

      std::vector<double> dx_n;
      std::vector<double> dy_n;
      std::vector<double> dz_n;
      
      std::shared_ptr<PanNDE::Mesh> mesh=nullptr;
  };
};