/*! \headerfile RightHex.hpp "modules/HostSolver/internal/RightHex.hpp"
* "RightHex.hpp" contains a utility class for computing spatial derivatives on
* hexahedral mesh elements. It provides methods for calculating gradients and
* divergence operators using finite difference schemes on structured hexahedral grids.
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
#include <cstdint>
#include <vector>
#include <stdexcept>

#include "Field.hpp"
#include "Mesh.hpp"

namespace HostSolver {
  /*! \class RightHex RightHex.hpp "modules/HostSolver/internal/RightHex.hpp"
  *
  * Implements finite difference schemes for computing spatial derivatives on hexahedral meshes.
  * This class provides methods to calculate derivatives in x, y, and z directions using 
  * staggered grid approaches for both cell-centered and node-centered fields. It automatically
  * computes mesh spacing information needed for accurate finite differencing.
  *
  */
  class RightHex /*: public PanNDE::Element*/ {
    public:
      /*!
      * Constructs a RightHex solver for the specified hexahedral mesh.
      * Computes and stores all necessary geometric information for finite differencing.
      * 
      * \param mesh std::shared_ptr<PanNDE::Mesh> The mesh to operate on
      * \throw std::logic_error If the mesh has invalid spacing or degenerate cells
      */
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
      
      /*!
      * Computes the x-derivative of a field at a cell center using node values.
      * Uses the average of differences between east and west node pairs.
      * 
      * \param idx int Index of the cell
      * \param field std::shared_ptr<PanNDE::Field> The field to differentiate
      * \return double The derivative value at the cell center
      */
      double DDx_cc(int idx,std::shared_ptr<PanNDE::Field> field)/*override*/{
        int32_t box[8];
        mesh->cell(idx,box);
        int32_t fwd[4]={box[1],box[2],box[5],box[6]};
        int32_t bwd[4]={box[0],box[3],box[4],box[7]};
        double value=ccDifferencer(fwd,bwd,field);
        return (value*0.25/dx_c.at(idx));
      };
      
      /*!
      * Computes the y-derivative of a field at a cell center using node values.
      * Uses the average of differences between north and south node pairs.
      * 
      * \param idx int Index of the cell
      * \param field std::shared_ptr<PanNDE::Field> The field to differentiate
      * \return double The derivative value at the cell center
      */
      double DDy_cc(int idx,std::shared_ptr<PanNDE::Field> field)/*override*/{
        int32_t box[8];
        mesh->cell(idx,box);
        int32_t fwd[4]={box[2],box[3],box[6],box[7]};
        int32_t bwd[4]={box[0],box[1],box[4],box[5]};
        double value=ccDifferencer(fwd,bwd,field);
        return (value*0.25/dy_c.at(idx));
      };
      
      /*!
      * Computes the z-derivative of a field at a cell center using node values.
      * Uses the average of differences between top and bottom node pairs.
      * 
      * \param idx int Index of the cell
      * \param field std::shared_ptr<PanNDE::Field> The field to differentiate
      * \return double The derivative value at the cell center
      */
      double DDz_cc(int idx,std::shared_ptr<PanNDE::Field> field)/*override*/{
        int32_t box[8];
        mesh->cell(idx,box);
        int32_t fwd[4]={box[4],box[5],box[6],box[7]};
        int32_t bwd[4]={box[0],box[1],box[2],box[3]};
        double value=ccDifferencer(fwd,bwd,field);
        return (value*0.25/dz_c.at(idx));
      };
      
      /*!
      * Computes the x-derivative of a field at a node using connected cell values.
      * Uses the average of differences between east and west cell pairs.
      * 
      * \param idx int Index of the node
      * \param field std::shared_ptr<PanNDE::Field> The field to differentiate
      * \return double The derivative value at the node
      */
      double DDx_nc(int idx,std::shared_ptr<PanNDE::Field> field)/*override*/{
        int32_t box[8];
        mesh->connectedCells(idx,box);
        int32_t fwd[4]={box[1],box[2],box[5],box[6]};
        int32_t bwd[4]={box[0],box[3],box[4],box[7]};
        double value=ncDifferencer(fwd,bwd,field);
        return (value*0.25/dx_n.at(idx));
      };
      
      /*!
      * Computes the y-derivative of a field at a node using connected cell values.
      * Uses the average of differences between north and south cell pairs.
      * 
      * \param idx int Index of the node
      * \param field std::shared_ptr<PanNDE::Field> The field to differentiate
      * \return double The derivative value at the node
      */
      double DDy_nc(int idx,std::shared_ptr<PanNDE::Field> field)/*override*/{
        int32_t box[8];
        mesh->connectedCells(idx,box);
        int32_t fwd[4]={box[2],box[3],box[6],box[7]};
        int32_t bwd[4]={box[0],box[1],box[4],box[5]};
        double value=ncDifferencer(fwd,bwd,field);
        return (value*0.25/dy_n.at(idx));
      };
      
      /*!
      * Computes the z-derivative of a field at a node using connected cell values.
      * Uses the average of differences between top and bottom cell pairs.
      * 
      * \param idx int Index of the node
      * \param field std::shared_ptr<PanNDE::Field> The field to differentiate
      * \return double The derivative value at the node
      */
      double DDz_nc(int idx,std::shared_ptr<PanNDE::Field> field)/*override*/{
        int32_t box[8];
        mesh->connectedCells(idx,box);
        int32_t fwd[4]={box[4],box[5],box[6],box[7]};
        int32_t bwd[4]={box[0],box[1],box[2],box[3]};
        double value=ncDifferencer(fwd,bwd,field);
        return (value*0.25/(dz_n.at(idx)));
      };

      /*!
      * Gets the x-direction spacing at a specified node.
      * \param nidx int Node index
      * \return double The effective x-spacing for this node
      */
      double getNodeDx(int nidx){return dx_n.at(nidx);};
      
      /*!
      * Gets the y-direction spacing at a specified node.
      * \param nidx int Node index
      * \return double The effective y-spacing for this node
      */
      double getNodeDy(int nidx){return dy_n.at(nidx);};
      
      /*!
      * Gets the z-direction spacing at a specified node.
      * \param nidx int Node index
      * \return double The effective z-spacing for this node
      */
      double getNodeDz(int nidx){return dz_n.at(nidx);};
    
    private:
      /*!
      * Computes differences for node-to-cell derivatives.
      * Handles boundary cases where cells might not exist (-1 indices).
      * 
      * \param fwd int[4] Forward cell indices
      * \param bwd int[4] Backward cell indices
      * \param field std::shared_ptr<PanNDE::Field> The field to differentiate
      * \return double Sum of differences between forward and backward cells
      */
      double ncDifferencer(int fwd[4],int bwd[4],std::shared_ptr<PanNDE::Field> field){
        double value=0.0;
        for(int k=0;k<4;k++){
          value+=(((-1==fwd[k])?0.:field->atCell(fwd[k]))-((-1==bwd[k])?0.:field->atCell(bwd[k])));
        };
        return value;
      };
      
      /*!
      * Computes differences for cell-to-cell derivatives using node values.
      * 
      * \param fwd int[4] Forward node indices
      * \param bwd int[4] Backward node indices
      * \param field std::shared_ptr<PanNDE::Field> The field to differentiate
      * \return double Sum of differences between forward and backward nodes
      */
      double ccDifferencer(int fwd[4],int bwd[4],std::shared_ptr<PanNDE::Field> field){
        double value=0.0;
        for(int k=0;k<4;k++){
          value+=(field->atNode(fwd[k])-field->atNode(bwd[k]));
        };
        return value;
      };
      
      /*!
      * Computes the cell spacing in each direction.
      * For each cell, finds the distance between opposite corners.
      */
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
      
      /*!
      * Computes effective nodal spacing in each direction.
      * For boundary nodes, uses extrapolated spacing to maintain consistency.
      * 
      * \throw std::logic_error If any computed spacing is not positive
      */
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

      //! Cell x-direction spacing
      std::vector<double> dx_c;
      
      //! Cell y-direction spacing
      std::vector<double> dy_c;
      
      //! Cell z-direction spacing
      std::vector<double> dz_c;

      //! Node x-direction effective spacing
      std::vector<double> dx_n;
      
      //! Node y-direction effective spacing
      std::vector<double> dy_n;
      
      //! Node z-direction effective spacing
      std::vector<double> dz_n;
      
      //! Reference to the mesh being operated on
      std::shared_ptr<PanNDE::Mesh> mesh=nullptr;
  };
};