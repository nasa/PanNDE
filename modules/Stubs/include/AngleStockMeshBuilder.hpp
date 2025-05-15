/*! \headerfile AngleStockMeshBuilder.hpp "modules/Stubs/include/AngleStockMeshBuilder.hpp"
* "AngleStockMeshBuilder.hpp" contains a utility class for creating L-shaped angle stock meshes
* for testing and demonstration purposes.
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
#include <vector>
#include <unordered_map>
#include <array>
#include <memory>

#include "Mesh.hpp"

namespace Stubs {
  /*! \class AngleStockMeshBuilder AngleStockMeshBuilder.hpp "modules/Stubs/include/AngleStockMeshBuilder.hpp"
  *
  * Implements a mesh generator for creating L-shaped angle stock meshes with
  * hexahedral elements. This utility is primarily used for testing and examples,
  * providing a simple L-profile geometry with configurable dimensions.
  *
  */
  class AngleStockMeshBuilder {
    public:
      /*!
      * Creates a hexahedral mesh representing an L-shaped angle stock.
      * \param meshFactory std::shared_ptr<PanNDE::MeshFactory> Factory for creating mesh elements
      * \param Length double Length of the angle stock in meters (default: 50.8e-3)
      * \param thickness double Thickness of the angle stock in meters (default: 3.175e-3)
      * \param outerHeight double Height of the vertical arm in meters (default: 25.4e-3)
      * \param outerWidth double Width of the horizontal arm in meters (default: 25.4e-3)
      * \param Nthickness int Number of cells through the thickness (default: 16)
      * \return std::shared_ptr<PanNDE::Mesh> The created angle stock mesh
      */
      std::shared_ptr<PanNDE::Mesh> makeAngleMesh(
                          std::shared_ptr<PanNDE::MeshFactory> meshFactory,
                          double Length=50.8e-3,double thickness=3.175e-3,
                          double outerHeight=25.4e-3,double outerWidth=25.4e-3, 
                          int Nthickness=16){
        this->Nthickness=Nthickness;
        this->Length=Length;
        this->thickness=thickness;
        this->outerWidth=outerWidth;
        this->outerHeight=outerHeight;
        this->meshFactory=meshFactory;

        setGridSpacingAndDimensions();
        setAxes();
        makeNodes();
        makeCells();
        return std::move(meshFactory->makeManagedMesh(nodes.data(),nodes.size(),
                                                      cells.data(),cells.size()));
      };

    private:
      /*!
      * Calculates the grid spacing and dimensions based on input parameters.
      */
      void setGridSpacingAndDimensions(){
        ds=thickness/((double)(Nthickness-1));
        N[0]=int(Length/ds+1);
        N[1]=int(outerWidth/ds+1);
        N[2]=int(outerHeight/ds+1);
      };

      /*!
      * Creates the coordinate vectors for the x, y, and z axes.
      */
      void setAxes(){
        for(int kx=0;kx<N[0];kx++){x.push_back(kx*ds);};
        for(int ky=0;ky<N[1];ky++){y.push_back(ky*ds);};
        for(int kz=0;kz<N[2];kz++){z.push_back(kz*ds);};
      };

      /*!
      * Creates all nodes in the mesh for the L-shaped profile.
      * First creates nodes for the horizontal arm, then for the vertical arm.
      */
      void makeNodes(){
        int nidx=0;
        // Horizontal arm nodes
        for(int kx=0;kx<N[0];kx++){
          for(int ky=0;ky<N[1];ky++){
            for(int kz=0;kz<=Nthickness;kz++){
              makeNode(kx,ky,kz,nidx);
            };
          };
        };
        // Vertical arm nodes
        for(int kx=0;kx<N[0];kx++){
          for(int kz=(Nthickness+1);kz<N[2];kz++){
            for(int ky=0;ky<=Nthickness;ky++){
              makeNode(kx,ky,kz,nidx);
            };
          };
        };
      };

      /*!
      * Creates a single node and adds it to the mesh.
      * \param kx int X-index of the node
      * \param ky int Y-index of the node
      * \param kz int Z-index of the node
      * \param nidx int& Reference to the running node counter, incremented after use
      */
      void makeNode(int kx,int ky,int kz,int& nidx){
        double coord[3]={x.at(kx),y.at(ky),z.at(kz)};
        nodes.push_back(meshFactory->makeNode(coord,nidx));
        volmap.emplace(kz+N[2]*(ky+N[1]*kx),nidx);
        nidx++;
      };

      /*!
      * Creates all hexahedral cells in the mesh by connecting nodes.
      * First creates cells for the horizontal arm, then for the vertical arm.
      */
      void makeCells(){
        // Horizontal arm cells
        for(int kx=0;kx<(N[0]-1);kx++){
          for(int ky=0;ky<(N[1]-1);ky++){
            for(int kz=0;kz<Nthickness;kz++){
              makeCell(kx,ky,kz);
            };
          };
        };
        // Vertical arm cells
        for(int kx=0;kx<(N[0]-1);kx++){
          for(int kz=Nthickness;kz<(N[2]-1);kz++){
            for(int ky=0;ky<Nthickness;ky++){
              makeCell(kx,ky,kz);
            };
          };
        };
      };

      /*!
      * Creates a single hexahedral cell and adds it to the mesh.
      * \param kx int X-index of the cell
      * \param ky int Y-index of the cell
      * \param kz int Z-index of the cell
      */
      void makeCell(int kx,int ky,int kz){
        int32_t box[8];
        int64_t gid=cells.size();
        int idx=kz+N[2]*(ky+N[1]*kx);
        box[0]=volmap.at(idx);
        box[1]=volmap.at(idx+N[2]*N[1]);
        box[2]=volmap.at(idx+N[2]*N[1]+N[2]);
        box[3]=volmap.at(idx+N[1]);
        box[4]=volmap.at(idx+1);
        box[5]=volmap.at(idx+1+N[2]*N[1]);
        box[6]=volmap.at(idx+1+N[2]*N[1]+N[2]);
        box[7]=volmap.at(idx+1+N[1]);
        cells.push_back(meshFactory->makeCell(box,8,gid));
      };

      //! Number of cells through the thickness
      int32_t Nthickness=16;
      //! Length of the angle stock in meters
      double Length=50.8e-3;
      //! Thickness of the angle stock in meters
      double thickness=3.175e-3;
      //! Height of the vertical arm in meters
      double outerHeight=25.4e-3;
      //! Width of the horizontal arm in meters
      double outerWidth=25.4e-3;

      //! Cell spacing (uniform in all directions)
      double ds;
      //! Number of nodes in each dimension [x,y,z]
      int N[3];

      //! Vector of x-coordinates
      std::vector<double> x;
      //! Vector of y-coordinates
      std::vector<double> y;
      //! Vector of z-coordinates
      std::vector<double> z;

      //! Collection of all nodes in the mesh
      std::vector<PanNDE::Mesh::Node> nodes;
      //! Mapping from 3D index to node ID
      std::unordered_map<int64_t,int32_t> volmap;
      //! Collection of all cells in the mesh
      std::vector<PanNDE::Mesh::Cell> cells;

      //! Factory for creating mesh elements
      std::shared_ptr<PanNDE::MeshFactory> meshFactory;
  };
};