// Material stub
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

#include <cmath>
#include <map>
#include <vector>
#include <stdexcept>

namespace Stubs {
  class MaterialStubs {
    public:
      virtual int Nplys(){return 1;};
      virtual double density() =0;
      virtual double CIJ(int ki, int kj, int ply_idx=0) =0;
      virtual double LongitudinalVelocity() =0;
      virtual double TransverseVelocity() =0;
  };
  class TransverseIsoPlate : public Stubs::MaterialStubs {
    public:
      TransverseIsoPlate(double cij_along_x[6][6],double density,
        std::vector<double> layup={0.,45.,135.,90.,90.,135.,45.,0.}){
        for(int ki=0;ki<6;ki++){for(int kj=0;kj<6;kj++){cij[ki][kj]=cij_along_x[ki][kj];};};
        ply_angles=layup;
        ply_stiffnesses.resize(ply_angles.size());
        toIJKL();
        for(int k=0;k<ply_angles.size();k++){
          rotate(k);
          copyCij(cij_r,ply_stiffnesses.at(k));
        };
        rotate(0);
      };
      //IM7 default
      TransverseIsoPlate(std::vector<double> layup={0.,45.,135.,90.,90.,135.,45.,0.}){
        ply_angles=layup;
        toIJKL();
        rotate(0);
      };
      int Nplys()override{return ply_angles.size();};
      double density()override{return rho;};
      double CIJ(int ki,int kj,int ply_idx=0)override{
        if(ki>5 || kj>5 || ki<0 || kj<0){throw std::out_of_range("invalid entry");};
        //if(ply_idx!=current_ply){rotate(ply_idx);};
        //return cij_r[ki][kj];
        return ply_stiffnesses.at(ply_idx).cij[ki][kj];
      };
      double LongitudinalVelocity()override{return sqrt(cij[0][0]/density());};
      double TransverseVelocity()override{return sqrt(cij[1][1]/density());};
    private:
      struct Stiff{
        double cij[6][6];
      };
      void copyCij(double src[6][6],Stiff& dest){
        for(int ki=0;ki<6;ki++){for(int kj=0;kj<6;kj++){dest.cij[ki][kj]=src[ki][kj];};};
      };

      void rotate(int ply){
        rotationMatrix(ply_angles.at(ply)*M_PI/180.);
        rotateCij();
      };
      void rotationMatrix(double alpha){
        R[0][0]=cos(alpha);
        R[1][1]=cos(alpha);
        R[0][1]=-sin(alpha);
        R[1][0]=sin(alpha);
      };
      void rotateCij(){
        for(int k1=0;k1<3;k1++){
          for(int k2=0;k2<3;k2++){
            for(int k3=0;k3<3;k3++){
              for(int k4=0;k4<3;k4++){
                cijkl_r[k1][k2][k3][k4]=0.;
        for(int ki=0;ki<3;ki++){
          for(int kj=0;kj<3;kj++){
            for(int kk=0;kk<3;kk++){
              for(int kl=0;kl<3;kl++){
                cijkl_r[k1][k2][k3][k4]+=R[k1][ki]*R[k2][kj]*R[k3][kk]*R[k4][kl]*cijkl[ki][kj][kk][kl];
              };
            };
          };
        };
              };
            };
          };
        };
        toIJ();
      };
      void toIJKL(){
        for(int k1=0;k1<3;k1++){
          for(int k2=0;k2<3;k2++){
            for(int k3=0;k3<3;k3++){
              for(int k4=0;k4<3;k4++){
                cijkl[k1][k2][k3][k4]=cij[ijtable[k1][k2]][ijtable[k3][k4]];
              };
            };
          };
        };
      };
      void toIJ(){
        for(int k1=0;k1<3;k1++){
          for(int k2=0;k2<3;k2++){
            for(int k3=0;k3<3;k3++){
              for(int k4=0;k4<3;k4++){
                cij_r[ijtable[k1][k2]][ijtable[k3][k4]]=cijkl_r[k1][k2][k3][k4];
              };
            };
          };
        };
      };
      int ijtable[3][3]={{0,5,4},
                         {5,1,3},
                         {4,3,2}};

      double R[3][3]={{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
      std::vector<double> ply_angles;
      double cij[6][6]={
                          {174.89e9,  4.09e9,  4.09e9,     0.0,     0.0,     0.0},
                          {  4.09e9, 15.03e9,  5.01e9,     0.0,     0.0,     0.0},
                          {  4.09e9,  5.01e9, 15.03e9,     0.0,     0.0,     0.0},
                          {     0.0,     0.0,     0.0,  5.01e9,     0.0,     0.0},
                          {     0.0,     0.0,     0.0,     0.0,  6.26e9,     0.0},
                          {     0.0,     0.0,     0.0,     0.0,     0.0,  6.26e9}
                       };
      double cijkl[3][3][3][3];
      double cij_r[6][6];
      double cijkl_r[3][3][3][3];
      int current_ply=0;
      double rho=1570.;
      std::vector<Stiff> ply_stiffnesses;
  };
  class Aluminum : public Stubs::MaterialStubs {
    public:
      Aluminum(){
        double c_ij[6][6]={
                            {lam+2*mu,     lam,     lam,0.,0.,0.},
                            {lam,     lam+2*mu,     lam,0.,0.,0.},
                            {lam,          lam,lam+2*mu,0.,0.,0.},
                            {0.,            0.,      0.,mu,0.,0.},
                            {0.,            0.,      0.,0.,mu,0.},
                            {0.,            0.,      0.,0.,0.,mu}
                          };
        for(int ki=0;ki<6;ki++){for(int kj=0;kj<6;kj++){cij[ki][kj]=c_ij[ki][kj];};};
      };
      double density()override{return rho;};
      double CIJ(int ki,int kj,int ply_idx=0)override{
        if(ki>5 || kj>5 || ki<0 || kj<0){throw std::out_of_range("invalid entry");};
        return cij[ki][kj];
      };
      double LongitudinalVelocity()override{return sqrt(cij[0][0]/density());};
      double TransverseVelocity()override{return sqrt(cij[1][1]/density());};
    private:
      double rho=2780.0;
      double lam=51.749939080e9;
      double mu=26.664117020e9;
      double cij[6][6];
  };
};