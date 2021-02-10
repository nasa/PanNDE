// App Utilities
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
#include <vector>
#include <string>
#include <getopt.h>
#include <stdexcept>

#include "Mesh.hpp"
#include "MaterialStubs.hpp"

namespace Utilities {
  template<class T>
  inline std::shared_ptr<T> makeShared(){return std::move(std::make_shared<T>(T()));};

  inline std::string ParseParameterFile(int& argc, char**& argv, int rank=0){
    int c=getopt(argc, argv, "p:");
    if(-1==c && 0==rank){
      printf("required -p <parameter_file> argument\n");
      throw std::runtime_error("required -p <parameter_file> argument");
    };
    return std::string(optarg);
  };

  inline std::string ParseWriteFile(int& argc, char**& argv, int rank=0){
    int c=getopt(argc, argv, "o:");
    if(-1==c && 0==rank){
      printf("required -o <output_file> argument (no ext.)\n");
      throw std::runtime_error("required -o <output_file> argument (no ext.)");
    };
    return std::string(optarg);
  };

  inline std::shared_ptr<PanNDE::Array<double>> MakeWriteTimes(
                                std::shared_ptr<PanNDE::ArrayFactory<double>> f64_maker,
                                double dt,double t_sim,int Ntqw=10){
    auto writeTimes=f64_maker->makeManagedArray();
    int Nwrite=int(t_sim/(Ntqw*dt))+1;
    writeTimes->resize(Nwrite);
    for(int ktw=0;ktw<writeTimes->size();ktw++){
      writeTimes->at(ktw)=dt*Ntqw*ktw;
    };
    return std::move(writeTimes);
  };

  inline void PrintBoundingBox(std::shared_ptr<PanNDE::Mesh> mesh){
    double aabb[6]={0.,0.,0.,0.,0.,0.};
    double pt[3];
    mesh->nodeCoordinate(0,pt);
    for(int kp=0;kp<3;kp++){aabb[kp]=pt[kp];aabb[kp+3]=pt[kp];};
    for(int k=1;k<mesh->nodeCount();k++){
      mesh->nodeCoordinate(k,pt);
      for(int kp=0;kp<3;kp++){
        aabb[kp]=std::min(pt[kp],aabb[kp]);
        aabb[kp+3]=std::max(pt[kp],aabb[kp+3]);
      };
    };
    printf("bounding box:\n");
    printf("  %e %e %e\n",aabb[0],aabb[1],aabb[2]);
    printf("  %e %e %e\n",aabb[3],aabb[4],aabb[5]);
  };

  inline void cellCenter(int kc,std::shared_ptr<PanNDE::Mesh> mesh,double pt_c[3]){
    int32_t box[8];
    double pt[3];
    mesh->cell(kc,box);
    pt_c[0]=0.;pt_c[1]=0.;pt_c[2]=0.;
    for(int kb=0;kb<8;kb++){
      mesh->nodeCoordinate(box[kb],pt);
      for(int kd=0;kd<3;kd++){pt_c[kd]+=pt[kd];};
    };
    for(int kd=0;kd<3;kd++){pt_c[kd]=pt_c[kd]/8.;};
  };

  inline std::shared_ptr<Stubs::MaterialStubs> makeAluminum(){
    return makeShared<Stubs::Aluminum>();
  };
  inline std::shared_ptr<Stubs::MaterialStubs> makeIM7_8ply(){
    double cij[6][6]={{174.89e9,  4.09e9,  4.09e9,     0.0,     0.0,     0.0},
                      {  4.09e9, 15.03e9,  5.01e9,     0.0,     0.0,     0.0},
                      {  4.09e9,  5.01e9, 15.03e9,     0.0,     0.0,     0.0},
                      {     0.0,     0.0,     0.0,  5.01e9,     0.0,     0.0},
                      {     0.0,     0.0,     0.0,     0.0,  6.26e9,     0.0},
                      {     0.0,     0.0,     0.0,     0.0,     0.0,  6.26e9}};
    double density=1570.;
    std::vector<double>angles={0.,45.,-45.,90.,90.,-45.,45.,0.};

    return std::make_shared<Stubs::TransverseIsoPlate>(
                    Stubs::TransverseIsoPlate(cij,density,angles));
  };

  inline std::shared_ptr<PanNDE::FieldBundle> BuildDemoDomain(
                              std::shared_ptr<PanNDE::Mesh>& mesh,
                              std::shared_ptr<Stubs::MaterialStubs> matl,
                              std::shared_ptr<PanNDE::FieldFactory> field_maker,
                              double thickness=1.0,int Ncells_per_ply=1,int Nthickness=0){
    PanNDE::ElasticMaterialNames CoefficientNames;
    auto domain=field_maker->makeEmptyManagedFieldBundle();
    domain->mesh()=mesh;
    auto constfield=field_maker->makeManagedField(mesh,PanNDE::Field::CONSTANT);
    domain->emplaceField(CoefficientNames.density,
                         field_maker->makeManagedField(mesh,PanNDE::Field::NODE));
    constfield->at(0)=matl->density();
    domain->field(CoefficientNames.density)->mapFrom(constfield);

    for(int k1=0;k1<6;k1++){
      for(int k2=k1;k2<6;k2++){
        domain->emplaceField(CoefficientNames.CIJ[k1][k2],
                             field_maker->makeManagedField(mesh,PanNDE::Field::CELL));
        auto cij_field=domain->field(CoefficientNames.CIJ[k1][k2]);
        //#pragma omp parallel for
        for(int kc=0;kc<mesh->cellCount();kc++){
          double pt[3];
          Utilities::cellCenter(kc,mesh,pt);
          int ply_idx=int((pt[2]/thickness)*double(Nthickness)/double(Ncells_per_ply));
          cij_field->at(kc)=matl->CIJ(k1,k2,ply_idx);
        };
      };
    };
    return std::move(domain);
  };
};
