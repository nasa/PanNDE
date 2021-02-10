// Parameterized Demo Case Construction
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

#include <memory>
#include <cmath>
#include <map>
#include <getopt.h>

#include "StandardNames.hpp"

#include "Array.hpp"
#include "Mesh.hpp"
#include "Field.hpp"

#include "HostArray.hpp"
#include "HexMesh.hpp"
#include "HostField.hpp"

#include "Gateway.hpp"

#include "VTKGateway.hpp"

#include "PlateBuilder.hpp"

#include "MaterialStubs.hpp"

#include "simple_parameter_reader.hpp"
#include "Utilities.hpp"

//////---------------------------

int main(int argc, char** argv){
  auto param_file=Utilities::ParseParameterFile(argc,argv);
  auto tree=ParamUtils::BuildTreeFromFile(param_file);

  auto plate_maker=Utilities::makeShared<Stubs::PlateBuilder>();
  auto mesh_maker=Utilities::makeShared<HostData::HexMeshFactory>();
  auto field_maker=Utilities::makeShared<HostData::HostFieldFactory>();
  auto f64_maker=Utilities::makeShared<HostData::HostArrayFactory<double>>();
  auto gateway=Utilities::makeShared<VTKIO::VTKGateway>();

  auto root_node=tree.at(0);
  auto filename_node=root_node->FindChildByName("case filename");

  auto plate_node=root_node->FindChildByName("Square Plate");
  auto length_node=plate_node->FindChildByName("Side Length");
  auto dxy_node=plate_node->FindChildByName("Side spacing");
  auto thickness_node=plate_node->FindChildByName("Ply Thickness");
  auto ply_cell_cnt_node=plate_node->FindChildByName("Cells per ply");
  auto ply_angles_node=plate_node->FindChildByName("Ply Angles");

  auto xd_node=root_node->FindChildByName("Centered Transducer");
  auto xd_radius_node=xd_node->FindChildByName("radius");
  auto xd_freq_node=xd_node->FindChildByName("center frequency");
  auto xd_cycle_node=xd_node->FindChildByName("cycle count");
  auto xd_phase_node=xd_node->FindChildByName("phase");

  auto timing_node=root_node->FindChildByName("Timing");
  auto dt_node=timing_node->FindChildByName("dt");
  auto tqwrite_node=timing_node->FindChildByName("Time steps per write");
  auto sim_time_node=timing_node->FindChildByName("simulation time");

  auto matl_node=root_node->FindChildByName("Material");
  auto density_node=matl_node->FindChildByName("Density");
  auto c11_node=matl_node->FindChildByName("C11");
  auto c22_node=matl_node->FindChildByName("C22");
  auto c12_node=matl_node->FindChildByName("C12");
  auto c23_node=matl_node->FindChildByName("C23");
  auto c44_node=matl_node->FindChildByName("C44");
  auto c55_node=matl_node->FindChildByName("C55");

  std::string fname_base=filename_node->getContents().content;

  double Length=(length_node->GetScalarD())*1e-3;
  double ds=(dxy_node->GetScalarD())*1e-3;
  double ply_thickness=(thickness_node->GetScalarD())*1e-3;
  int cells_q_ply=ply_cell_cnt_node->GetScalarI();
  std::vector<double> angles=ply_angles_node->getContents().fdata;
  double Width=Length;
  int Nlength=int(ceil(Length/ds));
  int Nthickness=angles.size()*cells_q_ply;
  double thickness=ply_thickness*angles.size();
  //double dz=ply_thickness/double(cells_q_ply);

  double xducer_radius=(xd_radius_node->GetScalarD())*1e-3;
  double xducer_center[3]={Length/2.0,Width/2.0,0.};
  double freq=(xd_freq_node->GetScalarD())*1e3;
  int Ncycle=xd_cycle_node->GetScalarI();
  double phase=xd_phase_node->GetScalarD();

  double dt=(dt_node->GetScalarD())*1e-6;
  int N_tqw=tqwrite_node->GetScalarI();
  double t_sim=(sim_time_node->GetScalarD())*1e-6;
  auto write_times=Utilities::MakeWriteTimes(f64_maker,dt,t_sim,N_tqw);

  double density=density_node->GetScalarD();
  double cij[6][6]={{0.,0.,0.,0.,0.,0.},
                    {0.,0.,0.,0.,0.,0.},
                    {0.,0.,0.,0.,0.,0.},
                    {0.,0.,0.,0.,0.,0.},
                    {0.,0.,0.,0.,0.,0.},
                    {0.,0.,0.,0.,0.,0.}};
  cij[0][0]=(c11_node->GetScalarD())*1e9;
  cij[1][1]=(c22_node->GetScalarD())*1e9;
  cij[0][1]=(c12_node->GetScalarD())*1e9;
  cij[1][2]=(c23_node->GetScalarD())*1e9;
  cij[3][3]=(c44_node->GetScalarD())*1e9;
  cij[4][4]=(c55_node->GetScalarD())*1e9;
  cij[2][2]=cij[1][1];
  cij[1][0]=cij[0][1];cij[2][0]=cij[0][1];cij[0][2]=cij[0][1];
  cij[2][1]=cij[1][2];
  cij[5][5]=cij[4][4];

  printf("density: %e\n",density);
  for(int k1=0;k1<6;k1++){
    for(int k2=0;k2<6;k2++){
      printf("  %e",cij[k1][k2]);
    };
    printf("\n");
  };
  for(auto it=angles.begin();it!=angles.end();it++){
    printf("%e  ",*it);
  };
  printf("\n");

  auto matl=std::make_shared<Stubs::TransverseIsoPlate>(
                    Stubs::TransverseIsoPlate(cij,density,angles));

  auto mesh=plate_maker->makeMesh(mesh_maker,Length,Width,thickness,
                                  Nlength,Nthickness);

  Utilities::PrintBoundingBox(mesh);


  auto domain=Utilities::BuildDemoDomain(mesh,matl,field_maker,thickness,
                                         cells_q_ply,Nthickness);

  auto metadata=f64_maker->makeManagedDataBundle();
  //timing and transducer arrays need to be loaded to metadata
  PanNDE::TimeDomainMetadataNames TimingNames;
  metadata->emplaceScalar(TimingNames.dt,dt);
  metadata->emplaceArray(TimingNames.write_times,write_times);

  //updated xd naming scheme
  PanNDE::TransducerWindowedSineParameterNames xdNames;
  metadata->emplaceScalar(xdNames.transducerCount(),2);
  for(int k=0;k<2;k++){
    metadata->emplaceScalar(xdNames.centerX(k),Length/2.);
    metadata->emplaceScalar(xdNames.centerY(k),Width/2.);
    metadata->emplaceScalar(xdNames.frequency(k),freq);
    metadata->emplaceScalar(xdNames.cycleCount(k),Ncycle);
    metadata->emplaceScalar(xdNames.phase(k),phase);
    metadata->emplaceScalar(xdNames.radius(k),xducer_radius);
  };
  metadata->emplaceScalar(xdNames.centerZ(0),0.);
  metadata->emplaceScalar(xdNames.centerZ(1),thickness);

  //gateway->writeSolution(fname_base,domain,metadata,0);
  gateway->writeSolution(fname_base,domain,metadata);
//*/
  return 0;
};
