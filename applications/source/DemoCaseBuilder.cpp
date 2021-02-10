// Demo Case Construction -- use yaml
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

#include "StandardNames.hpp"

#include "Array.hpp"
#include "Mesh.hpp"
#include "Field.hpp"

#include "HostArray.hpp"
#include "HexMesh.hpp"
#include "HostField.hpp"

#include "Gateway.hpp"

#include "VTKGateway.hpp"

#include "AngleStockMeshBuilder.hpp"
#include "MaterialStubs.hpp"

#include "Utilities.hpp"

//////---------------------------

int main(int argc, char** argv){
  auto angle_maker=Utilities::makeShared<Stubs::AngleStockMeshBuilder>();
  auto mesh_maker=Utilities::makeShared<HostData::HexMeshFactory>();
  auto field_maker=Utilities::makeShared<HostData::HostFieldFactory>();
  auto f64_maker=Utilities::makeShared<HostData::HostArrayFactory<double>>();
  auto gateway=Utilities::makeShared<VTKIO::VTKGateway>();
  
  double Length=50.8e-3;
  double thickness=3.175e-3;
  double outerHeight=25.4e-3;
  double outerWidth=25.4e-3;
  int Nthickness=32;
  double xducer_center[3]={Length/4.0,outerWidth/2.0,0.};
  double xducer_radius=1.25*thickness;
  double frequency=200e3;
  int Ncycle=3;
  double phase=0.;

  auto matl=Utilities::makeAluminum();
  
  double CFL=0.85;
  double ds=thickness/(double(Nthickness));
  double Vbulk=sqrt((matl->CIJ(0,0,0))/(matl->density()));//m/s
  double dt=CFL*ds/Vbulk;//sec
  int Ntqwrite=50;
  int Nwrite=200;
  double t_sim=Ntqwrite*dt*Nwrite;
  auto write_times=Utilities::MakeWriteTimes(f64_maker,dt,t_sim,Ntqwrite);

  auto mesh=angle_maker->makeAngleMesh(mesh_maker,Length,thickness,
                                       outerHeight,outerWidth,Nthickness);

  auto domain=Utilities::BuildDemoDomain(mesh,matl,field_maker,thickness,
                                         Nthickness,Nthickness);

  auto metadata=f64_maker->makeManagedDataBundle();
  //timing and transducer arrays need to be loaded to metadata
  PanNDE::TimeDomainMetadataNames TimingNames;
  metadata->emplaceScalar(TimingNames.dt,dt);
  metadata->emplaceArray(TimingNames.write_times,write_times);


  PanNDE::TransducerWindowedSineParameterNames xdNames;
  metadata->emplaceScalar(xdNames.transducerCount(),1);
  metadata->emplaceScalar(xdNames.centerX(0),xducer_center[0]);
  metadata->emplaceScalar(xdNames.centerY(0),xducer_center[1]);
  metadata->emplaceScalar(xdNames.centerZ(0),xducer_center[2]);
  metadata->emplaceScalar(xdNames.frequency(0),frequency);
  metadata->emplaceScalar(xdNames.cycleCount(0),Ncycle);
  metadata->emplaceScalar(xdNames.phase(0),phase);
  metadata->emplaceScalar(xdNames.radius(0),xducer_radius);

  gateway->writeSolution("data/demo_case_",domain,metadata,0);

  return 0;
};