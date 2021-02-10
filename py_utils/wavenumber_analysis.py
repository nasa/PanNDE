# Python based vtk data post processing driver with frequency domain post processing. 
#  - included for the V&V work performed in the paper, and as an example  
#    post processing driver
#
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

import numpy as np
import scipy.interpolate as spinterp
import scipy.optimize as spopt
import vtk
import recon

#user supplied data
save_path="data/UDCFRP"
data_path="%s/vol_data"%save_path
file_root="cfrp"
data_to_process="Vz"
t_index_start=50
number_of_sequential_datasets=75
target_frequency=200e3
z_target_slice=0.1e-3
xd_radius=0.002

# note, if files are indexed 0,1,2,..., 
#   - this requires a time_step_multiplier of 1
#   - this requires timestep_seconds to be the time interval between files
# if files are indexed e.g. 0,50,100,150..., 
#   - this requires a time_step_multiplier of 50
#   - this requires a timestep_seconds to be the simulation time step 
#         (e.g., the time interval between files is 50x timestep_seconds)
timestep_seconds=0.1/target_frequency
time_step_multiplier=1

#post processing
rc=recon.reconstructor("%s/%s"%(data_path,file_root))
Result,kx,ky,Acc=rc.generateWavenumberDataAtFrequency(number_of_sequential_datasets,t_index_start,
                                                      target_frequency,timestep_seconds,
                                                      z_target_slice,data_name=data_to_process,
                                                      time_step_multiplier=time_step_multiplier,
                                                      xdrad=xd_radius);


FVzM=np.abs(Result)
FVzM=FVzM/np.max(np.max(FVzM))

FV=spinterp.RectBivariateSpline(kx,ky,FVzM.T)
theta=np.linspace(0,90,91)*np.pi/180.
kr_max=np.zeros(theta.size)

mid=np.array([int(kx.size/2),int(ky.size/2)])
idx_maxx=mid[1]+np.argmax(FVzM[mid[1],mid[0]:])
idx_maxy=mid[0]+np.argmax(FVzM[mid[1]:,mid[0]])
guess=0.5*(kx[idx_maxx]+ky[idx_maxy])

def f(kr,theta):
  return FV.ev(kr*np.cos(theta),kr*np.sin(theta))

for k in range(theta.size):
  kr_max[k]=spopt.fmin(lambda x: -f(x,theta[k]), guess)

#saving
np.savez("%s/wavenumber_results.npz"%save_path,Result=Result,kx=kx,ky=ky,theta=theta,kr=kr_max)

csv_data=np.array([theta,kr_max]).T

np.savetxt("%s/wavenumber_results.csv"%save_path,csv_data,delimiter=",")

