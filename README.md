# Welcome to PanNDE!

PanNDE is an attempt at building a modular, developable infrastructure for doing field simulations (initially elastodynamic) for Non-Destructive Evaluation applications. Performance as well as low-level access is a key concern to enable rapid, adaptive modelling for model inversion or for MAPOD like applications. As with all research, work begets work, and development will hopefully continue. For now, if you find use for this, great, and if not, hopefully future development will bring features that are useful, or that the code base provides a jumping-off point for other development by other researchers. For now, good luck, and have fun! Science is a grand adventure!

# Notices:

Copyright 2021 United States Government as represented by the Administrator of the National Aeronautics and Space Administration. No copyright is claimed in the United States under Title 17, U.S. Code. All Other Rights Reserved.
 
Googletest is a product of Google Inc. and is subject to the following:
 
Copyright 2008, Google Inc. All rights reserved.
           
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* Neither the name of Google Inc. nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
# Disclaimers

No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS." 
 
Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.

# Building PanNDE

1. Download latest release
2. Ensure correct dependencies are installed, and are in the search path
  - MPI
  - VTK >= 9.0 (built and linked with the same MPI being used for PanNDE)
  - Metis
  - gcc >= 7.3
  - cmake >= 3.13.5
3. Ensure Python3 is installed for some included post processing utilities
  - this is optional, if the utilities are not desired
4. cd to `PanNDE/root/directory`
5. generate compilation instructions using `cmake -B build -S .`
6. compile with `make -C build`

# Running the Test Suites

1. Ensure PanNDE is built on the system
2. Run serial tests:
  - `./bin/HostDataTests > HostDataTestResults.log`
3. Run parallel tests:
  - `mpirun -n 4 -outfile-pattern=NetMPITestResult_%r.log ./bin/NetMPITests`
  - `mpirun -n 4 -outfile-pattern=VTKIOTestResult_%r.log ./bin/VTKIOTests`
  - `mpirun -n 4 -outfile-pattern=HostSolverTestResult_%r.log ./bin/HostSolverTests`


# Creating a Demonstration Case

1. Ensure PanNDE is built on the system
2. There are three demonstration case builders provided:
  - Aluminum angle stock with single transducer:
    1. Run `./bin/DemoCaseBuilder`
    2. The resulting case file will be `./data/demo_case_0.vtu`
  - Quasi-isotropic 8-ply IM7 CFRP layup with single transducer:
    1. Run `./bin/PlateCaseBuilder -o $JOB_FILENAME_NO_EXTENSION`
    2. The resulting case file will be named with the filename provided
  - Parametric plate with two transducers:
    1. Run `./bin/ParameterizedDemoPlateCase -p $PARAMETER_FILE_YAML`
    2. The resulting case file will be generated based on the parameters contained in the provided parameter file
    3. A sample parameter file that produces the previous quasi-isotropic 8-ply CFRP case is provided in `demo_quasi_iso_cfrp.yaml`. The case can be changed to suit user's use cases, however stability criteria are the user's responsibility

# Running a Case

1. Ensure PanNDE is built on the system
2. Produce a `*.vtu` file that provides all required input data
  - For required data, see Job File Requirements, below
3. Run `mpirun -n $NPROCS ./bin/DemoModelNxd -f $JOB_FILE -o $RESULT_FILE_NAME_ROOT`
  - results will be written to a `*.pvtu`/`*.vtu` set with appropriate file numbering schema for reconstruction in scientific visualization software

# Job File Requirements

1. Elastodynamic parameter fields:
  - cell centered stiffness coefficients named:
    - `C11 C12 C13 C14 C15 C16 C22 C23 C24 C25 C26 C33 C34 C35 C36 C44 C45 C46 C55 C56 C66`
  - cell centered density values named: 
    - `density`
2. Time step and file write frequency information:
  - simulation time step named:
    - `dt`
  - simulation write times array named:
    - `write_times`
3. Hanning windowed Ncycle sine transducer excitation information:
  - Number of transducers named:
    - `NTransducers`
  - Transducer parameters:
    - `XD<transducer index>/XCenter`
    - `XD<transducer index>/YCenter`
    - `XD<transducer index>/ZCenter`
    - `XD<transducer index>/Frequency`
    - `XD<transducer index>/NCycle`
    - `XD<transducer index>/Phase`
    - `XD<transducer index>/Radius`

# Future Work/Features

1. Generic transducers, or more transducer type options
2. PZT model
3. Thermal conduction model
4. GPU data allocator and solver

