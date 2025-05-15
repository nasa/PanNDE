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
  - VTK >= 9.0 (linked with the same MPI being used for PanNDE)
  - Metis
  - gcc >= 7.3
  - cmake >= 3.13.5
3. Ensure Python3 is installed for some included post processing utilities
  - this is optional, if the utilities are not desired
  - recommended to use a python virtual environement 
  - required packages are in py_utils/requirements.txt
4. cd to `PanNDE/root/directory`
5. Ensure that MPI, VTK, and Metis are in your search path
  - example `export PATH=/path/to/package/bin:$PATH`
6. generate compilation instructions using `cmake -B build -S .`
  - Available flags (default is ON) are:
    - `DGO_TEST:bool=False` turns off building tests
    - `DGO_PROF:bool=False` turns off compiling with profiler hooks
7. compile with `make -C build`

# Running the Test Suites

1. Ensure PanNDE is built on the system
2. Run serial tests:
  - `./bin/HostDataTests > HostDataTestResults.log`
3. Run parallel tests:
  - `mpirun -n 4 -outfile-pattern=NetMPITestResult_%r.log ./bin/NetMPITests`
  - `mpirun -n 4 -outfile-pattern=VTKIOTestResult_%r.log ./bin/VTKIOTests`
  - `mpirun -n 4 -outfile-pattern=HostSolverTestResult_%r.log ./bin/HostSolverTests`


# Python Utilities 

All python utilities can be found in ./py_utils.  Python >= 3.11.11 is required along with the following packages:
- matplotlib>=3.10.0
- numpy>=2.2.1
- scipy>=1.14.1
- vtk>=9.4.1 


# Creating a Demonstration Case

1. Ensure PanNDE is built on the system
2. All examples can bebuilt using the files in ./py_utils/Case_Makers
  -  Composite_Laminate_Transducer_Case_Maker.py
    - run with `TransducerModel`
  -  Small_Composite_Laminate_Transducer_Case_Maker.py
    - run with `TransducerModel`
  -  Large_Plate_wHoles_Transducer_Case_Maker.py
    - run with `TransducerModel`
  -  Small_Plate_wHoles_Transducer_Case_Maker.py
    - run with `TransducerModel`
  -  Ti64_Angle_Bracket_Water_Column_Case_Maker.py
    - run with `WaterColumnModel`
  -  AL_Angle_Bracket_Crack_Transducer_Case_Maker_Parametric.py
    - run with `TransducerModel`
  -  AL_Angle_Bracket_Water_Column_Case_Maker.py
    - run with `WaterColumnModel`

# Running a Case

1. Ensure PanNDE is built on the system
2. Produce a `*.vtu` file that provides all required input data
  - For required data, see Job File Requirements, below
3. Run `mpirun -n $NPROCS ./bin/$MODEL -f $JOB_FILE -o $RESULT_FILE_NAME_ROOT`
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
3. Input signals (one of these is necessary)
  - a. Transducers
    - Number of transducers: `NTransducers`
    - X-location of each transducer: `XD$TRANSDUCER_NUMBER/XCenter`
    - Y-location of each transducer: `XD$TRANSDUCER_NUMBER/YCenter`
    - Z-location of each transducer: `XD$TRANSDUCER_NUMBER/ZCenter`
    - Radius of each transducer: `XD$TRANSDUCER_NUMBER/Radius`
    - Signal Times of each transducer: `XD$TRANSDUCER_NUMBER/SignalTimes`
    - Signal Values of each transducer: `XD$TRANSDUCER_NUMBER/SignalValues`
  - b. Water Columns
    - Number of water columns: `NWaterColumns`
    - X-location of each water column: `WC$WATER_COlUMN_NUMBER/XCenter`
    - Y-location of each water column: `WC$WATER_COlUMN_NUMBER/YCenter`
    - Z-location of each water column: `WC$WATER_COlUMN_NUMBER/ZCenter`
    - Radius of each water column: `WC$WATER_COlUMN_NUMBER/Radius`
    - Signal Times of each water column: `WC$WATER_COlUMN_NUMBER/SignalTimes`
    - Signal Values of each water column in the XX direction: `XD$TRANSDUCER_NUMBER/SignalSignalSigmaXX`
    - Signal Values of each water column in the YY direction: `XD$TRANSDUCER_NUMBER/SignalSignalSigmaYY`
    - Signal Values of each water column in the YY direction: `XD$TRANSDUCER_NUMBER/SignalSignalSigmaZZ`




