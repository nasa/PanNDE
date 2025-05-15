#!/bin/bash

source build_and_run_resources/dome1_env/env.sh

mpirun -n 4 -outfile-pattern=NetMPITestResult_%r.log ./bin/NetMPITests
mpirun -n 4 -outfile-pattern=VTKIOTestResult_%r.log ./bin/VTKIOTests
mpirun -n 4 -outfile-pattern=HostSolverTestResult_%r.log ./bin/HostSolverTests
