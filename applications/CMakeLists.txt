# CMake file for setting up build of PanNDE application binaries
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

cmake_minimum_required (VERSION 3.13.5)

project(PanNDE_binaries C CXX)#CUDA

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib) #for static lib
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib) # for shared lib

include_directories("${PROJECT_SOURCE_DIR}/include")

add_executable(TransducerModel
    "${PROJECT_SOURCE_DIR}/source/CaseRunners/TransducerModel.cpp"
    )
target_link_libraries(TransducerModel PanNDE_interfaces Controller_interfaces)
target_link_libraries(TransducerModel HostData HostSolver NetMPI VTKIO)
target_link_libraries(TransducerModel Stubs)

add_executable(WaterColumnModel
    "${PROJECT_SOURCE_DIR}/source/CaseRunners/WaterColumnModel.cpp"
    )
target_link_libraries(WaterColumnModel PanNDE_interfaces Controller_interfaces)
target_link_libraries(WaterColumnModel HostData HostSolver NetMPI VTKIO)
target_link_libraries(WaterColumnModel Stubs)

if(GO_LEGACY)
    add_executable(DemoCaseBuilder
        "${PROJECT_SOURCE_DIR}/source/Legacy/CaseMakers/DemoCaseBuilder.cpp"
        )
    target_link_libraries(DemoCaseBuilder PanNDE_interfaces Controller_interfaces)
    target_link_libraries(DemoCaseBuilder Stubs HostData VTKIO)

    add_executable(TTTDemoCaseBuilder
        "${PROJECT_SOURCE_DIR}/source/Legacy/CaseMakers/TTTDemoCaseBuilder.cpp"
        )
    target_link_libraries(TTTDemoCaseBuilder PanNDE_interfaces Controller_interfaces)
    target_link_libraries(TTTDemoCaseBuilder Stubs HostData VTKIO)

    add_executable(PlateCaseBuilder
        "${PROJECT_SOURCE_DIR}/source/Legacy/CaseMakers/PlateCaseBuilder.cpp"
        )
    target_link_libraries(PlateCaseBuilder PanNDE_interfaces Controller_interfaces)
    target_link_libraries(PlateCaseBuilder Stubs HostData VTKIO)

    add_executable(ParameterizedDemoPlateCase
        "${PROJECT_SOURCE_DIR}/source/Legacy/CaseMakers/ParameterizedDemoPlateCase.cpp"
        )
    target_link_libraries(ParameterizedDemoPlateCase PanNDE_interfaces Controller_interfaces)
    target_link_libraries(ParameterizedDemoPlateCase Stubs HostData VTKIO)

    add_executable(DemoModelNxd
        "${PROJECT_SOURCE_DIR}/source/Legacy/DemoModelNxd.cpp"
        )
    target_link_libraries(DemoModelNxd PanNDE_interfaces Controller_interfaces)
    target_link_libraries(DemoModelNxd HostData HostSolver NetMPI VTKIO)
    target_link_libraries(DemoModelNxd Stubs)
    add_executable(DemoModel
        "${PROJECT_SOURCE_DIR}/source/Legacy/DemoModel.cpp"
        )
    target_link_libraries(DemoModel PanNDE_interfaces Controller_interfaces)
    target_link_libraries(DemoModel HostData HostSolver NetMPI VTKIO)
    target_link_libraries(DemoModel Stubs)

    add_executable(DemoSPModel
        "${PROJECT_SOURCE_DIR}/source/Legacy/DemoSingleProcessModel.cpp"
        )
    target_link_libraries(DemoSPModel PanNDE_interfaces Controller_interfaces)
    target_link_libraries(DemoSPModel HostData HostSolver VTKIO)
    target_link_libraries(DemoSPModel Stubs)

    add_executable(DemoModel2xd
        "${PROJECT_SOURCE_DIR}/source/Legacy/DemoModel2xd.cpp"
        )
    target_link_libraries(DemoModel2xd PanNDE_interfaces Controller_interfaces)
    target_link_libraries(DemoModel2xd HostData HostSolver NetMPI VTKIO)
    target_link_libraries(DemoModel2xd Stubs)
endif()

