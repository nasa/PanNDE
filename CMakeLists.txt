# PanNDE top-level CMake build file
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

#Plugin Architecture for Numerical NDE (Pan NDE)
project(PanNDE C CXX)#CUDA

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib) #for static lib
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib) # for shared lib
set(CMAKE_BUILD_TYPE RelWithDebInfo)
file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/data")
file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/testlogs")

option(GO_TEST "Enable Tests" ON)
option(GO_PROF "Enable Profiling" ON)
option(GO_LEGACY "Build Legacy Model Executors" OFF)
option(NO_WARN "Prevent Compile Warnings" OFF)

if(NO_WARN)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
endif()

if(GO_TEST)
    enable_testing()
    message(STATUS "-- TESTING ENABLED --")
    set(GTEST_SERVICE_DIR "${CMAKE_SOURCE_DIR}/modules/services/gtest")
    add_subdirectory("${GTEST_SERVICE_DIR}/googletest")
    include_directories("${GTEST_SERVICE_DIR}/googletest/include"
                        "${GTEST_SERVICE_DIR}/googlemock/include")
    set(GTEST_DRIVERS_DIR "${CMAKE_SOURCE_DIR}/modules/testDrivers")
else()
    message(WARNING "-- TESTING DISABLED --")
endif()

include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(GO_PROF)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

add_subdirectory("${CMAKE_SOURCE_DIR}/modules")
add_subdirectory("${CMAKE_SOURCE_DIR}/applications")
