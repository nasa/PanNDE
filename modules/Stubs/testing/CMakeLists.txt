#stubs tests


cmake_minimum_required (VERSION 3.13.5)

project(StubsTests C CXX)#CUDA

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib) #for static lib
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib) # for shared lib

enable_testing()

include_directories("${GTEST_SERVICE_DIR}/googletest/include"
                    "${GTEST_SERVICE_DIR}/googlemock/include")

add_executable(StubsTests
    "${PROJECT_SOURCE_DIR}/source/HostNormalXDZTests.cpp"
    "${GTEST_DRIVERS_DIR}/source/TestMain.cpp"
    )

target_link_libraries(StubsTests Stubs HostData PanNDE_interfaces)
#target_link_libraries(StubsTests PanNDE_interface_tests)
target_link_libraries(StubsTests gtest gtest_main)

target_include_directories(StubsTests PUBLIC "${PROJECT_SOURCE_DIR}/include")
