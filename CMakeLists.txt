cmake_minimum_required(VERSION 2.8)

project(suffosion)

include_directories(${CMAKE_SOURCE_DIR})

#----------------------------------------------------------------------------
# Find Boost
#----------------------------------------------------------------------------
set(BOOST_COMPONENTS
    chrono
    system
    program_options
    filesystem
    )
find_package(Boost COMPONENTS ${BOOST_COMPONENTS})
include_directories(${Boost_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Find AMGCL
#----------------------------------------------------------------------------
set(AMGCL_ROOT $ENV{AMGCL_ROOT} CACHE STRING "AMGCL root")
message(STATUS "Found AMGCL at ${AMGCL_ROOT}")
include_directories(${AMGCL_ROOT})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

#----------------------------------------------------------------------------
# Find OpenMP
#----------------------------------------------------------------------------
find_package(OpenMP)

if (OpenMP_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

set(VEXCL_BACKEND "CUDA" CACHE STRING "Select VexCL backend (OpenCL/CUDA)")
set_property(CACHE VEXCL_BACKEND PROPERTY STRINGS "OpenCL" "CUDA" "Compute")

#----------------------------------------------------------------------------
# Find Backend
#----------------------------------------------------------------------------
if ("${VEXCL_BACKEND}" STREQUAL "OpenCL")
    set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
    find_package(OpenCL REQUIRED)
    include_directories( ${OPENCL_INCLUDE_DIRS} )
    set(BACKEND_LIBS ${OPENCL_LIBRARIES})
    add_definitions(-DVEXCL_BACKEND_OPENCL)
elseif ("${VEXCL_BACKEND}" STREQUAL "Compute")
    find_boost_compute()
    set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
    find_package(OpenCL REQUIRED)
    include_directories( ${OPENCL_INCLUDE_DIRS} )
    set(BACKEND_LIBS ${OPENCL_LIBRARIES})
    add_definitions(-DVEXCL_BACKEND_COMPUTE)
elseif ("${VEXCL_BACKEND}" STREQUAL "CUDA")
    find_package(CUDA REQUIRED)
    include_directories( ${CUDA_INCLUDE_DIRS} )
    set(BACKEND_LIBS ${CUDA_CUDA_LIBRARY})
    add_definitions(-DVEXCL_BACKEND_CUDA)
endif()

# add the executable
add_executable(suffosion 
    suffosion.cpp
    pr_util.cpp
    saver.cpp
    )

target_link_libraries(
    suffosion
    ${Boost_LIBRARIES}
    hdf5 
    hdf5_cpp 
    pugixml   
    cuda
    )
