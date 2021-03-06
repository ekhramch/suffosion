cmake_minimum_required(VERSION 2.8)

project(suffosion)

include_directories(${CMAKE_SOURCE_DIR})

#----------------------------------------------------------------------------
# Find AMGCL
#----------------------------------------------------------------------------
find_package(amgcl)
find_package(VexCL)

# add the executable
vexcl_add_executables(suffosion 
    suffosion.cpp
    utils.cpp
    parameters.cpp
    saver.cpp
    pressure.cpp
    stress.cpp
    chemistry.cpp
    )

target_link_libraries(
    suffosion INTERFACE
    ${Boost_LIBRARIES}
    hdf5 
    hdf5_cpp 
    pugixml   
    amgcl::amgcl
    )
