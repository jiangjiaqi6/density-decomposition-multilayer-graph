# Minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Project name and version
# project(den_decom_vsht VERSION 1.0)
# project(den_decom_hpdd VERSION 1.0)
project(den_decom_hpdd_plus VERSION 1.0)



# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_BUILD_TYPE Debug)




# Enable warnings and debug symbols
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-w -g)
endif()



if(PROJECT_NAME STREQUAL "den_decom_vsht")
    add_definitions(-DVSHT)
endif()


if(PROJECT_NAME STREQUAL "den_decom_hpdd")
    add_definitions(-DHPDD)
endif()

if(PROJECT_NAME STREQUAL "den_decom_hpdd_plus")
    add_definitions(-DHPDDPLUS)
endif()



# Add source directories
include_directories(utilities)
include_directories(graph)

# List of source files
set(SOURCES
    main.cpp
    utilities/Log.cpp
    graph/Graph.cpp
    graph/MultilayerGraph.cpp
    core/Reorientation.cpp
)

# Add the executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Optional: Link libraries (if any external libraries are used)
# target_link_libraries(${PROJECT_NAME} ...)
