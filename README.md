# Datasets
## Type:  Multi-layers Graph

## example graph name: AUCS

In this graph, the multiple layers represent relationships between 61 employees of a University department in five different aspects: (i) coworking, (ii) having lunch together, (iii) Facebook friendship, (iv) offline friendship (having fun together), and (v) coauthor-ship.




## Preprocess
goal:Write Binary file From edge list
- format: edge list without n or m or other informations (graph.txt)
- function: write_graph_binary(graph_name)


# Info for files

## sequential_influential_community

This folder contains:
- ClimbStairs algorithm
- Divide and conquer algorithm

## parallel_influential_community

This folder contains:
- Basic parallel algorithm based on block
- Advanced parallel algorithm based on thread pool


# Structure of directory

以下是为你的项目 density_decomposition 编写的 CMakeLists.txt 文件。它支持组织清晰的项目结构，并正确包含和链接子目录中的文件：
CMakeLists.txt

# Minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Project name and version
project(density_decomposition VERSION 1.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Enable warnings and debug symbols
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-Wall -Wextra -Wpedantic -g)
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
)

# Add the executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Optional: Link libraries (if any external libraries are used)
# target_link_libraries(${PROJECT_NAME} ...)

目录结构说明

假设你的项目目录结构如下：

density_decomposition/
├── CMakeLists.txt
├── main.cpp
├── utilities/
│   ├── Defines.h
│   ├── Log.cpp
│   ├── Log.h
│   ├── ThreadPool.h
│   ├── Time.h
│   ├── Tools.h
│   ├── Utility.h
├── graph/
│   ├── Graph.cpp
│   ├── Graph.h
│   ├── MultilayerGraph.cpp
│   ├── MultilayerGraph.h
