cmake_minimum_required(VERSION 3.8)
project(pspy)

set(CMAKE_CXX_STANDARD 14)

find_package(Eigen3 REQUIRED)

add_subdirectory(parasolid)

add_library(pspy   "disjointset.cpp" "disjointset.h" "lsh.h" "lsh.cpp" "eclass.h" "eclass.cpp" "types.h" "face.h" "face.cpp" "loop.h" "loop.cpp" "edge.h" "edge.cpp" "body.h" "body.cpp" "topology.h" "vertex.h" "vertex.cpp" "part.h" "part.cpp")

target_link_libraries(pspy PUBLIC parasolid Eigen3::Eigen)

target_include_directories(pspy INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
