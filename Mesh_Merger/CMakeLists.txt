cmake_minimum_required(VERSION 3.12...3.29)

project(test5)

# Find CGAL
find_package(CGAL QUIET COMPONENTS)
if (NOT CGAL_FOUND)
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()
endif()

# Find Boost
find_package(Boost REQUIRED)
if (NOT Boost_FOUND)
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

# Include directories (Add mesh_lib folder to search path for headers)
include_directories(${CMAKE_SOURCE_DIR}/mesh_lib)

# Add mesh_lib as a library
add_library(mesh_lib STATIC mesh_lib/mesh_lib.cpp)

# Add test5 executable
add_executable(test5 test5.cpp)

# Link test5 with mesh_lib and CGAL
target_link_libraries(test5 mesh_lib CGAL::CGAL)
