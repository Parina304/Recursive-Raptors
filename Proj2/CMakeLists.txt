cmake_minimum_required(VERSION 3.12...3.29)

Project(mesh_demo)


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


# add executable
add_executable(mesh_demo mesh_demo.cpp)

# Include directories (Add mesh_lib folder to search path for headers)
# include_directories(${CMAKE_SOURCE_DIR}/mesh_lib)
# include_directories(${CMAKE_SOURCE_DIR}/mesh_merge)

include_directories("${PROJECT_BINARY_DIR}"
                    "${PROJECT_SOURCE_DIR}/mesh_lib"
                    "${PROJECT_SOURCE_DIR}/mesh_merge"
                    )

add_subdirectory(mesh_merge)
add_subdirectory("${CMAKE_SOURCE_DIR}/mesh_lib")

# link with mesh_lib mesh_merge and CGAL
target_link_libraries(mesh_demo PUBLIC mesh_lib mesh_merge CGAL::CGAL)
