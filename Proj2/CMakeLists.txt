cmake_minimum_required(VERSION 3.12...3.29)

Project(mesh_demo)
Project(mesh_demo_union)
Project(mesh_demo_union)
Project(vizmesh)
Project(VolMesh)
Project(VTKtoOBJ)
Project(VolMesh_old)

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

# find GLUT
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(cppcolormap REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )


# add executable
add_executable(mesh_demo mesh_demo.cpp)
add_executable(mesh_demo_union mesh_demo_union.cpp)
add_executable(mesh_demo_vol mesh_demo_vol.cpp)
add_executable(vizmesh vizmesh.cpp)
add_executable(VolMesh VolMesh.cpp)
add_executable(VTKtoOBJ VTKtoOBJ.cpp)
add_executable(VolMesh_old VolMesh_old.cpp)

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
target_link_libraries(mesh_demo_union PUBLIC mesh_lib mesh_merge CGAL::CGAL)
target_link_libraries(mesh_demo_vol PUBLIC mesh_lib mesh_merge CGAL::CGAL)

target_link_libraries(vizmesh ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${cppcolormap_LIBRARY})
