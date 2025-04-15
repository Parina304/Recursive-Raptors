cmake_minimum_required(VERSION 3.12...3.29)

project(glutviz)
project(demo_glut)

# Find OpenGL
find_package(OpenGL REQUIRED)

# Find FreeGLUT (instead of standard GLUT)
find_package(GLUT REQUIRED)

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

# Include directories for OpenGL and FreeGLUT
include_directories(${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})

# Include additional directories
include_directories("${PROJECT_BINARY_DIR}"
                    "${PROJECT_SOURCE_DIR}/../mesh_lib"
                    "${PROJECT_SOURCE_DIR}/../mesh_merge"
                    )
add_subdirectory("${CMAKE_SOURCE_DIR}/../mesh_merge" mesh_merge)
add_subdirectory("${CMAKE_SOURCE_DIR}/../mesh_lib" mesh_lib)

# Add the executable
add_executable(glutviz glutviz.cpp)
add_executable(demo_glut demo_glut.cpp)

# Link FreeGLUT and OpenGL libraries
target_link_libraries(glutviz ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} mesh_lib mesh_merge CGAL::CGAL)
target_link_libraries(demo_glut ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
