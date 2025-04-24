cmake_minimum_required(VERSION 3.12...3.29)

project(glutviz)
project(demo_glut)
project(visuals_3)

# Find OpenGL
find_package(OpenGL REQUIRED)

# Find FreeGLUT (instead of standard GLUT)
find_package(GLUT REQUIRED)

# find_package(GLFW REQUIRED)
find_package(glfw3 3.4 REQUIRED)

find_package(glm REQUIRED)

# find_package(imgui REQUIRED)

find_package(cppcolormap REQUIRED)

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

add_library(imgui
	imgui/imgui.cpp
	imgui/imgui_draw.cpp
	imgui/imgui_tables.cpp
	imgui/imgui_widgets.cpp
	imgui/backends/imgui_impl_glfw.cpp
	imgui/backends/imgui_impl_opengl3.cpp
)

add_library(implot
	implot/implot.cpp
	implot/implot_items.cpp)

# Include directories for OpenGL and FreeGLUT
include_directories(${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS} ${GLM_INCLUDE_DIRS})
include_directories("/home/heartman/miniconda3/include/")
include_directories("${PROJECT_SOURCE_DIR}/glad/")
include_directories("${PROJECT_SOURCE_DIR}/imgui/")
include_directories("${PROJECT_SOURCE_DIR}/imgui/backends/")
include_directories("${PROJECT_SOURCE_DIR}/implot/")



# Include additional directories
include_directories("${PROJECT_SOURCE_DIR}/")
include_directories("${PROJECT_BINARY_DIR}"
                    "${PROJECT_SOURCE_DIR}/../mesh_lib"
                    "${PROJECT_SOURCE_DIR}/../mesh_merge"
                    )
add_subdirectory("${CMAKE_SOURCE_DIR}/../mesh_merge" mesh_merge)
add_subdirectory("${CMAKE_SOURCE_DIR}/../mesh_lib" mesh_lib)

# Add the executable
# add_executable(glutviz glutviz.cpp)
# add_executable(demo_glut demo_glut.cpp)
add_executable(visuals_3 visuals_3.cpp "${CMAKE_SOURCE_DIR}/glad/glad.c")

# Link FreeGLUT and OpenGL libraries
# target_link_libraries(glutviz ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${cppcolormap_LIBRARY} mesh_lib mesh_merge CGAL::CGAL)
# target_link_libraries(demo_glut ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${cppcolormap_LIBRARY})
target_link_libraries(visuals_3 ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${GLFW_LIBRARIES} ${GLM_LIBRARIES} ${cppcolormap_LIBRARY} glfw imgui::imgui implot mesh_lib mesh_merge CGAL::CGAL)
