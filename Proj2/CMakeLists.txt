cmake_minimum_required(VERSION 3.10)
Project(mesh_demo)
add_executable(mesh_demo mesh_demo.cpp)
target_link_libraries(mesh_demo PUBLIC mesh_lib)
add_subdirectory(mesh_lib)
target_include_directories(mesh_demo PUBLIC
	"${PROJECT_BINARY_DIR}"
	"${PROJECT_SOURCE_DIR}/mesh_lib")
