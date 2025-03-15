cmake_minimum_required(VERSION 3.10)
Project(operDemo)
add_executable(operDemo operDemo.cpp)
target_link_libraries(operDemo PUBLIC mesh_lib)
add_subdirectory(src)
target_include_directories(operDemo PUBLIC
	"${PROJECT_BINARY_DIR}"
	"${PROJECT_SOURCE_DIR}/inc")
