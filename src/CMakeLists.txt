cmake_minimum_required(VERSION 3.1)
project(3DA_project_CageBasedDef)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use IMGUI"          ON)
#option(LIBIGL_WITH_PNG "Use PNG"          ON)

find_package(LIBIGL REQUIRED QUIET)

# Add project files
set(SRC_FILES ${PROJECT_SOURCE_DIR}/main.cpp
              MeanValueCoordController.hpp MeanValueCoordController.cpp
MeshProcessor.hpp MeshProcessor.cpp CageGenerator.hpp CageGenerator.cpp
DeformCageViewerPlugin.hpp)
              

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_executable(${PROJECT_NAME}_bin ${SRC_FILES})
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw igl::opengl_glfw_imgui)
