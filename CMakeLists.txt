cmake_minimum_required(VERSION 3.19)
cmake_policy(VERSION 3.19[...3.19])

project(origami-simulator)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS Off)
set(VORPALINE_PLATFORM Linux64-gcc-dynamic)

# libigl
option(LIBIGL_WITH_OPENGL                  "Use OpenGL"          ON)
option(LIBIGL_WITH_OPENGL_GLFW             "Use GLFW"            ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI       "Use IMGUI"           ON)
option(LIBIGL_WITH_PNG                     "Use PNG"             ON)

# include directories
include_directories(${CMAKE_SOURCE_DIR}/autodiff)
include_directories(${CMAKE_SOURCE_DIR}/include)

message(STATUS ${CMAKE_MODULE_PATH})
include(libigl)

find_package(OpenMP)

# add external libraries
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")

# Add your project files
file(GLOB SRCFILES *.cpp)
file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB CPP_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# build executable
add_executable(${PROJECT_NAME} main.cpp ${SRCFILES} ${H_FILES} ${CPP_FILES})
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui igl::png OpenMP::OpenMP_CXX)