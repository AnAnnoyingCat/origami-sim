cmake_minimum_required(VERSION 3.19)
cmake_policy(VERSION 3.19[...3.19])

project(origami-simulator)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS Off)
set(VORPALINE_PLATFORM Linux64-gcc-dynamic)

set(LIBIGL_WITH_OPENGL ON CACHE BOOL "" FORCE)
set(LIBIGL_WITH_OPENGL_GLFW ON CACHE BOOL "" FORCE)
set(LIBIGL_WITH_PREDICATES ON CACHE BOOL "" FORCE)

include(FetchContent)

set(IPC_TOOLKIT_GIT_TAG v1.3.1)

FetchContent_Declare(
    ipc_toolkit
    GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
    GIT_TAG ${IPC_TOOLKIT_GIT_TAG}
)
FetchContent_MakeAvailable(ipc_toolkit)
include_directories(${ipc_toolkit_SOURCE_DIR}/src)

# include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/external/include)

find_package(OpenMP)

# add external libraries
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")

# Add your project files
file(GLOB SRCFILES *.cpp)
file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB CPP_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# build executable
add_executable(${PROJECT_NAME} main.cpp ${SRCFILES} ${H_FILES} ${CPP_FILES})

target_link_libraries(${PROJECT_NAME} ipc::toolkit glfw glad OpenMP::OpenMP_CXX)

