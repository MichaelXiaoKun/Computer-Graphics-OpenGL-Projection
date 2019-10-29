cmake_minimum_required(VERSION 3.10.2)
project(ecs175hw2)
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})

find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIR})

set(CMAKE_CXX_FLAG "-g -Wall")
set(CMAKE_CXX_STANDARD 14)

add_executable(ecs175hw2 main.cpp)
target_link_libraries(ecs175hw2 ${GLUT_LIBRARY} ${OPENGL_LIBRARY} m)