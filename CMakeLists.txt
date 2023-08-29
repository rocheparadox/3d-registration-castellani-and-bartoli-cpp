cmake_minimum_required(VERSION 3.16)
project(3d_resigtration_castellani_and_bartoli_cpp)

set(CMAKE_CXX_STANDARD 17)

add_executable(3d_resigtration_castellani_and_bartoli_cpp main.cpp
        icp.cpp
        icp.h
        utils.cpp
        utils.h)
