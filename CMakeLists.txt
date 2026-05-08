cmake_minimum_required(VERSION 4.0)
project(magisterka)

set(CMAKE_CXX_STANDARD 17)

add_executable(magisterka main.cpp
        IQMConstellation.cpp
        IQMConstellation.h
        IQMTriangleConstellation.cpp
        IQMTriangleConstellation.h
        CircuitGenerator.cpp
        CircuitGenerator.h)
