cmake_minimum_required(VERSION 4.0)
project(magisterka)

set(CMAKE_CXX_STANDARD 17)

add_executable(magisterka main.cpp
        architecture/IQMConstellation.cpp
        architecture/IQMConstellation.h
        architecture/IQMTriangleConstellation.cpp
        architecture/IQMTriangleConstellation.h
        Generator/CircuitGenerator.cpp
        Generator/CircuitGenerator.h)
