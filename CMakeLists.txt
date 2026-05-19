cmake_minimum_required(VERSION 3.25)
project(magisterka LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")

# --- Zależności (Pakiety) ---
find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(boost_random CONFIG REQUIRED)
find_package(boost_align CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(absl CONFIG REQUIRED)

# --- Biblioteka Interfejsowa (Twoje zależności) ---
add_library(project_deps INTERFACE)

target_link_libraries(project_deps INTERFACE
        OpenMP::OpenMP_CXX
        Boost::random
        Boost::align
        fmt::fmt-header-only
        absl::linked_hash_set
)

target_include_directories(project_deps INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# --- Plik Wykonywalny ---
set(SOURCES
        main.cpp
        architecture/IQMConstellation.cpp
        architecture/IQMConstellation.h
        architecture/IQMTriangleConstellation.cpp
        architecture/IQMTriangleConstellation.h
        Generator/CircuitGenerator.cpp
        Generator/CircuitGenerator.h
        geneticAlgorithm/toSort/swap_evaluator.hpp
        structures/physical_qubit.hpp
        geneticAlgorithm/toSort/resonator_mute.hpp
        geneticAlgorithm/toSort/gate_mute.hpp
        geneticAlgorithm/toSort/resonator_gate_mute.hpp
        geneticAlgorithm/toSort/ox_pair.hpp
        structures/routing_record.hpp
        differentialEvolution/differential_evolution.cpp
        differentialEvolution/differential_evolution.h
        architecture/BaseTopology.cpp
        architecture/BaseTopology.h
)

add_executable(magisterka ${SOURCES})
if(MSVC AND CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(magisterka PRIVATE /openmp:llvm)
endif()

# Podpięcie zależności do projektu
target_link_libraries(magisterka PRIVATE project_deps)
target_compile_options(magisterka PRIVATE -fopenmp)
target_link_options(magisterka PRIVATE -fopenmp)