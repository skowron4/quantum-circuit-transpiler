cmake_minimum_required(VERSION 3.25)
project(magisterka LANGUAGES CXX)

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
        structures/swap_record.hpp
)

add_executable(magisterka ${SOURCES})

# Podpięcie zależności do projektu
target_link_libraries(magisterka PRIVATE project_deps)