cmake_minimum_required(VERSION 3.25)
project(magisterka LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")

# 1. WYŁĄCZENIE MAKR MIN/MAX Z WINDOWS.H
add_compile_definitions(NOMINMAX)

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
        greedy/GreedyASAP.cpp
        greedy/GreedyASAP.h
        sabre/SabreRouter.cpp
        sabre/SabreRouter.h
        dEHybridSabre/SabreDEHybridRouter.cpp
        dEHybridSabre/SabreDEHybridRouter.h
)

add_executable(magisterka ${SOURCES})

# 2. POPRAWNA KONFIGURACJA OPENMP DLA MSVC
if(MSVC AND CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(magisterka PRIVATE /openmp:llvm)
endif()

# Podpięcie zależności do projektu
target_link_libraries(magisterka PRIVATE project_deps)

# 1. Wymuś optymalizacje czasu konsolidacji (LTO / IPO) - łączy cały kod z plików .cpp i inlinuje cross-modułowo
include(CheckIPOSupported)
check_ipo_supported(RESULT result)
if(result)
    set_target_properties(magisterka PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# 2. Agresywne flagi dla profilu Release
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # /O2 - Max speed, /Ob2 - Aggressive inlining, /arch:AVX2 - wektoryzacja SIMD, /fp:fast - szybsza matma zmiennoprzecinkowa
    target_compile_options(magisterka PRIVATE
            "$<$<CONFIG:RELEASE>:/O2;/Ob2;/Oi;/Ot;/arch:AVX2;/fp:fast>"
    )
else()
    # Flagi dla GCC / Clang
    target_compile_options(magisterka PRIVATE
            "$<$<CONFIG:RELEASE>:-O3;-march=native;-ffast-math;-funroll-loops>"
    )
endif()