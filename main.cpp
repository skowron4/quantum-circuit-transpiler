#include <iostream>

#include "generator/CircuitGenerator.h"
#include "architecture/IQMConstellation.h"
#include "geneticAlgorithm/toSort/swap_evaluator.hpp"
#include "architecture/IQMTriangleConstellation.h"
#include "differentialEvolution/differential_evolution.h"
#include <windows.h>

#include "greedy/GreedyASAP.h"


void run_test_for_topology(const BaseTopology &qpu, const std::string &name) {
    std::cout << "\n============================================\n";
    std::cout << "URUCHAMIAM TEST DLA TOPOLOGII: " << name << "\n";
    std::cout << "============================================\n";
    qpu.print_stats();

    int logical_qubits_count = qpu.get_num_qubits() * 1.0;
    CircuitGenerator generator(qpu, logical_qubits_count, 1234);

    // Generowanie obwodu testowego
    std::vector<Gate> circuit = generator.generate(50, 0.6, 0.4);
    generator.print_gate_statistics(circuit);
    generator.print_circuit(circuit);

    // 1. Wykonanie algorytmu bazowego (Benchmark zachłanny ASAP)
    GreedyASAP greedy_router(qpu, circuit, qpu.get_num_qubits(), name);
    greedy_router.run_benchmark();

    // 2. Wykonanie algorytmu optymalizacyjnego (Differential Evolution)
    DifferentialEvolution de_router(qpu, circuit, qpu.get_num_qubits(), name);
    de_router.optimize(100, 500);
}

int main() {
    SetConsoleOutputCP(CP_UTF8);
    // Topologia siatki 2x2
    IQMConstellation gridQpu(3, 3);
    run_test_for_topology(gridQpu, "Siatka 3x3");

    // Topologia trójkątna
    IQMTriangleConstellation triangleQpu;
    run_test_for_topology(triangleQpu, "Trójkąt 3-Rezonatory");

    return 0;
}
