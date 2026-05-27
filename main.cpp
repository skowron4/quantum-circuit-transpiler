#include <iostream>
#include "generator/CircuitGenerator.h"
#include "geneticAlgorithm/toSort/swap_evaluator.hpp"
#include "architecture/IQMTriangleConstellation.h"
#include "dEHybridSabre/SabreDEHybridRouter.h"
#include "newDESABRE/newSabreDEHybridRouter.h"
#include "differentialEvolution/differential_evolution.h"
#include "sabre/SabreRouter.h"
#ifdef _WIN32
    #include <windows.h>
#endif

#include "greedy/GreedyASAP.h"


void run_test_for_topology(const BaseTopology &qpu, const std::string &name) {
    std::cout << "\n============================================\n";
    std::cout << "URUCHAMIAM TEST DLA TOPOLOGII: " << name << "\n";
    std::cout << "============================================\n";
    qpu.print_stats();

    int logical_qubits_count = qpu.get_num_qubits() * 1.0;
    CircuitGenerator generator(qpu, logical_qubits_count, 1234);

    // Generowanie obwodu testowego
    std::vector<Gate> circuit = generator.generate(2000, 0.6, 0.4);
    generator.print_gate_statistics(circuit);
    generator.print_circuit(circuit);

    // 1. Wykonanie algorytmu bazowego (Benchmark zachłanny ASAP)
    // GreedyASAP greedy_router(qpu, circuit, qpu.get_num_qubits(), name);
    // greedy_router.run_benchmark();

    // 2. Wykonanie algorytmu optymalizacyjnego (Differential Evolution)
    // DifferentialEvolution de_router(qpu, circuit, qpu.get_num_qubits(), name);
    // de_router.optimize(100, 2000000);
    // de_router.optimize_windowed(100, 30000, 200, 1.0, qpu.get_num_qubits()); //30, 1000, 100, 0.35

    // 3. Wykonanie algorytmu SABRE
    SabreRouter sabre(qpu, circuit, qpu.get_num_qubits(), name);
    sabre.run_benchmark();

    // 4. Wykonanie algorytmu DEHybridSabre
    SabreDEHybridRouter SabreDERouter(qpu, circuit, qpu.get_num_qubits(), name, 50);
    SabreDERouter.optimize_windowed(100, 5000, 150, 1.0, false, 0.9);

    // 5. Wykonanie algorytmu newDEHybridSabre
    // newSabreDEHybridRouter newSabreDERouter(qpu, circuit, qpu.get_num_qubits(), name, 4);
    // newSabreDERouter.optimize_windowed(30, 100, 500, 0.7, true, false, 0.8, 0.9);

}
void run_test_for_topologyII(const BaseTopology &qpu, const std::string &name) {
    std::cout << "\n============================================\n";
    std::cout << "URUCHAMIAM TEST DLA TOPOLOGII: " << name << "\n";
    std::cout << "============================================\n";
    qpu.print_stats();

    int logical_qubits_count = qpu.get_num_qubits() * 1.0;
    CircuitGenerator generator(qpu, logical_qubits_count, 1234);

    std::vector<Gate> circuit = generator.generate(200, 0.6, 0.4);
    generator.print_gate_statistics(circuit);
    generator.print_circuit(circuit);

    SabreRouter sabre(qpu, circuit, qpu.get_num_qubits(), name);
    sabre.run_benchmark();

    std::cout << "\n--- Rozpoczynam 10 iteracji dla DEHybridSabre ---\n";

    const int NUM_RUNS = 5;
    std::vector<OptimizationResult> results;
    results.reserve(NUM_RUNS);

    for (int i = 0; i < NUM_RUNS; ++i) {
        SabreDEHybridRouter SabreDERouter(qpu, circuit, qpu.get_num_qubits(), name, false, 15);

        // Odbieramy całą strukturę wyników
        OptimizationResult res = SabreDERouter.optimize_windowed(30, 300, 500, 0.7,  false, 0.9);
        results.push_back(res);

        std::cout << "Iteracja " << (i + 1) << " zakonczona. "
                  << "Wynik: " << res.score << ", "
                  << "SWAPy: " << res.swaps << ", "
                  << "Makespan: " << res.makespan << "\n";
    }

    // --- OBLICZANIE STATYSTYK ---
    OptimizationResult min_val = results[0];
    OptimizationResult max_val = results[0];
    OptimizationResult sum_val = results[0];

    for (int i = 1; i < NUM_RUNS; ++i) {
        min_val.score = std::min(min_val.score, results[i].score);
        max_val.score = std::max(max_val.score, results[i].score);
        sum_val.score += results[i].score;

        min_val.makespan = std::min(min_val.makespan, results[i].makespan);
        max_val.makespan = std::max(max_val.makespan, results[i].makespan);
        sum_val.makespan += results[i].makespan;

        min_val.depth = std::min(min_val.depth, results[i].depth);
        max_val.depth = std::max(max_val.depth, results[i].depth);
        sum_val.depth += results[i].depth;

        min_val.swaps = std::min(min_val.swaps, results[i].swaps);
        max_val.swaps = std::max(max_val.swaps, results[i].swaps);
        sum_val.swaps += results[i].swaps;

        min_val.moves = std::min(min_val.moves, results[i].moves);
        max_val.moves = std::max(max_val.moves, results[i].moves);
        sum_val.moves += results[i].moves;
    }

    std::cout << "\n============================================\n";
    std::cout << " PODSUMOWANIE DEHybridSabre (" << NUM_RUNS << " uruchomien)\n";
    std::cout << "============================================\n";

    std::cout << " [SCORE]\n";
    std::cout << " Min: " << min_val.score << " | Max: " << max_val.score << " | Średnia: " << (sum_val.score / NUM_RUNS) << "\n\n";

    std::cout << " [MAKESPAN]\n";
    std::cout << " Min: " << min_val.makespan << " | Max: " << max_val.makespan << " | Średnia: " << (sum_val.makespan / NUM_RUNS) << "\n\n";

    std::cout << " [DEPTH]\n";
    std::cout << " Min: " << min_val.depth << " | Max: " << max_val.depth << " | Średnia: " << (static_cast<double>(sum_val.depth) / NUM_RUNS) << "\n\n";

    std::cout << " [SWAPS]\n";
    std::cout << " Min: " << min_val.swaps << " | Max: " << max_val.swaps << " | Średnia: " << (static_cast<double>(sum_val.swaps) / NUM_RUNS) << "\n\n";

    std::cout << " [MOVES]\n";
    std::cout << " Min: " << min_val.moves << " | Max: " << max_val.moves << " | Średnia: " << (static_cast<double>(sum_val.moves) / NUM_RUNS) << "\n";
    std::cout << "============================================\n";
}


int main() {
#ifdef _WIN32
    // Windows potrzebuje instrukcji zmiany kodowania na UTF-8
    SetConsoleOutputCP(CP_UTF8);
    // Opcjonalnie: upewnia się, że wejście (cin) też używa UTF-8
    SetConsoleCP(CP_UTF8);
#endif
    // Topologia siatki 2x2
    // IQMConstellation gridQpu(5, 5);
    // run_test_for_topology(gridQpu, "Siatka 3x3");
    // run_test_for_topologyII(gridQpu, "Siatka 5x5");

    // Topologia trójkątna
    IQMTriangleConstellation triangleQpu;
    // run_test_for_topology(triangleQpu, "Trójkąt 3-Rezonatory");
    run_test_for_topologyII(triangleQpu, "Trójkąt 3-Rezonatory");

    return 0;
}
