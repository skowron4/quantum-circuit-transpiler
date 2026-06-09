#ifndef MAGISTERKA_SABREROUTER_H
#define MAGISTERKA_SABREROUTER_H

#include <chrono>
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <string>
#include <sstream>
#include <random>
#include <cstdint> // <- DODANE DLA uint8_t

#include "architecture/BaseTopology.h"
#include "differentialEvolution/differential_evolution.h"

inline double fast_max_sabre(double a, double b) {
    return (a > b) ? a : b;
}

struct LogicalPair {
    int lq1;
    int lq2;
    double weight; // Waga zależna od głębokości warstwy
};

struct SabreDAG {
    std::vector<GateNode> nodes;
    std::vector<int> initial_deps;
};

struct SwapCandidate {
    int p1, p2, r_shared;
    double score;
};

struct SabreBenchmarkResult {
    double score = 0.0;
    double makespan = 0.0;
    double execution_time = 0.0;
    int depth = 0;
    int swaps = 0;
    int moves = 0;
};

// NOWE: Struktura prekompilowanych krawędzi
struct stockSwapEdge {
    int neighbor;
    int resonator;
};

class SabreRouter {
public:
    SabreRouter(const BaseTopology& topology, const std::vector<Gate>& circuit, int qubits_cnt, std::string name);

    SabreBenchmarkResult run_benchmark(bool print_details = true);
    void run_benchmark_multiple(int num_runs = 30);
    void set_random_seed(uint32_t seed);

private:
    struct SabrePassResult {
        std::vector<int> mapping;
        SabreBenchmarkResult benchmark;
    };

    static constexpr int RES_OFFSET = 1000;
    static constexpr double EXTENDED_LAYER_WEIGHT = 0.5;
    static constexpr double MOVE_COST = 1.0;
    static constexpr double SWAP_COST = 3.14;
    static constexpr double DECAY_RESET = 1.0;
    static constexpr double DECAY_INCREMENT = 0.2;

    const BaseTopology& topo;
    std::vector<Gate> original_circuit;
    int num_qubits;
    std::string name;

    mutable std::vector<int> dist_cache;
    mutable std::mt19937 rng;
    std::vector<std::vector<stockSwapEdge>> precomputed_swaps; // NOWE
    SabreDAG dag;
    SabreDAG reversed_dag;

    int get_distance(int p1, int p2) const;
    std::vector<int> create_random_initial_mapping();
    SabreDAG build_dag(const std::vector<Gate>& circuit) const;
    SabrePassResult run_pass(const SabreDAG& current_dag, const std::vector<int>& initial_l2p, bool is_final, bool print_details = true);

    bool validate_logical_circuit(const std::vector<Gate>& orig_circ, const std::vector<RoutingRecord>& log, bool print_success = true) const;
    void print_optimized_circuit(const std::vector<RoutingRecord>& log, int num_logical_qubits) const;

    void force_qubit_down(int p, int current_gate_id, MachineState& state, std::vector<int>& q_logical_depth,
                      std::vector<RoutingRecord>& log, int& out_moves, double& total_routing_cost, bool is_final,
                      bool allow_relocation = false, const SabreDAG* current_dag = nullptr, const std::vector<int>* front_layer = nullptr) const;

    void process_executable_gates(const SabreDAG& current_dag, const std::vector<int>& executable_gates,
                                  MachineState& state, std::vector<int>& q_logical_depth, std::vector<RoutingRecord>& log,
                                  int& out_moves, double& total_routing_cost, bool is_final,
                                  int& gates_routed, std::vector<int>& deps, std::vector<int>& front_layer) const;

    // ZMIENIONE: uzycie uint8_t
    void update_active_and_extended_layers(const SabreDAG& current_dag, const std::vector<int>& front_layer,
                                           const MachineState& state, std::vector<int>& active_phys_qubits,
                                           std::vector<uint8_t>& in_active, std::vector<int>& extended_layer,
                                           std::vector<uint8_t>& in_extended) const;

    SwapCandidate find_best_swap_candidate(
    const SabreDAG& current_dag,
    const MachineState& state,
    const std::vector<int>& front_layer,
    const std::vector<std::vector<int>>& lookahead_layers,
    const std::vector<int>& active_phys_qubits,
    const std::vector<double>& decay,
    int max_lookahead_depth = 5,
    double lookahead_weight_base = 0.5,
    double lookahead_decay_rate = 0.5
) const;

    void apply_swap_candidate(const SwapCandidate& candidate, MachineState& state, std::vector<int>& q_logical_depth,
                              std::vector<RoutingRecord>& log, std::vector<double>& decay, int& out_moves, int& out_swaps,
                              double& total_routing_cost, bool is_final) const;
};

#endif //MAGISTERKA_SABREROUTER_H