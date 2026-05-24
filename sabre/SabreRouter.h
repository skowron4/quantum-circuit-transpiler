#ifndef MAGISTERKA_SABREROUTER_H
#define MAGISTERKA_SABREROUTER_H

#include <chrono>
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <string>
#include <sstream> // <- DODANE DO WIZUALIZACJI

#include "architecture/BaseTopology.h"
#include "differentialEvolution/differential_evolution.h"

inline double fast_max_sabre(double a, double b) {
    return (a > b) ? a : b;
}

struct SabreDAG {
    std::vector<GateNode> nodes;
    std::vector<int> initial_deps;
};

struct SwapCandidate {
    int p1, p2, r_shared;
    double score;
};

class SabreRouter {
public:
    SabreRouter(const BaseTopology& topology, const std::vector<Gate>& circuit, int qubits_cnt, std::string name);

    void run_benchmark();

private:
    // --- STAŁE KOSZTÓW I PARAMETRÓW ---
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
    SabreDAG dag;
    SabreDAG reversed_dag;

    // --- METODY GŁÓWNE ---
    int get_distance(int p1, int p2) const;
    SabreDAG build_dag(const std::vector<Gate>& circuit) const;
    std::vector<int> run_pass(const SabreDAG& current_dag, const std::vector<int>& initial_l2p, bool is_final);

    // --- METODY WALIDACJI I WYPISYWANIA (NOWE) ---
    bool validate_logical_circuit(const std::vector<Gate>& orig_circ, const std::vector<RoutingRecord>& log) const;
    void print_optimized_circuit(const std::vector<RoutingRecord>& log, int num_logical_qubits) const;

    // --- METODY POMOCNICZE DLA run_pass ---
    void force_qubit_down(int p, int current_gate_id, MachineState& state, std::vector<int>& q_logical_depth,
                          std::vector<RoutingRecord>& log, int& out_moves, double& total_routing_cost, bool is_final) const;

    void process_executable_gates(const SabreDAG& current_dag, const std::vector<int>& executable_gates,
                                  MachineState& state, std::vector<int>& q_logical_depth, std::vector<RoutingRecord>& log,
                                  int& out_moves, double& total_routing_cost, bool is_final,
                                  int& gates_routed, std::vector<int>& deps, std::vector<int>& front_layer) const;

    void update_active_and_extended_layers(const SabreDAG& current_dag, const std::vector<int>& front_layer,
                                           const MachineState& state, std::vector<int>& active_phys_qubits,
                                           std::vector<bool>& in_active, std::vector<int>& extended_layer,
                                           std::vector<bool>& in_extended) const;

    SwapCandidate find_best_swap_candidate(const SabreDAG& current_dag, const MachineState& state,
                                           const std::vector<int>& front_layer, const std::vector<int>& extended_layer,
                                           const std::vector<int>& active_phys_qubits, const std::vector<double>& decay) const;

    void apply_swap_candidate(const SwapCandidate& candidate, MachineState& state, std::vector<int>& q_logical_depth,
                              std::vector<RoutingRecord>& log, std::vector<double>& decay, int& out_moves, int& out_swaps,
                              double& total_routing_cost, bool is_final) const;
};

#endif //MAGISTERKA_SABREROUTER_H