//
// Created by skowr on 23.05.2026.
//

#ifndef MAGISTERKA_SABREDEHYBRIDROUTER_H
#define MAGISTERKA_SABREDEHYBRIDROUTER_H

#pragma once
#include "architecture/BaseTopology.h"
#include "structures/Gate.h"
#include "structures/routing_record.hpp"
#include <vector>
#include <string>
#include <limits>
#include <boost/random/xoshiro.hpp>

// --- Struktury wewnątrz hybrydy ---
struct HybridGateNode {
    int id;
    int lq1, lq2;
    int q1_1q_count, q2_1q_count;
    int deps_count = 0;
    std::vector<int> next_gates;
};

struct HybridSabreDAG {
    std::vector<HybridGateNode> nodes;
    std::vector<int> initial_deps;
};

// Zoptymalizowany HybridWorkspace (Usunięto SwapCand, dodano flagi O(1) i pre-alokowane pary)
struct HybridWorkspace {
    std::vector<int> l2p, p2l, q_on_res, res_occ, q_depth, deps, ready_queue, shared_res;
    std::vector<double> q_time, r_time;

    // Zbuforowane kontenery
    std::vector<int> front_layer, exec_gates, ext_layer;

    // --- NOWE: Optymalizacja BFS (Licznik Epok) ---
    // Zastępuje std::vector<bool> in_ext;
    int current_bfs_epoch = 0;
    std::vector<int> visited_epoch;

    std::vector<double> decay;

    std::vector<uint8_t> in_active;
    std::vector<int> active_phys;

    // NOWE POLA: optymalizacja O(1) dla krawędzi w routingu i wyciągnięte do góry lokacje
    int route_step_counter = 0;
    std::vector<int> seen_edges;
    std::vector<std::pair<int, int>> f_locs;
    std::vector<int> flat_lookahead_layers;
    std::vector<int> lookahead_offsets;
    std::vector<std::pair<int, int>> flat_lookahead_locs;
    std::vector<int> lookahead_locs_offsets;
};

// Stan maszyny wyłącznie dla przebiegu SABRE
struct LocalMachineState {
    std::vector<int> phys_qubit_to_logical;
    std::vector<int> logical_to_phys_qubit;
    std::vector<int> q_on_resonator;
    std::vector<int> res_occupant;
    std::vector<double> q_ready_time;
    std::vector<double> r_ready_time;
};

// Kontekst dla heurystyki SABRE
struct SabreRunContext {
    LocalMachineState state;
    std::vector<int> q_logical_depth;
    std::vector<double> decay;
    int out_swaps = 0;
    int out_moves = 0;
    double total_routing_cost = 0.0;
};

// Kontekst dla ewaluacji DE oparty o pre-alokowany workspace
struct DERunContext {
    std::vector<int>& l2p;
    std::vector<int>& p2l;
    std::vector<int>& q_on_res;
    std::vector<int>& res_occ;
    std::vector<int>& q_depth;
    std::vector<double>& q_time;
    std::vector<double>& r_time;
    int out_swaps = 0;
    int out_moves = 0;
    double total_routing_cost = 0.0;
    HybridWorkspace& ws;

    DERunContext(HybridWorkspace& w)
        : l2p(w.l2p), p2l(w.p2l), q_on_res(w.q_on_res),
          res_occ(w.res_occ), q_depth(w.q_depth),
          q_time(w.q_time), r_time(w.r_time), ws(w) {}
};

// --- Deklaracja klasy ---
class SabreDEHybridRouter {
private:
    static constexpr int RES_OFFSET = 1000;
    static constexpr double MOVE_COST = 1.0;
    static constexpr double SWAP_COST = 3.14;
    static constexpr double EXTENDED_LAYER_WEIGHT = 0.5;

    const BaseTopology& topo;
    std::vector<Gate> original_circuit;
    int num_qubits;
    std::string name;
    int max_look_ahead_depth;

    mutable std::vector<int> dist_cache;
    HybridSabreDAG dag;
    HybridSabreDAG reversed_dag;
    HybridWorkspace workspace;
    std::vector<int> saved_initial_l2p;

    // --- Metody Bazowe ---
    int get_distance(int p1, int p2) const;
    HybridSabreDAG build_dag(const std::vector<Gate>& circuit) const;

    // --- Metody Pomocnicze SABRE ---
    std::vector<int> run_sabre_pass(const HybridSabreDAG& current_dag, const std::vector<int>& initial_l2p, bool is_final, std::vector<int>* out_logical_order);
    void sabre_force_qubit_down(int p, int current_gate_id, SabreRunContext& ctx);
    void sabre_execute_gates(const HybridSabreDAG& current_dag, std::vector<int>& executable_gates, std::vector<int>& front_layer, std::vector<int>& deps, SabreRunContext& ctx, std::vector<int>* out_logical_order);
    bool sabre_route_step(const HybridSabreDAG& current_dag, const std::vector<int>& front_layer, const std::vector<std::vector<int>>& lookahead_layers, SabreRunContext& ctx);

    // --- Metody Pomocnicze do ewaluacji (wydzielone z evaluate_fitness_impl) ---
    void reset_evaluation_state(DERunContext& ctx, std::vector<RoutingRecord>* out_log);
    void init_front_layer(int n_gates);
    void extract_executable_gates(DERunContext& ctx);
    void build_lookahead_layers();

    // --- Metody Pomocnicze DE ---
    template <bool LOGGING>
    double finalize_evaluation(DERunContext& ctx, const double* genes, int n_gates, int p_cnt,
                               int& out_swaps, int& out_moves, double& max_physical_time,
                               int& max_logical_depth, std::vector<RoutingRecord>* out_log);

    std::vector<int> prepare_sabre_seed();
    void initialize_de_population(int pop_size, int n_genes, int n_gates, int blocks, const std::vector<int>& logical_exec_order, std::vector<double>& pop, std::vector<double>& fitness, boost::random::xoshiro256pp& rng, bool can_lazy);

    template <bool LOGGING>
    void de_force_down(int& p, int gate_id, DERunContext& ctx, const double* genes, int n_gates, std::vector<RoutingRecord>* out_log, bool can_lazy);

    template <bool LOGGING>
    void de_execute_gates(const HybridSabreDAG& current_dag, std::vector<int>& exec_gates, std::vector<int>& front_layer, std::vector<int>& deps, DERunContext& ctx, const double* genes, int n_gates, std::vector<RoutingRecord>* out_log, bool can_lazy);

    template <bool LOGGING>
    bool de_route_step(const HybridSabreDAG& current_dag,
    const std::vector<int>& front_layer,
    std::vector<double>& decay,
    DERunContext& ctx,
    const double* genes, int n_gates,
    std::vector<RoutingRecord>* out_log);

    template <bool LOGGING>
    double evaluate_fitness_impl(const double* genes, int& out_swaps, int& out_moves, double& max_physical_time, int& max_logical_depth, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard);
    double evaluate_fitness(const double* genes, int& out_swaps, int& out_moves, double& max_physical_time, int& max_logical_depth, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard);

    void run_evolution_loop_windowed(std::vector<double>& pop, std::vector<double>& fitness, int pop_size, int generations, int n_genes, int n_gates, double f, double cr, boost::random::xoshiro256pp& rng, int w_start, int w_end);
    void generate_trial_vector_windowed(const std::vector<double>& pop, std::vector<double>& trial, int current_idx, int a, int b, int c, int n_genes, int n_gates, double f, double cr, boost::random::xoshiro256pp& rng, int w_start, int w_end);

    // --- Narzędzia ---
    void print_statistics(double best_score, double physical_time, int depth, int sw, int mv, double exec_time);
    void print_optimized_circuit(const std::vector<RoutingRecord>& log) const;
    void export_moves_to_file(const std::vector<RoutingRecord>& moves, const std::string& filename);
    bool validate_logical_circuit(const std::vector<RoutingRecord>& log) const;

public:
    SabreDEHybridRouter(const BaseTopology& topology, const std::vector<Gate>& circuit, int qubits_cnt, std::string name, int max_look_ahead_depth = 1);

    void optimize(int pop_size, int generations, bool allow_lazy_offboard = false, double f = 0.8, double cr = 0.9);
    void optimize_windowed(int pop_size, int generations_per_window, int window_size = 100, double accept_ratio = 0.5, bool allow_lazy_offboard = false, double f = 0.8, double cr = 0.9);
};

#endif //MAGISTERKA_SABREDEHYBRIDROUTER_H