//
// Created by skowr on 23.05.2026.
//

#ifndef MAGISTERKA_SABREDEHYBRIDROUTER_H
#define MAGISTERKA_SABREDEHYBRIDROUTER_H

#pragma once
#include <chrono>

#include "architecture/BaseTopology.h"
#include "structures/Gate.h"
#include "structures/routing_record.hpp"
#include <vector>
#include <string>
#include <limits>
#include <cstdint>
#include <boost/random/xoshiro.hpp>

struct OptimizationResult {
    double score;
    double makespan;
    int depth;
    int swaps;
    int moves;
};

struct SwapEdge { int neighbor; int resonator; };

struct QubitInteraction {
    int target;
    double weight;
};

struct DEStateSnapshot {
    bool is_valid = false;
    int target_gates_routed = 0;
    int actual_gates_routed = 0;

    // Dane z ctx
    std::vector<int> l2p, p2l, q_on_res, res_occ, q_depth;
    std::vector<double> q_time, r_time;
    int out_swaps = 0, out_moves = 0;
    double total_routing_cost = 0.0;

    // Dane z workspace
    std::vector<int> deps, front_layer;
    std::vector<double> decay;
};

// --- Struktury wewnątrz hybrydy ---
struct HybridGateNode {
    int id;
    int lq1, lq2;
    int q1_1q_count, q2_1q_count;
    int deps_count = 0;
    absl::InlinedVector<int, 4> q1_1q_ids;
    absl::InlinedVector<int, 4> q2_1q_ids;
    int orig_2q_id;
    // std::vector<int> next_gates;

    // 2. MAGIA: Podmieniamy std::vector.
    // Liczba '2' oznacza, że jeśli bramka ma 1 lub 2 następników (co w obwodach kwantowych jest regułą),
    // pamięć NIE zostanie zaalokowana na stercie, tylko "wklejona" w samą strukturę.
    absl::InlinedVector<int, 2> next_gates;

    HybridGateNode(int id, int lq1, int lq2, int q1_cnt, int q2_cnt,
                   absl::InlinedVector<int, 4> q1_ids,
                   absl::InlinedVector<int, 4> q2_ids,
                   int orig_2q)
        // POPRAWKA 2: Kolejność na liście inicjalizacyjnej idealnie zgadza się z deklaracjami
        : id(id),
          lq1(lq1),
          lq2(lq2),
          q1_1q_count(q1_cnt),
          q2_1q_count(q2_cnt),
          // deps_count pomijamy, bo ma domyślną wartość 0 w samej deklaracji
          q1_1q_ids(std::move(q1_ids)),
          q2_1q_ids(std::move(q2_ids)),
          orig_2q_id(orig_2q) {}
};

struct HybridSabreDAG {
    std::vector<HybridGateNode> nodes;
    std::vector<int> initial_deps;
};

// Zoptymalizowany HybridWorkspace (Usunięto SwapCand, dodano flagi O(1) i pre-alokowane pary)
struct HybridWorkspace {
    std::vector<int> l2p, p2l, q_on_res, res_occ, q_depth, deps, ready_queue, shared_res;
    std::vector<double> q_time, r_time, phys_gene_weights;;

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
    static constexpr int MAX_INTER_PER_Q = 64; // Zapas, którego praktycznie nigdy nie przekroczymy
    std::vector<int> inter_targets;
    std::vector<double> inter_weights;
    std::vector<int> interaction_counts;
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

struct SabreSeedProfile {
    std::vector<int> logical_exec_order;
    std::vector<int> resonator_choice_index;
    std::vector<int> resonator_choice_count;
    std::vector<int> q1_goes_in;
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
    static constexpr double DECAY_RESET = 1.0;
    static constexpr double DECAY_INCREMENT = 0.2;

    const BaseTopology& topo;
    std::vector<Gate> original_circuit;
    int num_qubits;
    std::string name;
    int max_look_ahead_depth;
    int p_cnt_cache;

    bool slow;
    uint64_t random_seed = 1337;
    int sabre_seed_trials = 1;
    mutable boost::random::xoshiro256pp sabre_rng;

    std::vector<SwapEdge> flat_precomputed_swaps;
    std::vector<int> precomputed_swaps_offsets;

    mutable std::vector<int> dist_cache;
    HybridSabreDAG dag;
    HybridSabreDAG reversed_dag;
    HybridWorkspace workspace;
    std::vector<int> saved_initial_l2p;
    SabreSeedProfile sabre_seed_profile;

    // --- Metody Bazowe ---
    int get_distance_old(int p1, int p2) const;
    int get_distance(int p1, int p2) const;
    HybridSabreDAG build_dag(const std::vector<Gate>& circuit) const;
    void add_2q_node_to_dag(HybridSabreDAG& local_dag, const Gate& g, std::vector<absl::InlinedVector<int, 4>>& pending_1q, std::vector<int>& last_gate, int& group_id) const;
    void flush_trailing_1q_nodes(HybridSabreDAG& local_dag, std::vector<absl::InlinedVector<int, 4>>& pending_1q, std::vector<int>& last_gate, int& group_id) const;
    void compute_initial_deps(HybridSabreDAG& local_dag) const;

    // --- Metody Pomocnicze SABRE ---
    std::vector<int> run_sabre_pass(const HybridSabreDAG& current_dag, const std::vector<int>& initial_l2p, bool is_final, SabreSeedProfile* out_seed_profile, double* out_routing_cost = nullptr, bool print_details = true);
    void sabre_force_qubit_down(int p, int current_gate_id, SabreRunContext& ctx);
    void sabre_execute_gates(const HybridSabreDAG& current_dag, std::vector<int>& executable_gates, std::vector<int>& front_layer, std::vector<int>& deps, SabreRunContext& ctx, SabreSeedProfile* out_seed_profile);
    bool sabre_route_step(const HybridSabreDAG& current_dag, const std::vector<int>& front_layer, const std::vector<std::vector<int>>& lookahead_layers, SabreRunContext& ctx);
    inline void sabre_execute_1q_gate(const HybridGateNode& node, int p1, SabreRunContext& ctx);
    inline void sabre_execute_2q_gate(const HybridGateNode& node, int p1, int p2, SabreRunContext& ctx, SabreSeedProfile* out_seed_profile);
    inline void extract_sabre_executables(const HybridSabreDAG& current_dag, std::vector<int>& front_layer, std::vector<int>& executable_gates, SabreRunContext& ctx);
    inline void build_sabre_lookahead(const HybridSabreDAG& current_dag, const std::vector<int>& front_layer, std::vector<std::vector<int>>& lookahead_layers, std::vector<bool>& in_extended);
    std::vector<int> create_random_initial_mapping();

    // --- Metody Pomocnicze do ewaluacji (wydzielone z evaluate_fitness_impl) ---
    void reset_evaluation_state(DERunContext& ctx, std::vector<RoutingRecord>* out_log);
    void init_front_layer(int n_gates);
    void extract_executable_gates(DERunContext& ctx);
    void build_lookahead_layers();
    // Nowe metody pomocnicze do resetowania stanu
    inline void reset_mappings(DERunContext& ctx);
    inline void reset_hardware_states(DERunContext& ctx);
    // Nowe metody pomocnicze do inicjalizacji front_layer
    inline void load_initial_dependencies();
    inline void populate_initial_front_layer(int n_gates);
    // Nowa metoda pomocnicza do sprawdzania gotowości bramki
    inline bool check_execution_ready(const HybridGateNode& node, const DERunContext& ctx) const;
    // Nowe metody pomocnicze do BFS (Lookahead)
    inline void setup_lookahead_bfs();
    inline void process_lookahead_layer();

    // --- Metody Pomocnicze DE ---
    template <bool LOGGING>
    bool de_naive_route_step(
        const HybridSabreDAG& current_dag,
        const std::vector<int>& front_layer,
        std::vector<double>& decay,
        DERunContext& ctx,
        const double* genes, int n_gates,
        std::vector<RoutingRecord>* out_log);

    template <bool LOGGING>
    double finalize_evaluation(DERunContext& ctx, const double* genes, int n_gates, int p_cnt,
                               int& out_swaps, int& out_moves, double& max_physical_time,
                               int& max_logical_depth, std::vector<RoutingRecord>* out_log);

    std::vector<int> prepare_sabre_seed(bool can_lazy, int n_genes, int blocks);
    std::vector<int> create_time_decay_mapping() const;
    std::vector<int> create_identity_mapping() const;
    std::vector<int> create_degree_based_mapping() const;
    void initialize_de_population(int pop_size, int n_genes, int n_gates, int blocks, const SabreSeedProfile& seed_profile, std::vector<double>& pop, std::vector<double>& fitness, boost::random::xoshiro256pp& rng, bool can_lazy);
    double calculate_fitness_score(double routing_cost, int swaps, int moves, double max_time, int max_depth);
    template <bool LOGGING>
    void force_remaining_qubits_down(DERunContext& ctx, const double* genes, int n_gates, int p_cnt, std::vector<RoutingRecord>* out_log);
    void extract_max_metrics(const DERunContext& ctx, double& max_physical_time, int& max_logical_depth);

    template <bool LOGGING, bool CAN_LAZY>
    void de_force_down(int& p, int gate_id, DERunContext& ctx, const double* genes, int n_gates, std::vector<RoutingRecord>* out_log);

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
    bool de_route_step_better(
        const HybridSabreDAG& current_dag,
        const std::vector<int>& front_layer,
        std::vector<double>& decay,
        DERunContext& ctx,
        const double* genes, int n_gates,
        std::vector<RoutingRecord>* out_log);

    inline int de_select_lazy_target(int p, int r, int gate_id, const DERunContext& ctx, const double* genes, int n_gates) const;

    template <bool LOGGING>
    void de_execute_1q_logic(const HybridGateNode& node, int p1, DERunContext& ctx, const double* genes, int n_gates, std::vector<RoutingRecord>* out_log, bool can_lazy);
    template <bool LOGGING>
    bool de_execute_2q_logic(const HybridGateNode& node, int p1, int p2, DERunContext& ctx, const double* genes, int n_gates, std::vector<RoutingRecord>* out_log, bool can_lazy);
    inline void de_update_dependencies(const HybridGateNode& node, std::vector<int>& front_layer, std::vector<int>& deps);
    void de_collect_active_qubits(const HybridSabreDAG& current_dag, const std::vector<int>& front_layer, DERunContext& ctx);
    template <bool LOGGING>
    void de_execute_route_candidate(int p1, int p2, int r_step, DERunContext& ctx, std::vector<double>& decay, const double* genes, int n_gates, std::vector<RoutingRecord>* out_log, bool can_lazy);

    template <bool LOGGING>
    double evaluate_fitness_impl(const double* genes, int& out_swaps, int& out_moves, double& max_physical_time, int& max_logical_depth, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard, DEStateSnapshot* load_snapshot = nullptr, DEStateSnapshot* save_snapshot = nullptr);
    double evaluate_fitness(const double* genes, int& out_swaps, int& out_moves, double& max_physical_time, int& max_logical_depth, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard);

    void run_evolution_loop_windowed(std::vector<double>& pop, std::vector<double>& fitness, int pop_size, int generations, int n_genes, int n_gates, double f, double cr, boost::random::xoshiro256pp& rng, int w_start, int w_end, DEStateSnapshot* snapshot, const std::vector<int>& gate_time_map);
    void generate_trial_vector_windowed(const std::vector<double>& pop, std::vector<double>& trial, int current_idx, int a, int b, int c, int n_genes, int n_gates, double f, double cr, boost::random::xoshiro256pp& rng, int w_start, int w_end, const std::vector<int>& gate_time_map, DEStateSnapshot* snapshot);

    inline void load_state_from_snapshot(DERunContext& ctx, const DEStateSnapshot* load_snapshot, int& gates_routed, std::vector<RoutingRecord>* out_log);
    inline void save_state_to_snapshot(const DERunContext& ctx, DEStateSnapshot* save_snapshot, int gates_routed);
    inline void select_distinct_parents(int pop_size, int current_idx, int& a, int& b, int& c, boost::random::xoshiro256pp& rng);

    // --- Narzędzia ---
    void print_statistics(double best_score, double physical_time, int depth, int sw, int mv, double exec_time);
    void print_optimized_circuit(const std::vector<RoutingRecord>& log) const;
    void export_moves_to_file(const std::vector<RoutingRecord>& moves, const std::string& filename);
    bool validate_logical_circuit(const std::vector<RoutingRecord>& log) const;

    // --- Metody pomocnicze do konstruktora ---
    void init_distance_cache();
    void precompute_physical_swaps();
    void allocate_workspace(int g_cnt);
    // --- Metody pomocnicze do głównej pętli DE ---
    void run_standard_evolution_loop(std::vector<double>& pop, std::vector<double>& fitness, int pop_size, int generations, int n_genes, int n_gates, double f, double cr, boost::random::xoshiro256pp& rng, bool can_lazy);
    // --- Metody pomocnicze do Windowed DE ---
    void propagate_window_success(std::vector<double>& pop, std::vector<double>& fitness, int pop_size, int n_genes, int n_gates, int blocks, int best_idx, int accept_end, const std::vector<int>& gate_time_map, DEStateSnapshot& current_snapshot, bool can_lazy);
    // --- Finalizacja i logowanie ---
    OptimizationResult finalize_optimization_result(const std::vector<double>& pop, int best_idx, int n_genes, bool can_lazy, std::chrono::time_point<std::chrono::high_resolution_clock> start_time);

public:
    SabreDEHybridRouter(const BaseTopology& topology, const std::vector<Gate>& circuit, int qubits_cnt, std::string name, bool slow = false, int max_look_ahead_depth = 1);

    void set_random_seed(uint64_t seed);
    void set_sabre_seed_trials(int trials);
    void optimize(int pop_size, int generations, bool allow_lazy_offboard = false, double f = 0.8, double cr = 0.9);
    OptimizationResult optimize_windowed(int pop_size, int generations_per_window, int window_size = 100, double accept_ratio = 0.5, bool allow_lazy_offboard = false, double f = 0.9, double cr = 0.9);
    OptimizationResult optimize_windowed_with_feedback(
    int pop_size, int generations_per_window, int window_size, double accept_ratio,
    int feedback_passes, bool allow_lazy_offboard, double f, double cr);
    void optimize_with_feedback(int pop_size, int generations_per_pass, int feedback_passes, bool allow_lazy_offboard, double f, double cr);
    OptimizationResult evaluate_qiskit_routed_circuit(const std::vector<Gate>& qiskit_circuit);
};

#endif //MAGISTERKA_SABREDEHYBRIDROUTER_H