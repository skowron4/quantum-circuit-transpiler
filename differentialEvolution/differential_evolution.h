#ifndef MAGISTERKA_DIFFERENTIAL_EVOLUTION_H
#define MAGISTERKA_DIFFERENTIAL_EVOLUTION_H

#include "architecture/IQMConstellation.h"
#include "structures/Gate.h"
#include <vector>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/xoshiro.hpp>

#include "structures/routing_record.hpp"

struct GroupedGate {
    int id;
    int q1, q2;
    int q1_1q_count;
    int q2_1q_count;
    int type;
    absl::InlinedVector<int, 10> shared_resonators;
};

struct MachineState {
    std::vector<int> phys_qubit_to_logical;
    std::vector<int> logical_to_phys_qubit;

    // --- NOWE ZMIENNE DO LAZY OFF-BOARDING ---
    std::vector<int> q_on_resonator; // [fizyczny_q] -> ID rezonatora (-1 jeśli jest "na ziemi")
    std::vector<int> res_occupant;   // [ID_rezonatora] -> fizyczny_q, który aktualnie na nim wisi (-1 jeśli pusty)
    // -----------------------------------------

    std::vector<double> q_ready_time;
    std::vector<double> r_ready_time;
};

struct GateNode {
    int id;
    int lq1, lq2;
    int q1_1q_count, q2_1q_count;

    int deps_count = 0;
    absl::InlinedVector<int, 8> next_gates;
};

class DifferentialEvolution {
private:
    static constexpr double MOVE_COST = 1.0;
    static constexpr double SWAP_COST = 3.14;

    const BaseTopology& topo;
    std::vector<GroupedGate> grouped_circuit;
    std::vector<GateNode> grouped_nodes;
    std::vector<Gate> original_circuit;
    int num_qubits;
    std::string name;

    MachineState state_workspace;
    std::vector<int> q_logical_depth_workspace;
    std::vector<int> deps_workspace;
    std::vector<int> shared_res_workspace;
    std::vector<int> ready_queue_workspace;

    void prepare_workspaces(int n_gates, std::vector<RoutingRecord>* out_log);

    template <bool LOGGING>
    void force_qubit_down(int& p, int current_gate_id, const std::vector<double>& genes, int n_gates, int& out_moves, double& total_routing_cost, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard);

    template <bool LOGGING>
    void apply_1q_gates(int gate_id, int lq, int& p, int q_1q_count, const std::vector<double>& genes, int n_gates, int& out_moves, double& total_routing_cost, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard);

    template <bool LOGGING>
    void execute_local_2q_gate(const GateNode& node, int p1, int p2, int gate_id, const std::vector<double>& genes, int n_gates, int& out_moves, double& total_routing_cost, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard);

    template <bool LOGGING>
    void execute_routed_2q_gate(const GateNode& node, int p1, int p2, int gate_id, const std::vector<double>& genes, int n_gates, int& out_swaps, int& out_moves, double& total_routing_cost, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard);

    template <bool LOGGING>
    double finalize_and_calculate_fitness(const std::vector<double>& genes, int n_gates, int& out_moves, double& total_routing_cost, double& max_physical_time, int& max_logical_depth, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard);

    template <bool LOGGING>
    void process_routing_path(const GateNode& node, const RoutingPath& path, int& current_p, const std::vector<double>& genes, int n_gates, int& out_swaps, int& out_moves, double& total_routing_cost, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard);

    void preprocess_circuit(const std::vector<Gate>& circuit);

    bool validate_logical_circuit(const std::vector<Gate>& orig_circ, const std::vector<RoutingRecord>& log) const; // Dodano metodę walidacyjną

    void initialize_population(std::vector<double>& pop, std::vector<double>& fitness, int pop_size, int n_genes, boost::random::xoshiro256pp& rng);

    void run_evolution_loop(std::vector<double>& pop, std::vector<double>& fitness, int pop_size, int generations, int n_genes, double f, double cr, boost::random::xoshiro256pp& rng);

    void select_distinct_parents(int current_idx, int pop_size, boost::random::xoshiro256pp& rng, boost::random::uniform_int_distribution<int>& dist_int, int& a, int& b, int& c);

    void generate_trial_vector(const std::vector<double>& pop, std::vector<double>& trial, int current_idx, int a, int b, int c, int n_genes, double f, double cr, boost::random::xoshiro256pp& rng);

    void process_and_export_best_solution(const std::vector<double>& pop, const std::vector<double>& fitness, int n_genes, int num_logical_qubits, double exec_time_sec);

    void print_statistics(double best_score, int swaps, int moves, double physical_time, int depth, double exec_time);

    void export_moves_to_file(const std::vector<RoutingRecord>& moves, const std::string& filename);

public:
    DifferentialEvolution(const BaseTopology& topology, const std::vector<Gate>& circuit, int qubits_cnt, std::string name);

    template <bool LOGGING>
    double evaluate_fitness_impl(const std::vector<double>& genes, int& out_swaps, int& out_moves, double &max_physical_time, int &max_logical_depth, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard);

    double evaluate_fitness(const std::vector<double>& genes, int& out_swaps, int& out_moves, double &max_physical_time, int &max_logical_depth, std::vector<RoutingRecord>* out_log);

    void print_optimized_circuit(const std::vector<RoutingRecord>& log, int num_logical_qubits) const;

    void optimize(int pop_size, int generations, int num_logical_qubits = -1. , double f = 0.8, double cr = 0.9);

};

#endif //MAGISTERKA_DIFFERENTIAL_EVOLUTION_H