#ifndef MAGISTERKA_DIFFERENTIAL_EVOLUTION_H
#define MAGISTERKA_DIFFERENTIAL_EVOLUTION_H

#include "architecture/IQMConstellation.h"
#include "structures/Gate.h"
#include <vector>
#include "structures/routing_record.hpp"

struct GroupedGate {
    int id;
    int q1, q2;
    int q1_1q_count;
    int q2_1q_count;
    int type;
    std::vector<int> shared_resonators;
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
    std::vector<int> next_gates;
};

class DifferentialEvolution {
private:
    const BaseTopology& topo;
    std::vector<GroupedGate> grouped_circuit;
    std::vector<GateNode> grouped_nodes;
    std::vector<Gate> original_circuit;
    int num_qubits;
    std::string name;

    void preprocess_circuit(const std::vector<Gate>& circuit);
    bool validate_logical_circuit(const std::vector<Gate>& orig_circ, const std::vector<RoutingRecord>& log) const; // Dodano metodę walidacyjną

public:
    DifferentialEvolution(const BaseTopology& topology, const std::vector<Gate>& circuit, int qubits_cnt, std::string name);

    double evaluate_fitness(const std::vector<double>& genes, int& out_swaps, int& out_moves, std::vector<RoutingRecord>* out_log = nullptr);
    void print_optimized_circuit(const std::vector<RoutingRecord>& log, int num_logical_qubits) const;
    void optimize(int pop_size, int generations, int num_logical_qubits = -1);

};

#endif //MAGISTERKA_DIFFERENTIAL_EVOLUTION_H