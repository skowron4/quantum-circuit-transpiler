//
// Created by skowr on 19.05.2026.
//

#ifndef MAGISTERKA_GREEDYASAP_H
#define MAGISTERKA_GREEDYASAP_H
#include <chrono>
#include <iostream>
#include <queue>

#include "architecture/BaseTopology.h"
#include "differentialEvolution/differential_evolution.h"


// Własny super szybki max omijający koszty std::initializer_list
inline double fast_max(double a, double b) {
    return (a > b) ? a : b;
}

class GreedyASAP {
private:
    static constexpr int RES_OFFSET = 1000;

    const BaseTopology& topo;
    std::vector<GateNode> grouped_nodes;
    std::vector<Gate> original_circuit;
    int num_qubits;
    std::string name;

    void preprocess_circuit(const std::vector<Gate>& circuit) {
        std::vector<int> last_gate_on_qubit(num_qubits, -1);
        std::vector<int> pending_1q(num_qubits, 0);
        int group_id = 0;

        for (const auto& g : circuit) {
            if (g.type == "1Q") {
                pending_1q[g.target_q1]++;
            } else {
                GateNode node;
                node.id = group_id++;
                node.lq1 = g.target_q1;
                node.lq2 = g.target_q2;
                node.q1_1q_count = pending_1q[g.target_q1];
                node.q2_1q_count = pending_1q[g.target_q2];

                pending_1q[g.target_q1] = 0;
                pending_1q[g.target_q2] = 0;

                if (last_gate_on_qubit[node.lq1] != -1) {
                    grouped_nodes[last_gate_on_qubit[node.lq1]].next_gates.push_back(node.id);
                    node.deps_count++;
                }
                if (last_gate_on_qubit[node.lq2] != -1 && last_gate_on_qubit[node.lq2] != last_gate_on_qubit[node.lq1]) {
                    grouped_nodes[last_gate_on_qubit[node.lq2]].next_gates.push_back(node.id);
                    node.deps_count++;
                }

                grouped_nodes.push_back(node);
                last_gate_on_qubit[node.lq1] = node.id;
                last_gate_on_qubit[node.lq2] = node.id;
            }
        }

        for (int i = 0; i < num_qubits; ++i) {
            if (pending_1q[i] > 0) {
                GateNode node;
                node.id = group_id++;
                node.lq1 = i;
                node.lq2 = -1;
                node.q1_1q_count = pending_1q[i];
                node.q2_1q_count = 0;

                if (last_gate_on_qubit[i] != -1) {
                    grouped_nodes[last_gate_on_qubit[i]].next_gates.push_back(node.id);
                    node.deps_count++;
                }
                grouped_nodes.push_back(node);
                last_gate_on_qubit[i] = node.id;
                pending_1q[i] = 0;
            }
        }
    }

public:
    GreedyASAP(const BaseTopology& topology, const std::vector<Gate>& circuit, int qubits_cnt, std::string name)
        : topo(topology), original_circuit(circuit), num_qubits(qubits_cnt), name(name) {
        preprocess_circuit(circuit);
    }

    void run_benchmark() {
        std::cout << "Uruchamiam benchmark Greedy ASAP... " << std::endl;
        auto start_time = std::chrono::high_resolution_clock::now();

        int out_swaps = 0, out_moves = 0, max_logical_depth = 0;
        double max_physical_time = 0.0, total_routing_cost = 0.0;
        std::vector<RoutingRecord> log;

        int n_gates = grouped_nodes.size();
        const double MOVE_COST = 1.0;
        const double SWAP_COST = 3.14;

        MachineState state;
        state.phys_qubit_to_logical.assign(topo.qubits.size(), -1);
        state.logical_to_phys_qubit.assign(num_qubits, -1);
        state.q_on_resonator.assign(topo.qubits.size(), -1);
        state.res_occupant.assign(topo.resonators.size(), -1);
        state.q_ready_time.assign(topo.qubits.size(), 0.0);
        state.r_ready_time.assign(topo.resonators.size(), 0.0);
        std::vector<int> q_logical_depth(topo.qubits.size(), 0);

        auto force_qubit_down = [&](int p, int current_gate_id) {
            if (state.q_on_resonator[p] != -1) {
                int r = state.q_on_resonator[p];
                double step_start = fast_max(state.q_ready_time[p], state.r_ready_time[r]);

                int lq = state.phys_qubit_to_logical[p];
                log.push_back({OpType::MOVE, current_gate_id, lq, -1, r + RES_OFFSET, p, -1, step_start});

                out_moves++;
                total_routing_cost += MOVE_COST;
                state.q_ready_time[p] = step_start + MOVE_COST;
                state.r_ready_time[r] = step_start + MOVE_COST;
                q_logical_depth[p] += 1;
                state.q_on_resonator[p] = -1;
                state.res_occupant[r] = -1;
            }
        };

        for (int i = 0; i < num_qubits; ++i) {
            state.logical_to_phys_qubit[i] = i;
            state.phys_qubit_to_logical[i] = i;
        }

        auto cmp = [](int left, int right) { return left > right; };
        std::priority_queue<int, std::vector<int>, decltype(cmp)> ready_queue(cmp);

        std::vector<int> deps(n_gates);
        for (int i = 0; i < n_gates; ++i) {
            deps[i] = grouped_nodes[i].deps_count;
            if (deps[i] == 0) ready_queue.push(i);
        }

        while (!ready_queue.empty()) {
            int gate_id = ready_queue.top();
            ready_queue.pop();
            const GateNode& node = grouped_nodes[gate_id];

            int p1 = state.logical_to_phys_qubit[node.lq1];

            if (node.lq2 == -1) {
                if (node.q1_1q_count > 0) force_qubit_down(p1, node.id);
                for (int i = 0; i < node.q1_1q_count; ++i) {
                    log.push_back({OpType::GATE_1Q, node.id, node.lq1, -1, p1, -1, -1, state.q_ready_time[p1] + i});
                }
                state.q_ready_time[p1] += node.q1_1q_count;
                q_logical_depth[p1] += node.q1_1q_count;

                for (int next_id : node.next_gates) {
                    if (--deps[next_id] == 0) ready_queue.push(next_id);
                }
                continue;
            }

            int p2 = state.logical_to_phys_qubit[node.lq2];
            if (node.q1_1q_count > 0) force_qubit_down(p1, node.id);
            if (node.lq2 != -1 && node.q2_1q_count > 0) force_qubit_down(p2, node.id);

            for (int i = 0; i < node.q1_1q_count; ++i) log.push_back({OpType::GATE_1Q, node.id, node.lq1, -1, p1, -1, -1, state.q_ready_time[p1] + i});
            for (int i = 0; i < node.q2_1q_count; ++i) log.push_back({OpType::GATE_1Q, node.id, node.lq2, -1, p2, -1, -1, state.q_ready_time[p2] + i});

            state.q_ready_time[p1] += node.q1_1q_count;
            state.q_ready_time[p2] += node.q2_1q_count;
            q_logical_depth[p1] += node.q1_1q_count;
            q_logical_depth[p2] += node.q2_1q_count;

            std::vector<int> shared_res;
            const auto& p2_resonators = topo.qubits[p2].connected_resonators;
            for (int r : topo.qubits[p1].connected_resonators) {
                // ZOPTYMALIZOWANE: Błyskawiczny std::ranges::find zamiast powolnego drzewiastego set::count()
                if (std::ranges::find(p2_resonators, r) != p2_resonators.end()) shared_res.push_back(r);
            }

            double start_time = fast_max(state.q_ready_time[p1], state.q_ready_time[p2]);
            int gate_duration = 1;

            bool q1_goes_in = (state.q_ready_time[p1] <= state.q_ready_time[p2]);

            int p_in = q1_goes_in ? p1 : p2;
            int p_stay = q1_goes_in ? p2 : p1;
            int lq_in = q1_goes_in ? node.lq1 : node.lq2;
            int lq_stay = q1_goes_in ? node.lq2 : node.lq1;

            if (!shared_res.empty()) {
                int chosen_res = shared_res[0];
                double earliest_res_time = state.r_ready_time[chosen_res];
                for (size_t i = 1; i < shared_res.size(); ++i) {
                    if (state.r_ready_time[shared_res[i]] < earliest_res_time) {
                        earliest_res_time = state.r_ready_time[shared_res[i]];
                        chosen_res = shared_res[i];
                    }
                }

                force_qubit_down(p_stay, node.id);
                if (state.res_occupant[chosen_res] != -1 && state.res_occupant[chosen_res] != p_in) force_qubit_down(state.res_occupant[chosen_res], node.id);

                if (state.q_on_resonator[p_in] != chosen_res) {
                    force_qubit_down(p_in, node.id);
                    // ZOPTYMALIZOWANE: fast_max niszczy narzut czasowy std::max z initializer_list
                    start_time = fast_max(state.q_ready_time[p_in], fast_max(state.q_ready_time[p_stay], state.r_ready_time[chosen_res]));

                    log.push_back({OpType::MOVE, node.id, lq_in, -1, p_in, chosen_res + RES_OFFSET, -1, start_time});
                    out_moves++;
                    total_routing_cost += MOVE_COST;
                    start_time += MOVE_COST;

                    state.q_on_resonator[p_in] = chosen_res;
                    state.res_occupant[chosen_res] = p_in;
                    q_logical_depth[p_in] += 1;
                } else {
                    start_time = fast_max(state.q_ready_time[p_in], fast_max(state.q_ready_time[p_stay], state.r_ready_time[chosen_res]));
                }

                log.push_back({OpType::GATE_2Q_LOCAL, node.id, lq_in, lq_stay, chosen_res + RES_OFFSET, p_stay, p_in, start_time});
                state.q_ready_time[p_in] = start_time + gate_duration;
                state.q_ready_time[p_stay] = start_time + gate_duration;
                state.r_ready_time[chosen_res] = start_time + gate_duration;

                int gate_layer = std::max(q_logical_depth[p_in], q_logical_depth[p_stay]);
                q_logical_depth[p_in] = gate_layer + 1;
                q_logical_depth[p_stay] = gate_layer + 1;

            } else {
                force_qubit_down(p1, node.id);
                force_qubit_down(p2, node.id);

                // KRYTYCZNA ZMIANA OGRANICZAJĄCA ALOKACJĘ PAMIĘCI (zwrócenie referencji a nie kopii z macierzy)
                const auto& path = topo.get_routing_path(p1, p2);
                int current_p = p1;

                if (path.resonators.empty()) continue;

                for (size_t i = 0; i < path.resonators.size() - 1; ++i) {
                    int r_step = path.resonators[i];
                    int next_p = path.qubits[i+1];

                    force_qubit_down(next_p, node.id);
                    double step_start = fast_max(state.q_ready_time[current_p], fast_max(state.r_ready_time[r_step], state.q_ready_time[next_p]));

                    bool is_move = (state.phys_qubit_to_logical[next_p] == -1);
                    double op_cost = is_move ? (2.0 * MOVE_COST) : SWAP_COST;

                    total_routing_cost += op_cost;
                    int lq_mov = state.phys_qubit_to_logical[current_p];
                    int lq_target = state.phys_qubit_to_logical[next_p];

                    if (is_move) {
                        log.push_back({OpType::MOVE, node.id, lq_mov, -1, current_p, r_step + RES_OFFSET, -1, step_start});
                        log.push_back({OpType::MOVE, node.id, lq_mov, -1, r_step + RES_OFFSET, next_p, -1, step_start + MOVE_COST});
                        out_moves += 2;
                    } else {
                        log.push_back({OpType::SWAP, node.id, lq_mov, lq_target, current_p, next_p, r_step, step_start});
                        out_swaps++;
                    }

                    int current_layer = std::max(q_logical_depth[current_p], q_logical_depth[next_p]);
                    int layers_added = is_move ? 2 : 1;
                    q_logical_depth[current_p] = current_layer + layers_added;
                    q_logical_depth[next_p] = current_layer + layers_added;

                    state.phys_qubit_to_logical[next_p] = lq_mov;
                    state.logical_to_phys_qubit[lq_mov] = next_p;

                    if (!is_move && lq_target != -1) {
                        state.phys_qubit_to_logical[current_p] = lq_target;
                        state.logical_to_phys_qubit[lq_target] = current_p;
                    } else {
                        state.phys_qubit_to_logical[current_p] = -1;
                    }

                    state.q_ready_time[current_p] = step_start + op_cost;
                    state.r_ready_time[r_step] = step_start + op_cost;
                    state.q_ready_time[next_p] = step_start + op_cost;

                    current_p = next_p;
                }

                int final_res = path.resonators.back();

                bool g_q1_goes_in = (state.q_ready_time[current_p] <= state.q_ready_time[p2]);
                int g_p_in = g_q1_goes_in ? current_p : p2;
                int g_p_stay = g_q1_goes_in ? p2 : current_p;
                int g_lq_in_final = state.phys_qubit_to_logical[g_p_in];
                int g_lq_stay_final = state.phys_qubit_to_logical[g_p_stay];

                force_qubit_down(g_p_stay, node.id);
                if (state.res_occupant[final_res] != -1 && state.res_occupant[final_res] != g_p_in) force_qubit_down(state.res_occupant[final_res], node.id);

                double exec_start = fast_max(state.q_ready_time[current_p], fast_max(state.q_ready_time[p2], state.r_ready_time[final_res]));

                log.push_back({OpType::MOVE, node.id, g_lq_in_final, -1, g_p_in, final_res + RES_OFFSET, -1, exec_start});
                out_moves++;
                total_routing_cost += MOVE_COST;
                exec_start += MOVE_COST;

                state.q_on_resonator[g_p_in] = final_res;
                state.res_occupant[final_res] = g_p_in;
                q_logical_depth[g_p_in] += 1;

                log.push_back({OpType::GATE_2Q_GLOBAL, node.id, g_lq_in_final, g_lq_stay_final, final_res + RES_OFFSET, g_p_stay, g_p_in, exec_start});

                state.q_ready_time[current_p] = exec_start + gate_duration;
                state.q_ready_time[p2] = exec_start + gate_duration;
                state.r_ready_time[final_res] = exec_start + gate_duration;

                int global_gate_layer = std::max(q_logical_depth[current_p], q_logical_depth[p2]);
                q_logical_depth[current_p] = global_gate_layer + 1;
                q_logical_depth[p2] = global_gate_layer + 1;
            }

            for (int next_id : node.next_gates) {
                if (--deps[next_id] == 0) ready_queue.push(next_id);
            }
        }

        for (int p = 0; p < topo.qubits.size(); ++p) force_qubit_down(p, -1);

        for (double t : state.q_ready_time) max_physical_time = std::max(max_physical_time, t);
        for (int d : q_logical_depth) max_logical_depth = std::max(max_logical_depth, d);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> execution_time = end_time - start_time;

        std::cout << "\n=== WYNIKI BENCHMARKU (GREEDY ASAP) ===" << std::endl;
        std::cout << "Wykorzystane operacje SWAP: " << out_swaps << std::endl;
        std::cout << "Wykorzystane operacje MOVE: " << out_moves << std::endl;
        std::cout << "Czas fizyczny obwodu (Makespan): " << max_physical_time << " [jednostek]" << std::endl;
        std::cout << "Głębokość warstwowa (Logical Depth): " << max_logical_depth << " [warstw]" << std::endl;
        std::cout << "Koszt routingu (Total Routing Cost): " << total_routing_cost << std::endl;
        std::cout << "Czas wykonania algorytmu: " << execution_time.count() << " sekund" << std::endl;
        std::cout << "=======================================\n" << std::endl;
    }
};


#endif //MAGISTERKA_GREEDYASAP_H