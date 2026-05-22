#ifndef MAGISTERKA_SABREROUTER_H
#define MAGISTERKA_SABREROUTER_H

#include <chrono>
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <string>

#include "architecture/BaseTopology.h"
#include "differentialEvolution/differential_evolution.h"

inline double fast_max_sabre(double a, double b) {
    return (a > b) ? a : b;
}

struct SabreDAG {
    std::vector<GateNode> nodes;
    std::vector<int> initial_deps;
};

class SabreRouter {
private:
    static constexpr int RES_OFFSET = 1000;
    static constexpr double EXTENDED_LAYER_WEIGHT = 0.5;

    const BaseTopology& topo;
    std::vector<Gate> original_circuit;
    int num_qubits;
    std::string name;

    mutable std::vector<std::vector<int>> dist_cache;
    SabreDAG dag; // Zbudowane RAZ w konstruktorze
    SabreDAG reversed_dag; // Zbudowane RAZ dla backward passu

    int get_distance(int p1, int p2) const {
        if (p1 == p2) return 0;
        if (dist_cache[p1][p2] != -1) return dist_cache[p1][p2];

        int dist = topo.get_routing_path(p1, p2).resonators.size();
        dist_cache[p1][p2] = dist;
        dist_cache[p2][p1] = dist;
        return dist;
    }

    SabreDAG build_dag(const std::vector<Gate>& circuit) const {
        SabreDAG local_dag;
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
                    local_dag.nodes[last_gate_on_qubit[node.lq1]].next_gates.push_back(node.id);
                    node.deps_count++;
                }
                if (last_gate_on_qubit[node.lq2] != -1 && last_gate_on_qubit[node.lq2] != last_gate_on_qubit[node.lq1]) {
                    local_dag.nodes[last_gate_on_qubit[node.lq2]].next_gates.push_back(node.id);
                    node.deps_count++;
                }

                local_dag.nodes.push_back(node);
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
                    local_dag.nodes[last_gate_on_qubit[i]].next_gates.push_back(node.id);
                    node.deps_count++;
                }
                local_dag.nodes.push_back(node);
                last_gate_on_qubit[i] = node.id;
                pending_1q[i] = 0;
            }
        }

        local_dag.initial_deps.resize(local_dag.nodes.size());
        for (size_t i = 0; i < local_dag.nodes.size(); ++i) {
            local_dag.initial_deps[i] = local_dag.nodes[i].deps_count;
        }

        return local_dag;
    }

    struct SwapCandidate {
        int p1, p2, r_shared;
        double score;
    };

    std::vector<int> run_pass(const SabreDAG& current_dag,
                              const std::vector<int>& initial_l2p,
                              bool is_final) {

        int n_gates = current_dag.nodes.size();
        const double MOVE_COST = 1.0;
        const double SWAP_COST = 3.14;

        MachineState state;
        state.logical_to_phys_qubit = initial_l2p;
        state.phys_qubit_to_logical.assign(topo.qubits.size(), -1);

        for (int i = 0; i < num_qubits; ++i) {
            if (initial_l2p[i] != -1) {
                state.phys_qubit_to_logical[initial_l2p[i]] = i;
            }
        }

        state.q_on_resonator.assign(topo.qubits.size(), -1);
        state.res_occupant.assign(topo.resonators.size(), -1);
        state.q_ready_time.assign(topo.qubits.size(), 0.0);
        state.r_ready_time.assign(topo.resonators.size(), 0.0);

        std::vector<int> q_logical_depth(topo.qubits.size(), 0);
        std::vector<RoutingRecord> log;
        if (is_final) log.reserve(n_gates * 3);

        int out_swaps = 0, out_moves = 0;
        double total_routing_cost = 0.0;

        auto force_qubit_down = [&](int p, int current_gate_id) {
            if (state.q_on_resonator[p] != -1) {
                int r = state.q_on_resonator[p];
                double step_start = fast_max_sabre(state.q_ready_time[p], state.r_ready_time[r]);
                int lq = state.phys_qubit_to_logical[p];

                if (is_final) log.push_back({OpType::MOVE, current_gate_id, lq, -1, r + RES_OFFSET, p, -1, step_start});

                out_moves++;
                total_routing_cost += MOVE_COST;
                state.q_ready_time[p] = step_start + MOVE_COST;
                state.r_ready_time[r] = step_start + MOVE_COST;
                q_logical_depth[p] += 1;
                state.q_on_resonator[p] = -1;
                state.res_occupant[r] = -1;
            }
        };

        std::vector<int> front_layer;
        front_layer.reserve(num_qubits);
        std::vector<int> deps = current_dag.initial_deps;
        for (int i = 0; i < n_gates; ++i) {
            if (deps[i] == 0) front_layer.push_back(i);
        }

        std::vector<double> decay(topo.qubits.size(), 1.0);
        const double DECAY_RESET = 1.0;
        const double DECAY_INCREMENT = 0.2;
        int gates_routed = 0;

        std::vector<int> executable_gates;
        executable_gates.reserve(num_qubits);

        std::vector<SwapCandidate> candidates;
        candidates.reserve(topo.qubits.size() * 4);

        std::vector<bool> in_active(topo.qubits.size(), false);
        std::vector<int> active_phys_qubits;
        active_phys_qubits.reserve(num_qubits);

        std::vector<bool> in_extended(n_gates, false);
        std::vector<int> extended_layer;
        extended_layer.reserve(num_qubits * 4);

        // DEBUG: Zabezpieczenie przed nieskończoną pętlą
        int debug_iterations = 0;
        const int MAX_ITERATIONS = n_gates * 100;

        while (gates_routed < n_gates) {

            if (++debug_iterations > MAX_ITERATIONS) {
                std::cerr << "\n[CRITICAL ERROR] SABRE ugrzązł w nieskończonej pętli!\n";
                std::cerr << "Bramek ułożonych: " << gates_routed << " / " << n_gates << "\n";
                std::cerr << "Front layer size: " << front_layer.size() << "\n";
                std::cerr << "Candidates generated in last step: " << candidates.size() << "\n";
                break; // Uratuj program przed całkowitym zawieszeniem
            }

            executable_gates.clear();

            for (size_t i = 0; i < front_layer.size(); ) {
                int gate_id = front_layer[i];
                const GateNode& node = current_dag.nodes[gate_id];
                bool is_executable = false;

                if (node.lq2 == -1) {
                    is_executable = true;
                } else {
                    int p1 = state.logical_to_phys_qubit[node.lq1];
                    int p2 = state.logical_to_phys_qubit[node.lq2];
                    if (get_distance(p1, p2) <= 1) is_executable = true;
                }

                if (is_executable) {
                    executable_gates.push_back(gate_id);
                    front_layer[i] = front_layer.back();
                    front_layer.pop_back();
                } else {
                    ++i;
                }
            }

            if (!executable_gates.empty()) {
                for (int gate_id : executable_gates) {
                    const GateNode& node = current_dag.nodes[gate_id];
                    int p1 = state.logical_to_phys_qubit[node.lq1];

                    if (node.lq2 == -1) {
                        if (node.q1_1q_count > 0) force_qubit_down(p1, node.id);
                        if (is_final) {
                            for (int i = 0; i < node.q1_1q_count; ++i) {
                                log.push_back({OpType::GATE_1Q, node.id, node.lq1, -1, p1, -1, -1, state.q_ready_time[p1] + i});
                            }
                        }
                        state.q_ready_time[p1] += node.q1_1q_count;
                        q_logical_depth[p1] += node.q1_1q_count;
                    } else {
                        int p2 = state.logical_to_phys_qubit[node.lq2];
                        if (node.q1_1q_count > 0) force_qubit_down(p1, node.id);
                        if (node.q2_1q_count > 0) force_qubit_down(p2, node.id);

                        if (is_final) {
                            for (int i = 0; i < node.q1_1q_count; ++i) log.push_back({OpType::GATE_1Q, node.id, node.lq1, -1, p1, -1, -1, state.q_ready_time[p1] + i});
                            for (int i = 0; i < node.q2_1q_count; ++i) log.push_back({OpType::GATE_1Q, node.id, node.lq2, -1, p2, -1, -1, state.q_ready_time[p2] + i});
                        }

                        state.q_ready_time[p1] += node.q1_1q_count;
                        state.q_ready_time[p2] += node.q2_1q_count;
                        q_logical_depth[p1] += node.q1_1q_count;
                        q_logical_depth[p2] += node.q2_1q_count;

                        int chosen_res = -1;
                        const auto& p2_res = topo.qubits[p2].connected_resonators;
                        for (int r : topo.qubits[p1].connected_resonators) {
                            if (std::ranges::find(p2_res, r) != p2_res.end()) {
                                chosen_res = r;
                                break;
                            }
                        }

                        if (chosen_res != -1) {
                            double start_time = fast_max_sabre(state.q_ready_time[p1], state.q_ready_time[p2]);
                            bool q1_goes_in = (state.q_ready_time[p1] <= state.q_ready_time[p2]);

                            int p_in = q1_goes_in ? p1 : p2;
                            int p_stay = q1_goes_in ? p2 : p1;
                            int lq_in = q1_goes_in ? node.lq1 : node.lq2;
                            int lq_stay = q1_goes_in ? node.lq2 : node.lq1;

                            force_qubit_down(p_stay, node.id);
                            if (state.res_occupant[chosen_res] != -1 && state.res_occupant[chosen_res] != p_in) {
                                force_qubit_down(state.res_occupant[chosen_res], node.id);
                            }

                            if (state.q_on_resonator[p_in] != chosen_res) {
                                force_qubit_down(p_in, node.id);
                                start_time = fast_max_sabre(state.q_ready_time[p_in], fast_max_sabre(state.q_ready_time[p_stay], state.r_ready_time[chosen_res]));
                                if (is_final) log.push_back({OpType::MOVE, node.id, lq_in, -1, p_in, chosen_res + RES_OFFSET, -1, start_time});

                                out_moves++;
                                total_routing_cost += MOVE_COST;
                                start_time += MOVE_COST;
                                state.q_on_resonator[p_in] = chosen_res;
                                state.res_occupant[chosen_res] = p_in;
                                q_logical_depth[p_in] += 1;
                            } else {
                                start_time = fast_max_sabre(state.q_ready_time[p_in], fast_max_sabre(state.q_ready_time[p_stay], state.r_ready_time[chosen_res]));
                            }

                            if (is_final) log.push_back({OpType::GATE_2Q_LOCAL, node.id, lq_in, lq_stay, chosen_res + RES_OFFSET, p_stay, p_in, start_time});
                            state.q_ready_time[p_in] = start_time + 1.0;
                            state.q_ready_time[p_stay] = start_time + 1.0;
                            state.r_ready_time[chosen_res] = start_time + 1.0;

                            int gate_layer = fast_max_sabre(q_logical_depth[p_in], q_logical_depth[p_stay]);
                            q_logical_depth[p_in] = gate_layer + 1;
                            q_logical_depth[p_stay] = gate_layer + 1;
                        }
                    }

                    gates_routed++;
                    for (int next_id : node.next_gates) {
                        if (--deps[next_id] == 0) front_layer.push_back(next_id);
                    }
                }

                std::fill(decay.begin(), decay.end(), DECAY_RESET);
                continue;
            }

            for (int p : active_phys_qubits) in_active[p] = false;
            for (int g : extended_layer) in_extended[g] = false;
            active_phys_qubits.clear();
            extended_layer.clear();
            candidates.clear();

            for (int gate_id : front_layer) {
                const GateNode& node = current_dag.nodes[gate_id];
                if (node.lq1 != -1) {
                    int p1 = state.logical_to_phys_qubit[node.lq1];
                    if (!in_active[p1]) { in_active[p1] = true; active_phys_qubits.push_back(p1); }
                }
                if (node.lq2 != -1) {
                    int p2 = state.logical_to_phys_qubit[node.lq2];
                    if (!in_active[p2]) { in_active[p2] = true; active_phys_qubits.push_back(p2); }
                }

                for (int nxt : node.next_gates) {
                    if (!in_extended[nxt]) { in_extended[nxt] = true; extended_layer.push_back(nxt); }
                }
            }

            // KRYTYCZNA POPRAWKA GENEROWANIA KANDYDATÓW
            for (int p : active_phys_qubits) {
                for (int r : topo.qubits[p].connected_resonators) {
                    for (int neighbor_p : topo.resonators[r].connected_qubits) {
                        if (p == neighbor_p) continue; // Unikaj łączenia z samym sobą

                        // Ujednolicenie pary do sprawdzenia duplikatów
                        int u = std::min(p, neighbor_p);
                        int v = std::max(p, neighbor_p);

                        bool already_added = false;
                        for (const auto& c : candidates) {
                            if (c.p1 == u && c.p2 == v) { already_added = true; break; }
                        }
                        if (already_added) continue;

                        double score_F = 0.0;
                        double score_E = 0.0;

                        for (int gate_id : front_layer) {
                            const GateNode& n = current_dag.nodes[gate_id];
                            if (n.lq2 != -1) {
                                int lq1_p = state.logical_to_phys_qubit[n.lq1];
                                int lq2_p = state.logical_to_phys_qubit[n.lq2];
                                if (lq1_p == p) lq1_p = neighbor_p; else if (lq1_p == neighbor_p) lq1_p = p;
                                if (lq2_p == p) lq2_p = neighbor_p; else if (lq2_p == neighbor_p) lq2_p = p;
                                score_F += get_distance(lq1_p, lq2_p);
                            }
                        }

                        for (int gate_id : extended_layer) {
                            const GateNode& n = current_dag.nodes[gate_id];
                            if (n.lq2 != -1) {
                                int lq1_p = state.logical_to_phys_qubit[n.lq1];
                                int lq2_p = state.logical_to_phys_qubit[n.lq2];
                                if (lq1_p == p) lq1_p = neighbor_p; else if (lq1_p == neighbor_p) lq1_p = p;
                                if (lq2_p == p) lq2_p = neighbor_p; else if (lq2_p == neighbor_p) lq2_p = p;
                                score_E += get_distance(lq1_p, lq2_p);
                            }
                        }

                        double f_norm = front_layer.empty() ? 0 : (score_F / front_layer.size());
                        double e_norm = extended_layer.empty() ? 0 : (score_E / extended_layer.size());

                        double score = (f_norm + EXTENDED_LAYER_WEIGHT * e_norm) * fast_max_sabre(decay[u], decay[v]);
                        // Dodajemy ujednoliconą parę do wektora
                        candidates.push_back({u, v, r, score});
                    }
                }
            }

            SwapCandidate* best_candidate = nullptr;
            if (!candidates.empty()) {
                best_candidate = &candidates[0];
                for(size_t i = 1; i < candidates.size(); ++i) {
                    if(candidates[i].score < best_candidate->score) best_candidate = &candidates[i];
                }
            }

            if (best_candidate != nullptr) {
                int p1 = best_candidate->p1;
                int p2 = best_candidate->p2;
                int r_step = best_candidate->r_shared;

                force_qubit_down(p1, -1);
                force_qubit_down(p2, -1);

                double step_start = fast_max_sabre(state.q_ready_time[p1], fast_max_sabre(state.r_ready_time[r_step], state.q_ready_time[p2]));

                bool is_move = (state.phys_qubit_to_logical[p2] == -1) || (state.phys_qubit_to_logical[p1] == -1);
                double op_cost = is_move ? (2.0 * MOVE_COST) : SWAP_COST;

                int lq_1 = state.phys_qubit_to_logical[p1];
                int lq_2 = state.phys_qubit_to_logical[p2];

                if (is_final) {
                    if (state.phys_qubit_to_logical[p2] == -1) {
                        log.push_back({OpType::MOVE, -1, lq_1, -1, p1, r_step + RES_OFFSET, -1, step_start});
                        log.push_back({OpType::MOVE, -1, lq_1, -1, r_step + RES_OFFSET, p2, -1, step_start + MOVE_COST});
                    } else if (state.phys_qubit_to_logical[p1] == -1) {
                        log.push_back({OpType::MOVE, -1, lq_2, -1, p2, r_step + RES_OFFSET, -1, step_start});
                        log.push_back({OpType::MOVE, -1, lq_2, -1, r_step + RES_OFFSET, p1, -1, step_start + MOVE_COST});
                    } else {
                        log.push_back({OpType::SWAP, -1, lq_1, lq_2, p1, p2, r_step, step_start});
                    }
                }

                if(is_move) out_moves += 2; else out_swaps++;
                total_routing_cost += op_cost;

                std::swap(state.phys_qubit_to_logical[p1], state.phys_qubit_to_logical[p2]);
                if (lq_1 != -1) state.logical_to_phys_qubit[lq_1] = p2;
                if (lq_2 != -1) state.logical_to_phys_qubit[lq_2] = p1;

                int current_layer = fast_max_sabre(q_logical_depth[p1], q_logical_depth[p2]);
                int layers_added = is_move ? 2 : 1;
                q_logical_depth[p1] = current_layer + layers_added;
                q_logical_depth[p2] = current_layer + layers_added;

                state.q_ready_time[p1] = step_start + op_cost;
                state.r_ready_time[r_step] = step_start + op_cost;
                state.q_ready_time[p2] = step_start + op_cost;

                decay[p1] += DECAY_INCREMENT;
                decay[p2] += DECAY_INCREMENT;
            }
        }

        if (is_final) {
            for (size_t p = 0; p < topo.qubits.size(); ++p) force_qubit_down(p, -1);

            double max_physical_time = 0;
            int max_logical_depth = 0;
            for (double t : state.q_ready_time) max_physical_time = fast_max_sabre(max_physical_time, t);
            for (int d : q_logical_depth) max_logical_depth = fast_max_sabre(max_logical_depth, d);

            std::cout << "\n=== WYNIKI BENCHMARKU (SABRE High-Perf BiDirectional) ===" << std::endl;
            std::cout << "Wykorzystane operacje SWAP: " << out_swaps << std::endl;
            std::cout << "Wykorzystane operacje MOVE: " << out_moves << std::endl;
            std::cout << "Czas fizyczny obwodu (Makespan): " << max_physical_time << " [jednostek]" << std::endl;
            std::cout << "Głębokość warstwowa (Logical Depth): " << max_logical_depth << " [warstw]" << std::endl;
            std::cout << "Koszt routingu (Total Routing Cost): " << total_routing_cost << std::endl;
        }

        return state.logical_to_phys_qubit;
    }

public:
    SabreRouter(const BaseTopology& topology, const std::vector<Gate>& circuit, int qubits_cnt, std::string name)
        : topo(topology), original_circuit(circuit), num_qubits(qubits_cnt), name(name) {

        int p_cnt = topo.qubits.size();
        dist_cache.assign(p_cnt, std::vector<int>(p_cnt, -1));

        // Budujemy DAG dla Forward Pass
        dag = build_dag(original_circuit);

        // Budujemy odwrócony DAG dla Backward Pass
        std::vector<Gate> reversed_circuit = original_circuit;
        std::reverse(reversed_circuit.begin(), reversed_circuit.end());
        reversed_dag = build_dag(reversed_circuit);
    }

    void run_benchmark() {
        std::cout << "Uruchamiam benchmark SABRE (High-Performance Engine)... " << std::endl;
        auto start_time = std::chrono::high_resolution_clock::now();

        std::vector<int> initial_mapping(num_qubits);
        for (int i = 0; i < num_qubits; ++i) initial_mapping[i] = i;

        std::vector<int> post_forward_mapping = run_pass(dag, initial_mapping, false);
        std::vector<int> pre_final_mapping = run_pass(reversed_dag, post_forward_mapping, false);
        run_pass(dag, pre_final_mapping, true);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> execution_time = end_time - start_time;
        std::cout << "Czas wykonania całego procesu (F-B-F): " << execution_time.count() << " sekund" << std::endl;
        std::cout << "=========================================================\n" << std::endl;
    }
};

#endif //MAGISTERKA_SABREROUTER_H