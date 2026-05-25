//
// Created by skowr on 22.05.2026.
//

#include "SabreRouter.h"

SabreRouter::SabreRouter(const BaseTopology& topology, const std::vector<Gate>& circuit, int qubits_cnt, std::string name)
    : topo(topology), original_circuit(circuit), num_qubits(qubits_cnt), name(std::move(name)) {

    int p_cnt = topo.qubits.size();
    dist_cache.assign(p_cnt * p_cnt, -1);

    // OPTYMALIZACJA 1: Prekalkulacja fizycznych połączeń rezonatorów
    precomputed_swaps.resize(p_cnt);
    for (int p = 0; p < p_cnt; ++p) {
        for (int r : topo.qubits[p].connected_resonators) {
            for (int neighbor_p : topo.resonators[r].connected_qubits) {
                if (p != neighbor_p) {
                    precomputed_swaps[p].push_back({neighbor_p, r});
                }
            }
        }
    }

    dag = build_dag(original_circuit);

    std::vector<Gate> reversed_circuit = original_circuit;
    std::reverse(reversed_circuit.begin(), reversed_circuit.end());
    reversed_dag = build_dag(reversed_circuit);
}

int SabreRouter::get_distance(int p1, int p2) const {
    if (p1 == p2) return 0;
    int p_cnt = topo.qubits.size();
    int idx1 = p1 * p_cnt + p2;

    if (dist_cache[idx1] != -1) return dist_cache[idx1];
    int dist = topo.get_routing_path(p1, p2).resonators.size();

    dist_cache[idx1] = dist;
    dist_cache[p2 * p_cnt + p1] = dist;
    return dist;
}

SabreDAG SabreRouter::build_dag(const std::vector<Gate>& circuit) const {
    // ... Wnętrze funkcji (pozostawione bez zmian) ...
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

void SabreRouter::force_qubit_down(int p, int current_gate_id, MachineState& state, std::vector<int>& q_logical_depth,
                                   std::vector<RoutingRecord>& log, int& out_moves, double& total_routing_cost, bool is_final) const {
    // ... Wnętrze funkcji (pozostawione bez zmian) ...
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
}

void SabreRouter::process_executable_gates(const SabreDAG& current_dag, const std::vector<int>& executable_gates,
                                           MachineState& state, std::vector<int>& q_logical_depth, std::vector<RoutingRecord>& log,
                                           int& out_moves, double& total_routing_cost, bool is_final,
                                           int& gates_routed, std::vector<int>& deps, std::vector<int>& front_layer) const {
    // ... Wnętrze funkcji (pozostawione bez zmian) ...
    for (int gate_id : executable_gates) {
        const GateNode& node = current_dag.nodes[gate_id];
        int p1 = state.logical_to_phys_qubit[node.lq1];

        if (node.lq2 == -1) {
            if (node.q1_1q_count > 0) force_qubit_down(p1, node.id, state, q_logical_depth, log, out_moves, total_routing_cost, is_final);
            if (is_final) {
                for (int i = 0; i < node.q1_1q_count; ++i) {
                    log.push_back({OpType::GATE_1Q, node.id, node.lq1, -1, p1, -1, -1, state.q_ready_time[p1] + i});
                }
            }
            state.q_ready_time[p1] += node.q1_1q_count;
            q_logical_depth[p1] += node.q1_1q_count;
        } else {
            int p2 = state.logical_to_phys_qubit[node.lq2];
            if (node.q1_1q_count > 0) force_qubit_down(p1, node.id, state, q_logical_depth, log, out_moves, total_routing_cost, is_final);
            if (node.q2_1q_count > 0) force_qubit_down(p2, node.id, state, q_logical_depth, log, out_moves, total_routing_cost, is_final);

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

                force_qubit_down(p_stay, node.id, state, q_logical_depth, log, out_moves, total_routing_cost, is_final);
                if (state.res_occupant[chosen_res] != -1 && state.res_occupant[chosen_res] != p_in) {
                    force_qubit_down(state.res_occupant[chosen_res], node.id, state, q_logical_depth, log, out_moves, total_routing_cost, is_final);
                }

                if (state.q_on_resonator[p_in] != chosen_res) {
                    force_qubit_down(p_in, node.id, state, q_logical_depth, log, out_moves, total_routing_cost, is_final);
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
}

// ZMIANA: std::vector<bool> -> std::vector<uint8_t>
void SabreRouter::update_active_and_extended_layers(const SabreDAG& current_dag, const std::vector<int>& front_layer,
                                                    const MachineState& state, std::vector<int>& active_phys_qubits,
                                                    std::vector<uint8_t>& in_active, std::vector<int>& extended_layer,
                                                    std::vector<uint8_t>& in_extended) const {
    for (int p : active_phys_qubits) in_active[p] = 0;
    for (int g : extended_layer) in_extended[g] = 0;
    active_phys_qubits.clear();
    extended_layer.clear();

    for (int gate_id : front_layer) {
        const GateNode& node = current_dag.nodes[gate_id];
        if (node.lq1 != -1) {
            int p1 = state.logical_to_phys_qubit[node.lq1];
            if (!in_active[p1]) { in_active[p1] = 1; active_phys_qubits.push_back(p1); }
        }
        if (node.lq2 != -1) {
            int p2 = state.logical_to_phys_qubit[node.lq2];
            if (!in_active[p2]) { in_active[p2] = 1; active_phys_qubits.push_back(p2); }
        }

        for (int nxt : node.next_gates) {
            if (!in_extended[nxt]) { in_extended[nxt] = 1; extended_layer.push_back(nxt); }
        }
    }
}

SwapCandidate SabreRouter::find_best_swap_candidate(const SabreDAG& current_dag, const MachineState& state,
                                                    const std::vector<int>& front_layer, const std::vector<int>& extended_layer,
                                                    const std::vector<int>& active_phys_qubits, const std::vector<double>& decay) const {
    std::vector<SwapCandidate> candidates;
    candidates.reserve(topo.qubits.size() * 4);

    int p_cnt = topo.qubits.size();
    std::vector<uint8_t> seen_pairs(p_cnt * p_cnt, 0); // O(1) tablica sprawdzen

    // OPTYMALIZACJA 2: Zamiast fizycznych pair, zapamiętujemy LOGICZNE
    // Dzięki temu możemy wyeliminować instrukcje 'if' podczas ewaluacji
    std::vector<std::pair<int, int>> front_logical_pairs;
    front_logical_pairs.reserve(front_layer.size());
    for (int gate_id : front_layer) {
        const GateNode& n = current_dag.nodes[gate_id];
        if (n.lq2 != -1) front_logical_pairs.emplace_back(n.lq1, n.lq2);
    }

    std::vector<std::pair<int, int>> ext_logical_pairs;
    ext_logical_pairs.reserve(extended_layer.size());
    for (int gate_id : extended_layer) {
        const GateNode& n = current_dag.nodes[gate_id];
        if (n.lq2 != -1) ext_logical_pairs.emplace_back(n.lq1, n.lq2);
    }

    // Ekstremalnie szybka lokalna kopia bieżącego mapowania dla O(1) dostępu w cache
    std::vector<int> current_l2p = state.logical_to_phys_qubit;

    // OPTYMALIZACJA 1: Korzystamy z prekompilowanego, płaskiego wektora ścieżek SwapEdge
    for (int p : active_phys_qubits) {
        for (const SwapEdge& edge : precomputed_swaps[p]) {
            int neighbor_p = edge.neighbor;
            int r = edge.resonator;

            int u = std::min(p, neighbor_p);
            int v = std::max(p, neighbor_p);

            int pair_idx = u * p_cnt + v;
            if (seen_pairs[pair_idx]) continue;
            seen_pairs[pair_idx] = 1;

            // Sprawdzamy jakie kubity logiczne zamierzamy podmienić
            int lq_p = state.phys_qubit_to_logical[p];
            int lq_neighbor = state.phys_qubit_to_logical[neighbor_p];

            // Wirtualny swap - operacja w pełni O(1) na tablicy intów
            if (lq_p != -1) current_l2p[lq_p] = neighbor_p;
            if (lq_neighbor != -1) current_l2p[lq_neighbor] = p;

            double score_F = 0.0;
            // Branchless pętla - procesor z łatwością wektoryzuje sumowanie
            for (const auto& pr : front_logical_pairs) {
                score_F += get_distance(current_l2p[pr.first], current_l2p[pr.second]);
            }

            double score_E = 0.0;
            // Branchless pętla
            for (const auto& pr : ext_logical_pairs) {
                score_E += get_distance(current_l2p[pr.first], current_l2p[pr.second]);
            }

            // Cofnięcie wirtualnego swapa
            if (lq_p != -1) current_l2p[lq_p] = p;
            if (lq_neighbor != -1) current_l2p[lq_neighbor] = neighbor_p;

            double f_norm = front_layer.empty() ? 0 : (score_F / front_layer.size());
            double e_norm = extended_layer.empty() ? 0 : (score_E / extended_layer.size());

            double score = (f_norm + EXTENDED_LAYER_WEIGHT * e_norm) * fast_max_sabre(decay[u], decay[v]);
            candidates.push_back({u, v, r, score});
        }
    }

    SwapCandidate best_candidate = {-1, -1, -1, std::numeric_limits<double>::max()};
    for (const auto& candidate : candidates) {
        if (candidate.score < best_candidate.score) {
            best_candidate = candidate;
        }
    }
    return best_candidate;
}

void SabreRouter::apply_swap_candidate(const SwapCandidate& candidate, MachineState& state, std::vector<int>& q_logical_depth,
                                       std::vector<RoutingRecord>& log, std::vector<double>& decay, int& out_moves,
                                       int& out_swaps, double& total_routing_cost, bool is_final) const {
    // ... Wnętrze funkcji (pozostawione bez zmian) ...
    int p1 = candidate.p1;
    int p2 = candidate.p2;
    int r_step = candidate.r_shared;

    force_qubit_down(p1, -1, state, q_logical_depth, log, out_moves, total_routing_cost, is_final);
    force_qubit_down(p2, -1, state, q_logical_depth, log, out_moves, total_routing_cost, is_final);

    if (state.res_occupant[r_step] != -1) {
        force_qubit_down(state.res_occupant[r_step], -1, state, q_logical_depth, log, out_moves, total_routing_cost, is_final);
    }

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

std::vector<int> SabreRouter::run_pass(const SabreDAG& current_dag, const std::vector<int>& initial_l2p, bool is_final) {
    int n_gates = current_dag.nodes.size();

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
    int gates_routed = 0;

    std::vector<int> front_layer;
    front_layer.reserve(num_qubits);
    std::vector<int> deps = current_dag.initial_deps;
    for (int i = 0; i < n_gates; ++i) {
        if (deps[i] == 0) front_layer.push_back(i);
    }

    std::vector<double> decay(topo.qubits.size(), 1.0);
    std::vector<int> executable_gates;
    executable_gates.reserve(num_qubits);

    // ZMIANA: std::vector<bool> zastąpione uint8_t
    std::vector<uint8_t> in_active(topo.qubits.size(), 0);
    std::vector<int> active_phys_qubits;
    active_phys_qubits.reserve(num_qubits);

    std::vector<uint8_t> in_extended(n_gates, 0);
    std::vector<int> extended_layer;
    extended_layer.reserve(num_qubits * 4);

    int debug_iterations = 0;
    const int MAX_ITERATIONS = n_gates * 100;

    while (gates_routed < n_gates) {
        if (++debug_iterations > MAX_ITERATIONS) {
            std::cerr << "\n[CRITICAL ERROR] SABRE ugrzązł w nieskończonej pętli!\n";
            std::cerr << "Bramek ułożonych: " << gates_routed << " / " << n_gates << "\n";
            std::cerr << "Front layer size: " << front_layer.size() << "\n";
            break;
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
            process_executable_gates(current_dag, executable_gates, state, q_logical_depth, log, out_moves, total_routing_cost, is_final, gates_routed, deps, front_layer);
            std::fill(decay.begin(), decay.end(), DECAY_RESET);
            continue;
        }

        update_active_and_extended_layers(current_dag, front_layer, state, active_phys_qubits, in_active, extended_layer, in_extended);

        SwapCandidate best_candidate = find_best_swap_candidate(current_dag, state, front_layer, extended_layer, active_phys_qubits, decay);

        if (best_candidate.p1 != -1) {
            apply_swap_candidate(best_candidate, state, q_logical_depth, log, decay, out_moves, out_swaps, total_routing_cost, is_final);
        }
    }

    if (is_final) {
        for (size_t p = 0; p < topo.qubits.size(); ++p) force_qubit_down(p, -1, state, q_logical_depth, log, out_moves, total_routing_cost, is_final);

        std::stable_sort(log.begin(), log.end(), [](const RoutingRecord& a, const RoutingRecord& b) {
            return a.start_time < b.start_time;
        });

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

        // print_optimized_circuit(log, num_qubits);

        std::cout << "=== WALIDACJA OBWODU LOGICZNEGO ===\n";
        if (!validate_logical_circuit(original_circuit, log)) {
            exit(-1);
        }
        std::cout << "===================================\n\n";
    }

    return state.logical_to_phys_qubit;
}

void SabreRouter::run_benchmark() {
    // ... Wnętrze funkcji (pozostawione bez zmian) ...
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

bool SabreRouter::validate_logical_circuit(const std::vector<Gate>& orig_circ, const std::vector<RoutingRecord>& log) const {
    // 1. Zbieramy oczekiwaną sekwencję ID z DAG-a
    std::vector<std::vector<int>> expected(num_qubits);
    for (const auto& node : dag.nodes) {
        if (node.lq1 != -1) expected[node.lq1].push_back(node.id);
        if (node.lq2 != -1) expected[node.lq2].push_back(node.id);
    }

    // 2. Zbieramy faktyczną sekwencję ID z logów rutingu
    std::vector<std::vector<int>> actual(num_qubits);
    for (const auto& op : log) {
        if (op.op_type == OpType::GATE_1Q || op.op_type == OpType::GATE_2Q_LOCAL || op.op_type == OpType::GATE_2Q_GLOBAL) {

            if (op.lq1 != -1 && (actual[op.lq1].empty() || actual[op.lq1].back() != op.gate_id)) {
                actual[op.lq1].push_back(op.gate_id);
            }
            if (op.lq2 != -1 && (actual[op.lq2].empty() || actual[op.lq2].back() != op.gate_id)) {
                actual[op.lq2].push_back(op.gate_id);
            }
        }
    }

    // 3. Weryfikacja
    bool is_valid = true;
    for (int i = 0; i < num_qubits; ++i) {
        if (expected[i] != actual[i]) {
            is_valid = false;
            std::cout << "[BŁĄD WALIDACJI] Niezgodność na kubicie logicznym lq" << i << "!\n";
            std::cout << "  Oczekiwano ID: ";
            for (int id : expected[i]) std::cout << id << " ";
            std::cout << "\n  Otrzymano ID:  ";
            for (int id : actual[i]) std::cout << id << " ";
            std::cout << "\n";
        }
    }

    if (is_valid) {
        std::cout << "[OK] Walidacja pomyślna. Czystość logiczna (ID bramek z DAG) zachowana na 100% (SABRE)!\n";
    }
    return is_valid;
}

void SabreRouter::print_optimized_circuit(const std::vector<RoutingRecord>& log, int num_logical_qubits) const {
    // ... Wnętrze funkcji (pozostawione bez zmian) ...
    if (num_logical_qubits <= 0) num_logical_qubits = num_qubits;
    if (log.empty()) { std::cout << "Zoptymalizowany obwód jest pusty!!\n"; return; }

    int n_q = topo.qubits.size();
    std::vector<std::stringstream> lines_q(n_q);

    for (int i = 0; i < n_q; ++i) {
        lines_q[i] << "q";
        if (i < 10) lines_q[i] << "0";
        lines_q[i] << i;
        if (i < num_logical_qubits) lines_q[i] << " [L]: ";
        else lines_q[i] << " [E]: ";
    }

    for (size_t idx = 0; idx < log.size(); ++idx) {
        const auto& move = log[idx];

        int max_len = 0;
        for (const auto& ss : lines_q) {
            int len = ss.str().length();
            if (len > max_len) max_len = len;
        }
        for (auto& ss : lines_q) {
            int pad = max_len - ss.str().length();
            for (int p = 0; p < pad; ++p) ss << "-";
        }

        int l1 = move.loc1;
        int l2 = move.loc2;

        if (move.op_type == OpType::GATE_1Q) {
            if (l1 >= 0 && l1 < n_q) lines_q[l1] << "-[ U3 ]-";
            for (int i = 0; i < n_q; ++i) {
                if (i != l1) lines_q[i] << "--------";
            }
        } else if (move.op_type == OpType::SWAP) {
            int min_q = std::min(l1, l2);
            int max_q = fast_max_sabre(l1, l2);
            for (int i = 0; i < n_q; ++i) {
                if (i == l1 || i == l2) lines_q[i] << "--x--";
                else if (i > min_q && i < max_q) lines_q[i] << "--|--";
                else lines_q[i] << "-----";
            }
        } else if (move.op_type == OpType::MOVE) {
            if (l1 < RES_OFFSET && l2 < RES_OFFSET && l2 >= 0) {
                int min_q = std::min(l1, l2);
                int max_q = fast_max_sabre(l1, l2);
                for (int i = 0; i < n_q; ++i) {
                    if (i == l1 || i == l2) lines_q[i] << "--m--";
                    else if (i > min_q && i < max_q) lines_q[i] << "--|--";
                    else lines_q[i] << "-----";
                }
            } else {
                int q_idx = (l1 < RES_OFFSET) ? l1 : l2;
                for (int i = 0; i < n_q; ++i) {
                    if (i == q_idx) lines_q[i] << "--m--";
                    else lines_q[i] << "-----";
                }
            }
        } else if (move.op_type == OpType::GATE_2Q_LOCAL || move.op_type == OpType::GATE_2Q_GLOBAL) {
            int r_idx = -1, q_stay = -1;
            if (l1 >= RES_OFFSET) { r_idx = l1 - RES_OFFSET; q_stay = l2; }
            else if (l2 >= RES_OFFSET) { r_idx = l2 - RES_OFFSET; q_stay = l1; }

            int q_in = move.res_used;

            int min_q = std::min(q_in, q_stay);
            int max_q = fast_max_sabre(q_in, q_stay);

            bool q_in_is_control = false;
            if (move.gate_id != -1 && move.gate_id < dag.nodes.size()) {
                q_in_is_control = (move.lq1 == dag.nodes[move.gate_id].lq1);
            }

            std::string ctrl_char = (move.op_type == OpType::GATE_2Q_LOCAL) ? "--o--" : "--*--";
            std::string targ_char = (move.op_type == OpType::GATE_2Q_LOCAL) ? "-(+)-" : "-[+]-";

            for (int i = 0; i < n_q; ++i) {
                if (i == q_in) {
                    lines_q[i] << (q_in_is_control ? ctrl_char : targ_char);
                }
                else if (i == q_stay) {
                    lines_q[i] << (q_in_is_control ? targ_char : ctrl_char);
                }
                else if (i > min_q && i < max_q && q_in != -1) {
                    lines_q[i] << "--|--";
                }
                else {
                    lines_q[i] << "-----";
                }
            }
        }
    }

    std::cout << "\n=== WIZUALIZACJA ZOPTYMALIZOWANEGO OBWODU (SABRE) ===\n";
    std::cout << "Legenda bramek: -[ U3 ]- [1Q] | --o--/-(+)- [2Q Lokalna] | --*--/-[+]- [2Q Globalna]\n";
    std::cout << "Operacje: --x-- [SWAP] | --m-- [MOVE / Wejscie na rezonator]\n\n";
    for (int i = 0; i < n_q; ++i) { std::cout << lines_q[i].str() << "\n"; }
    std::cout << "=====================================================\n\n";
}