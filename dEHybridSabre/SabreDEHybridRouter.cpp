#include "SabreDEHybridRouter.h"
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include <fstream>
#include <queue>
#include <sstream>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

inline double fast_max_hyb(double a, double b) { return (a > b) ? a : b; }

SabreDEHybridRouter::SabreDEHybridRouter(const BaseTopology& topology, const std::vector<Gate>& circuit, int qubits_cnt, std::string name)
    : topo(topology), original_circuit(circuit), num_qubits(qubits_cnt), name(std::move(name)) {
    int p_cnt = topo.qubits.size();
    int r_cnt = topo.resonators.size();

    // Spłaszczony bufor cache dla odległości (p_cnt x p_cnt)
    dist_cache.assign(p_cnt * p_cnt, -1);

    dag = build_dag(original_circuit);
    int g_cnt = dag.nodes.size();

    std::vector<Gate> reversed_circuit = original_circuit;
    std::reverse(reversed_circuit.begin(), reversed_circuit.end());
    reversed_dag = build_dag(reversed_circuit);

    // Pre-alokacja workspace, aby uniknąć re-alokacji podczas ewaluacji
    workspace.l2p.reserve(p_cnt);
    workspace.p2l.reserve(p_cnt);
    workspace.q_on_res.reserve(p_cnt);
    workspace.res_occ.reserve(r_cnt);
    workspace.q_depth.reserve(p_cnt);
    workspace.q_time.reserve(p_cnt);
    workspace.r_time.reserve(r_cnt);
    workspace.front_layer.reserve(g_cnt);
    workspace.exec_gates.reserve(p_cnt);
    workspace.ext_layer.reserve(g_cnt);
    workspace.in_ext.reserve(g_cnt);
    workspace.decay.reserve(p_cnt);

    // Nowe rezerwacje
    workspace.shared_res.reserve(r_cnt);
    workspace.in_active.reserve(p_cnt);
    workspace.active_phys.reserve(p_cnt);
    workspace.candidates.reserve(p_cnt * 4);
}

int SabreDEHybridRouter::get_distance(int p1, int p2) const {
    if (p1 == p2) return 0;
    int p_cnt = topo.qubits.size();
    int idx = p1 * p_cnt + p2;
    if (dist_cache[idx] != -1) return dist_cache[idx];

    int dist = topo.get_routing_path(p1, p2).resonators.size();
    dist_cache[idx] = dist;
    dist_cache[p2 * p_cnt + p1] = dist;
    return dist;
}

HybridSabreDAG SabreDEHybridRouter::build_dag(const std::vector<Gate>& circuit) const {
    HybridSabreDAG local_dag;
    std::vector<int> last_gate(num_qubits, -1), pending_1q(num_qubits, 0);
    int group_id = 0;

    for (const auto& g : circuit) {
        if (g.type == "1Q") pending_1q[g.target_q1]++;
        else {
            HybridGateNode node{group_id++, g.target_q1, g.target_q2, pending_1q[g.target_q1], pending_1q[g.target_q2], 0, {}};
            pending_1q[g.target_q1] = 0; pending_1q[g.target_q2] = 0;
            if (last_gate[node.lq1] != -1) { local_dag.nodes[last_gate[node.lq1]].next_gates.push_back(node.id); node.deps_count++; }
            if (last_gate[node.lq2] != -1 && last_gate[node.lq2] != last_gate[node.lq1]) {
                local_dag.nodes[last_gate[node.lq2]].next_gates.push_back(node.id); node.deps_count++;
            }
            local_dag.nodes.push_back(node);
            last_gate[node.lq1] = node.id; last_gate[node.lq2] = node.id;
        }
    }
    for (int i = 0; i < num_qubits; ++i) {
        if (pending_1q[i] > 0) {
            HybridGateNode node{group_id++, i, -1, pending_1q[i], 0, 0, {}};
            if (last_gate[i] != -1) { local_dag.nodes[last_gate[i]].next_gates.push_back(node.id); node.deps_count++; }
            local_dag.nodes.push_back(node);
            last_gate[i] = node.id; pending_1q[i] = 0;
        }
    }
    local_dag.initial_deps.resize(local_dag.nodes.size());
    for (size_t i = 0; i < local_dag.nodes.size(); ++i) local_dag.initial_deps[i] = local_dag.nodes[i].deps_count;
    return local_dag;
}

// =========================================================================
//                              SABRE SEED
// =========================================================================
// [Tutaj bez zmian funkcje sabre_force_qubit_down, sabre_execute_gates,
//  sabre_route_step i run_sabre_pass - pozostają identyczne z oryginałem]

// =========================================================================
//                              DE SEED
// =========================================================================

template <bool LOGGING>
void SabreDEHybridRouter::de_force_down(int& p, int gate_id, DERunContext& ctx, const double* genes, int n_gates, std::vector<RoutingRecord>* out_log, bool can_lazy) {
    int r = ctx.q_on_res[p];
    if (r != -1) {
        int target_p = p;

        if (can_lazy && gate_id >= 0) {
            int valid_count = 1;
            const auto& conn = topo.resonators[r].connected_qubits;
            for (int c : conn) if (c != p && ctx.p2l[c] == -1 && ctx.q_on_res[c] == -1) valid_count++;

            if (valid_count > 1) {
                int gene_idx = std::min(static_cast<int>(genes[n_gates * 3 + gate_id] * valid_count), valid_count - 1);
                if (gene_idx > 0) {
                    int c_idx = 1;
                    for (int c : conn) if (c != p && ctx.p2l[c] == -1 && ctx.q_on_res[c] == -1) {
                        if (c_idx == gene_idx) { target_p = c; break; }
                        c_idx++;
                    }
                }
            }
        }

        double start = fast_max_hyb(ctx.q_time[p], ctx.r_time[r]);
        if (target_p != p) start = fast_max_hyb(start, ctx.q_time[target_p]);

        if constexpr (LOGGING) out_log->push_back({OpType::MOVE, gate_id, ctx.p2l[p], -1, r + RES_OFFSET, target_p, -1, start});
        ctx.out_moves++; ctx.total_routing_cost += MOVE_COST;

        ctx.q_time[target_p] = start + MOVE_COST;
        ctx.r_time[r] = start + MOVE_COST;

        if (target_p != p) {
            ctx.p2l[target_p] = ctx.p2l[p]; ctx.l2p[ctx.p2l[p]] = target_p;
            ctx.p2l[p] = -1;
            ctx.q_depth[target_p] = ctx.q_depth[p] + 1; ctx.q_depth[p] = 0;
        } else {
            ctx.q_depth[p] += 1;
        }

        ctx.q_on_res[p] = -1; ctx.q_on_res[target_p] = -1; ctx.res_occ[r] = -1;
        p = target_p;
    }
}

template <bool LOGGING>
void SabreDEHybridRouter::de_execute_gates(const HybridSabreDAG& current_dag, std::vector<int>& exec_gates, std::vector<int>& front_layer, std::vector<int>& deps, DERunContext& ctx, const double* genes, int n_gates, std::vector<RoutingRecord>* out_log, bool can_lazy) {
    std::sort(exec_gates.begin(), exec_gates.end(), [&](int a, int b) { return genes[a] > genes[b]; });

    for (int gate_id : exec_gates) {
        const auto& node = current_dag.nodes[gate_id];
        int p1 = ctx.l2p[node.lq1];

        if (node.lq2 == -1) {
            if (node.q1_1q_count > 0) {
                de_force_down<LOGGING>(p1, node.id, ctx, genes, n_gates, out_log, can_lazy);
                if constexpr (LOGGING) {
                    for (int k = 0; k < node.q1_1q_count; ++k) out_log->push_back({OpType::GATE_1Q, node.id, node.lq1, -1, p1, -1, -1, ctx.q_time[p1] + k});
                }
                ctx.q_time[p1] += node.q1_1q_count; ctx.q_depth[p1] += node.q1_1q_count;
            }
        } else {
            int p2 = ctx.l2p[node.lq2];
            if (node.q1_1q_count > 0) de_force_down<LOGGING>(p1, node.id, ctx, genes, n_gates, out_log, can_lazy);
            if (node.q2_1q_count > 0) de_force_down<LOGGING>(p2, node.id, ctx, genes, n_gates, out_log, can_lazy);

            if constexpr (LOGGING) {
                for (int k = 0; k < node.q1_1q_count; ++k) out_log->push_back({OpType::GATE_1Q, node.id, node.lq1, -1, p1, -1, -1, ctx.q_time[p1] + k});
                for (int k = 0; k < node.q2_1q_count; ++k) out_log->push_back({OpType::GATE_1Q, node.id, node.lq2, -1, p2, -1, -1, ctx.q_time[p2] + k});
            }
            ctx.q_time[p1] += node.q1_1q_count; ctx.q_time[p2] += node.q2_1q_count;
            ctx.q_depth[p1] += node.q1_1q_count; ctx.q_depth[p2] += node.q2_1q_count;

            ctx.ws.shared_res.clear();
            const auto& p2_res = topo.qubits[p2].connected_resonators;
            for (int r : topo.qubits[p1].connected_resonators) {
                if (std::find(p2_res.begin(), p2_res.end(), r) != p2_res.end()) ctx.ws.shared_res.push_back(r);
            }

            if (!ctx.ws.shared_res.empty()) {
                int chosen_res = ctx.ws.shared_res[std::min(static_cast<int>(genes[n_gates + node.id] * ctx.ws.shared_res.size()), (int)ctx.ws.shared_res.size() - 1)];
                bool q1_in = (node.q1_1q_count < node.q2_1q_count) ? true : (node.q1_1q_count > node.q2_1q_count) ? false : (genes[n_gates * 2 + node.id] < 0.5);

                int p_in = q1_in ? p1 : p2; int p_stay = q1_in ? p2 : p1;

                de_force_down<LOGGING>(p_stay, node.id, ctx, genes, n_gates, out_log, can_lazy);
                if (ctx.res_occ[chosen_res] != -1 && ctx.res_occ[chosen_res] != p_in) {
                    int r_occ = ctx.res_occ[chosen_res];
                    de_force_down<LOGGING>(r_occ, node.id, ctx, genes, n_gates, out_log, can_lazy);
                }

                double start = fast_max_hyb(ctx.q_time[p1], ctx.q_time[p2]);
                if (ctx.q_on_res[p_in] != chosen_res) {
                    de_force_down<LOGGING>(p_in, node.id, ctx, genes, n_gates, out_log, can_lazy);
                    start = fast_max_hyb(ctx.q_time[p_in], fast_max_hyb(ctx.q_time[p_stay], ctx.r_time[chosen_res]));
                    if constexpr (LOGGING) out_log->push_back({OpType::MOVE, node.id, ctx.p2l[p_in], -1, p_in, chosen_res + RES_OFFSET, -1, start});
                    ctx.out_moves++; ctx.total_routing_cost += MOVE_COST; start += MOVE_COST;
                    ctx.q_on_res[p_in] = chosen_res; ctx.res_occ[chosen_res] = p_in; ctx.q_depth[p_in]++;
                } else {
                    start = fast_max_hyb(ctx.q_time[p_in], fast_max_hyb(ctx.q_time[p_stay], ctx.r_time[chosen_res]));
                }

                if constexpr (LOGGING) out_log->push_back({OpType::GATE_2Q_LOCAL, node.id, ctx.p2l[p_in], ctx.p2l[p_stay], chosen_res + RES_OFFSET, p_stay, p_in, start});

                ctx.q_time[p_in] = ctx.q_time[p_stay] = ctx.r_time[chosen_res] = start + 1.0;
                int glayer = fast_max_hyb(ctx.q_depth[p_in], ctx.q_depth[p_stay]);
                ctx.q_depth[p_in] = ctx.q_depth[p_stay] = glayer + 1;
            }
        }
        for (int nxt : node.next_gates) if (--deps[nxt] == 0) front_layer.push_back(nxt);
    }
}

template <bool LOGGING>
bool SabreDEHybridRouter::de_route_step(const HybridSabreDAG& current_dag, const std::vector<int>& front_layer, const std::vector<int>& ext_layer, std::vector<double>& decay, DERunContext& ctx, const double* genes, int n_gates, std::vector<RoutingRecord>* out_log) {
    ctx.ws.in_active.assign(topo.qubits.size(), 0);
    ctx.ws.active_phys.clear();

    for (int gid : front_layer) {
        const auto& node = current_dag.nodes[gid];
        if (node.lq1 != -1 && ctx.l2p[node.lq1] != -1 && !ctx.ws.in_active[ctx.l2p[node.lq1]]) {
            ctx.ws.in_active[ctx.l2p[node.lq1]] = 1;
            ctx.ws.active_phys.push_back(ctx.l2p[node.lq1]);
        }
        if (node.lq2 != -1 && ctx.l2p[node.lq2] != -1 && !ctx.ws.in_active[ctx.l2p[node.lq2]]) {
            ctx.ws.in_active[ctx.l2p[node.lq2]] = 1;
            ctx.ws.active_phys.push_back(ctx.l2p[node.lq2]);
        }
    }

    ctx.ws.candidates.clear();

    for (int p : ctx.ws.active_phys) {
        for (int r : topo.qubits[p].connected_resonators) {
            for (int nb : topo.resonators[r].connected_qubits) {
                if (p == nb) continue;
                int u = std::min(p, nb), v = std::max(p, nb);
                bool added = false;
                for (const auto& c : ctx.ws.candidates) if (c.p1 == u && c.p2 == v) { added = true; break; }
                if (added) continue;

                double sF = 0, sE = 0;
                for (int gid : front_layer) {
                    if (current_dag.nodes[gid].lq2 != -1) {
                        int l1 = ctx.l2p[current_dag.nodes[gid].lq1], l2 = ctx.l2p[current_dag.nodes[gid].lq2];
                        if (l1 == p) l1 = nb; else if (l1 == nb) l1 = p;
                        if (l2 == p) l2 = nb; else if (l2 == nb) l2 = p;
                        sF += get_distance(l1, l2);
                    }
                }
                for (int gid : ext_layer) {
                    if (current_dag.nodes[gid].lq2 != -1) {
                        int l1 = ctx.l2p[current_dag.nodes[gid].lq1], l2 = ctx.l2p[current_dag.nodes[gid].lq2];
                        if (l1 == p) l1 = nb; else if (l1 == nb) l1 = p;
                        if (l2 == p) l2 = nb; else if (l2 == nb) l2 = p;
                        sE += get_distance(l1, l2);
                    }
                }
                double mutation = 0.95 + (genes[u % n_gates] * 0.1);
                double score = ((front_layer.empty() ? 0 : sF / front_layer.size()) + EXTENDED_LAYER_WEIGHT * (ext_layer.empty() ? 0 : sE / ext_layer.size())) * fast_max_hyb(decay[u], decay[v]) * mutation;
                ctx.ws.candidates.push_back({u, v, r, score});
            }
        }
    }

    SwapCand* best = nullptr;
    if (!ctx.ws.candidates.empty()) {
        best = &ctx.ws.candidates[0];
        for (size_t i = 1; i < ctx.ws.candidates.size(); ++i) if (ctx.ws.candidates[i].score < best->score) best = &ctx.ws.candidates[i];
    }

    if (best) {
        int p1 = best->p1, p2 = best->p2, r_step = best->r;
        de_force_down<LOGGING>(p1, -1, ctx, genes, n_gates, out_log, false);
        de_force_down<LOGGING>(p2, -1, ctx, genes, n_gates, out_log, false);

        double step_start = fast_max_hyb(ctx.q_time[p1], fast_max_hyb(ctx.r_time[r_step], ctx.q_time[p2]));
        bool is_move = (ctx.p2l[p2] == -1) || (ctx.p2l[p1] == -1);
        double cost = is_move ? 2.0 * MOVE_COST : SWAP_COST;

        if constexpr (LOGGING) {
            if (is_move) {
                if (ctx.p2l[p2] == -1) {
                    out_log->push_back({OpType::MOVE, -1, ctx.p2l[p1], -1, p1, r_step + RES_OFFSET, -1, step_start});
                    out_log->push_back({OpType::MOVE, -1, ctx.p2l[p1], -1, r_step + RES_OFFSET, p2, -1, step_start + MOVE_COST});
                } else {
                    out_log->push_back({OpType::MOVE, -1, ctx.p2l[p2], -1, p2, r_step + RES_OFFSET, -1, step_start});
                    out_log->push_back({OpType::MOVE, -1, ctx.p2l[p2], -1, r_step + RES_OFFSET, p1, -1, step_start + MOVE_COST});
                }
            } else out_log->push_back({OpType::SWAP, -1, ctx.p2l[p1], ctx.p2l[p2], p1, p2, r_step, step_start});
        }

        if (is_move) ctx.out_moves += 2; else ctx.out_swaps++;
        ctx.total_routing_cost += cost;

        std::swap(ctx.p2l[p1], ctx.p2l[p2]);
        if (ctx.p2l[p1] != -1) ctx.l2p[ctx.p2l[p1]] = p1;
        if (ctx.p2l[p2] != -1) ctx.l2p[ctx.p2l[p2]] = p2;

        int curr_lay = fast_max_hyb(ctx.q_depth[p1], ctx.q_depth[p2]) + (is_move ? 2 : 1);
        ctx.q_depth[p1] = ctx.q_depth[p2] = curr_lay;
        ctx.q_time[p1] = ctx.q_time[p2] = ctx.r_time[r_step] = step_start + cost;
        decay[p1] += 0.2; decay[p2] += 0.2;
        return true;
    }
    return false;
}

template <bool LOGGING>
double SabreDEHybridRouter::evaluate_fitness_impl(const double* genes, int& out_swaps, int& out_moves, double& max_physical_time, int& max_logical_depth, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard) {
    int n_gates = dag.nodes.size();
    int p_cnt = topo.qubits.size();
    int r_cnt = topo.resonators.size();
    out_swaps = 0; out_moves = 0; max_physical_time = 0.0; max_logical_depth = 0;

    // Praca na pre-alokowanym obiekcie (recykling zasobów z workspace)
    DERunContext ctx(workspace);
    ctx.l2p = saved_initial_l2p;
    ctx.p2l.assign(p_cnt, -1);
    for (int i = 0; i < num_qubits; ++i) if (ctx.l2p[i] != -1) ctx.p2l[ctx.l2p[i]] = i;

    ctx.q_on_res.assign(p_cnt, -1);
    ctx.res_occ.assign(r_cnt, -1);
    ctx.q_time.assign(p_cnt, 0.0);
    ctx.r_time.assign(r_cnt, 0.0);
    ctx.q_depth.assign(p_cnt, 0);

    if (out_log) { out_log->clear(); }

    workspace.front_layer.clear();
    workspace.deps = dag.initial_deps;
    for (int i = 0; i < n_gates; ++i) if (workspace.deps[i] == 0) workspace.front_layer.push_back(i);

    workspace.decay.assign(p_cnt, 1.0);
    int gates_routed = 0;
    workspace.in_ext.assign(n_gates, false);

    while (gates_routed < n_gates) {
        workspace.exec_gates.clear();
        for (size_t i = 0; i < workspace.front_layer.size(); ) {
            int g_id = workspace.front_layer[i];
            const auto& node = dag.nodes[g_id];
            bool exec = (node.lq2 == -1) || (ctx.l2p[node.lq1] != -1 && ctx.l2p[node.lq2] != -1 && get_distance(ctx.l2p[node.lq1], ctx.l2p[node.lq2]) <= 1);

            if (exec) { workspace.exec_gates.push_back(g_id); workspace.front_layer[i] = workspace.front_layer.back(); workspace.front_layer.pop_back(); }
            else ++i;
        }

        if (!workspace.exec_gates.empty()) {
            de_execute_gates<LOGGING>(dag, workspace.exec_gates, workspace.front_layer, workspace.deps, ctx, genes, n_gates, out_log, can_lazy_offboard);
            gates_routed += workspace.exec_gates.size();
            std::fill(workspace.decay.begin(), workspace.decay.end(), 1.0);
            continue;
        }

        workspace.ext_layer.clear();
        std::fill(workspace.in_ext.begin(), workspace.in_ext.end(), false);
        for (int gid : workspace.front_layer) {
            for (int nxt : dag.nodes[gid].next_gates) if (!workspace.in_ext[nxt]) { workspace.in_ext[nxt] = true; workspace.ext_layer.push_back(nxt); }
        }

        if (!de_route_step<LOGGING>(dag, workspace.front_layer, workspace.ext_layer, workspace.decay, ctx, genes, n_gates, out_log)) {
            break;
        }
    }

    for (size_t p = 0; p < p_cnt; ++p) { int temp_p = p; de_force_down<LOGGING>(temp_p, -1, ctx, genes, n_gates, out_log, false); }
    if constexpr (LOGGING) { if (out_log && !out_log->empty()) std::stable_sort(out_log->begin(), out_log->end(), [](const RoutingRecord& a, const RoutingRecord& b) { return a.start_time < b.start_time; }); }

    for (double t : ctx.q_time) max_physical_time = fast_max_hyb(max_physical_time, t);
    for (int d : ctx.q_depth) max_logical_depth = fast_max_hyb(max_logical_depth, d);

    out_swaps = ctx.out_swaps; out_moves = ctx.out_moves;
    return (ctx.total_routing_cost * 1000.0) + (out_swaps * 5000.0) + (max_physical_time * 1.0) + (max_logical_depth * 0.5);
}

double SabreDEHybridRouter::evaluate_fitness(const double* genes, int& out_swaps, int& out_moves, double& max_physical_time, int& max_logical_depth, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard) {
    if (out_log) return evaluate_fitness_impl<true>(genes, out_swaps, out_moves, max_physical_time, max_logical_depth, out_log, can_lazy_offboard);
    return evaluate_fitness_impl<false>(genes, out_swaps, out_moves, max_physical_time, max_logical_depth, nullptr, can_lazy_offboard);
}

// =========================================================================
//                          HYBRID ALGORITHMS
// =========================================================================

void SabreDEHybridRouter::initialize_de_population(int pop_size, int n_genes, int n_gates, int blocks, const std::vector<int>& logical_exec_order, std::vector<double>& pop, std::vector<double>& fitness, boost::random::xoshiro256pp& rng, bool can_lazy) {
    boost::random::uniform_real_distribution<double> dist_real(0.0, 1.0);
    std::vector<double> sabre_genes(n_genes, 0.5);

    for (size_t i = 0; i < logical_exec_order.size(); ++i) {
        int gid = logical_exec_order[i];
        sabre_genes[gid] = 1.0 - (static_cast<double>(i) / logical_exec_order.size());
        for (int b = 1; b < blocks; ++b) sabre_genes[b * n_gates + gid] = dist_real(rng);
    }

    int d_s, d_m, d_d; double d_t;
    std::vector<double> trial(n_genes);
    for (int i = 0; i < pop_size; ++i) {
        for (int j = 0; j < n_genes; ++j) trial[j] = pop[i * n_genes + j] = (i == 0) ? sabre_genes[j] : dist_real(rng);

        // Zero-copy: przekazanie wskaźnika
        fitness[i] = evaluate_fitness(trial.data(), d_s, d_m, d_t, d_d, nullptr, can_lazy);

        if (i == 0) {
            std::cout << "\n=== OCENA ZIARNA SABRE PRZEZ SILNIK DE (Rozwiązanie Bazowe) ===\n";
            std::cout << "Wygenerowane operacje SWAP przez DE: " << d_s << "\n";
            std::cout << "Wygenerowane operacje MOVE przez DE: " << d_m << "\n";
            std::cout << "Czas fizyczny: " << d_t << " [jednostek] | Głębokość logiczna: " << d_d << " [warstw]\n";
            std::cout << "Koszt w silniku DE (z karą za SWAP): " << fitness[i] << "\n";
            std::cout << "===============================================================\n\n";
        }
    }
}

void SabreDEHybridRouter::optimize(int pop_size, int generations, bool allow_lazy_offboard, double f, double cr) {
    std::cout << "\n[Hybryda 3.0] Faza 1: Generowanie ziarna (Priority Seed) przez heurystykę SABRE...\n";
    std::vector<int> logical_exec_order = prepare_sabre_seed();
    int n_gates = dag.nodes.size();
    if (n_gates == 0) return;

    bool can_lazy = allow_lazy_offboard && (num_qubits < topo.qubits.size());
    if (can_lazy) std::cout << "Aktywowano Lazy Off-boarding (wolne kubity fizyczne dostępne).\n";

    int blocks = can_lazy ? 4 : 3;
    int n_genes = n_gates * blocks;

    std::cout << "[Hybryda 3.0] Odrzucono sztuczne SWAP-y. Wyciągnięto czysty harmonogram logiczny (" << logical_exec_order.size() << " bramek).\n";
    std::cout << "[Hybryda 3.0] Faza 2: Dynamiczne DE minimalizuje koszty na bazie ziarna...\n";
    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<double> pop(pop_size * n_genes), fitness(pop_size);
    boost::random::xoshiro256pp rng(1337);
    boost::random::uniform_real_distribution<double> dist_real(0.0, 1.0);
    boost::random::uniform_int_distribution<int> dist_int(0, pop_size - 1);

    initialize_de_population(pop_size, n_genes, n_gates, blocks, logical_exec_order, pop, fitness, rng, can_lazy);

    std::vector<double> trial(n_genes);
    int d_s, d_m, d_d; double d_t;

    for (int gen = 0; gen < generations; ++gen) {
        for (int i = 0; i < pop_size; ++i) {
            int a, b, c;
            do { a = dist_int(rng); } while (a == i);
            do { b = dist_int(rng); } while (b == i || b == a);
            do { c = dist_int(rng); } while (c == i || c == a || c == b);

            int R = boost::random::uniform_int_distribution<int>(0, n_genes - 1)(rng);
            for (int j = 0; j < n_genes; ++j) {
                bool mutate = (dist_real(rng) < cr || j == R);
                trial[j] = mutate ? std::clamp(pop[a * n_genes + j] + f * (pop[b * n_genes + j] - pop[c * n_genes + j]), 0.0, 1.0) : pop[i * n_genes + j];
            }

            double trial_fit = evaluate_fitness(trial.data(), d_s, d_m, d_t, d_d, nullptr, can_lazy);
            if (trial_fit < fitness[i]) {
                std::copy(trial.begin(), trial.end(), pop.begin() + i * n_genes);
                fitness[i] = trial_fit;
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto best_it = std::min_element(fitness.begin(), fitness.end());
    int best_idx = std::distance(fitness.begin(), best_it);

    std::vector<RoutingRecord> final_log;
    int final_swaps = 0, final_moves = 0, final_depth = 0; double final_makespan = 0.0;

    // Przekazanie adresu najlepszego genotypu z bufora wektora, unikanie niepotrzebnej kopii
    double final_score = evaluate_fitness(pop.data() + best_idx * n_genes, final_swaps, final_moves, final_makespan, final_depth, &final_log, can_lazy);

    std::cout << "\n=== WALIDACJA OBWODU LOGICZNEGO ===\n";
    if (!validate_logical_circuit(final_log)) exit(-1);
    std::cout << "===================================\n";

    print_optimized_circuit(final_log);
    print_statistics(final_score, final_makespan, final_depth, final_swaps, final_moves, std::chrono::duration<double>(end_time - start_time).count());
    export_moves_to_file(final_log, "../hybrid_schedule_moves_" + name + ".txt");
}

void SabreDEHybridRouter::optimize_windowed(int pop_size, int generations_per_window, int window_size, double accept_ratio, bool allow_lazy_offboard, double f, double cr) {
    std::cout << "\n[Hybryda 3.0 (WINDOWED)] Faza 1: Generowanie ziarna (Priority Seed) przez heurystykę SABRE...\n";
    std::vector<int> logical_exec_order = prepare_sabre_seed();
    int n_gates = dag.nodes.size();
    if (n_gates == 0) return;

    bool can_lazy = allow_lazy_offboard && (num_qubits < topo.qubits.size());
    if (can_lazy) std::cout << "Aktywowano Lazy Off-boarding (wolne kubity fizyczne dostępne).\n";

    int blocks = can_lazy ? 4 : 3;
    int n_genes = n_gates * blocks;

    std::cout << "[Hybryda 3.0] Faza 2: DE Windowed na bazie logicznej (Rozmiar okna: " << window_size << ")...\n";
    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<double> pop(pop_size * n_genes), fitness(pop_size);
    boost::random::xoshiro256pp rng(1337);

    initialize_de_population(pop_size, n_genes, n_gates, blocks, logical_exec_order, pop, fitness, rng, can_lazy);

    int step = std::max(1, static_cast<int>(window_size * accept_ratio));
    int d_s, d_m, d_d; double d_t;

    for (int w_start = 0; w_start < n_gates; w_start += step) {
        int w_end = std::min(w_start + window_size, n_gates);
        std::cout << ">>> Przetwarzanie okna: operacje [" << w_start << " - " << w_end - 1 << "] / " << n_gates << "\n";
        run_evolution_loop_windowed(pop, fitness, pop_size, generations_per_window, n_genes, n_gates, f, cr, rng, w_start, w_end);

        auto best_iter = std::min_element(fitness.begin(), fitness.end());
        int best_idx = std::distance(fitness.begin(), best_iter), accept_end = std::min(w_start + step, n_gates);

        for (int i = 0; i < pop_size; ++i) {
            if (i == best_idx) continue;
            for (int blk = 0; blk < blocks; ++blk) {
                for (int g = w_start; g < accept_end; ++g) pop[i * n_genes + blk * n_gates + g] = pop[best_idx * n_genes + blk * n_gates + g];
            }

            // Zero-copy: brak niepotrzebnego wektora ind, przekazujemy offset z surowej pamięci pop
            fitness[i] = evaluate_fitness(pop.data() + i * n_genes, d_s, d_m, d_t, d_d, nullptr, can_lazy);
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto best_it = std::min_element(fitness.begin(), fitness.end());
    int best_idx = std::distance(fitness.begin(), best_it);

    std::vector<RoutingRecord> final_log;
    int final_swaps = 0, final_moves = 0, final_depth = 0; double final_makespan = 0.0;

    // Podobnie tutaj – pracujemy na offsecie
    double final_score = evaluate_fitness(pop.data() + best_idx * n_genes, final_swaps, final_moves, final_makespan, final_depth, &final_log, can_lazy);

    std::cout << "\n=== WALIDACJA OBWODU LOGICZNEGO ===\n";
    if (!validate_logical_circuit(final_log)) exit(-1);
    std::cout << "===================================\n";

    print_optimized_circuit(final_log);
    print_statistics(final_score, final_makespan, final_depth, final_swaps, final_moves, std::chrono::duration<double>(end_time - start_time).count());
    export_moves_to_file(final_log, "../hybrid_schedule_moves_windowed_" + name + ".txt");
}

void SabreDEHybridRouter::run_evolution_loop_windowed(std::vector<double>& pop, std::vector<double>& fitness, int pop_size, int generations, int n_genes, int n_gates, double f, double cr, boost::random::xoshiro256pp& rng, int w_start, int w_end) {
    boost::random::uniform_int_distribution<int> dist_int(0, pop_size - 1);
    std::vector<double> trial(n_genes);
    int d_s, d_m, d_d; double d_t; bool can_lazy = (n_genes / n_gates == 4);

    for (int gen = 0; gen < generations; ++gen) {
        for (int i = 0; i < pop_size; ++i) {
            int a, b, c; do { a = dist_int(rng); } while (a == i); do { b = dist_int(rng); } while (b == i || b == a); do { c = dist_int(rng); } while (c == i || c == a || c == b);
            generate_trial_vector_windowed(pop, trial, i, a, b, c, n_genes, n_gates, f, cr, rng, w_start, w_end);
            double trial_fit = evaluate_fitness(trial.data(), d_s, d_m, d_t, d_d, nullptr, can_lazy);
            if (trial_fit < fitness[i]) { std::copy(trial.begin(), trial.end(), pop.begin() + i * n_genes); fitness[i] = trial_fit; }
        }
    }
}

void SabreDEHybridRouter::generate_trial_vector_windowed(const std::vector<double>& pop, std::vector<double>& trial, int current_idx, int a, int b, int c, int n_genes, int n_gates, double f, double cr, boost::random::xoshiro256pp& rng, int w_start, int w_end) {
    boost::random::uniform_real_distribution<double> dist_real(0.0, 1.0);
    int R_idx = boost::random::uniform_int_distribution<int>(0, (w_end - w_start) * (n_genes / n_gates) - 1)(rng);
    int R = (R_idx / (w_end - w_start)) * n_gates + w_start + (R_idx % (w_end - w_start));

    const double *p_a = pop.data() + a * n_genes, *p_b = pop.data() + b * n_genes, *p_c = pop.data() + c * n_genes, *p_curr = pop.data() + current_idx * n_genes;
    for (int j = 0; j < n_genes; ++j) {
        int gate_id = j % n_gates;
        trial[j] = (gate_id >= w_start && gate_id < w_end && (dist_real(rng) < cr || j == R)) ? std::clamp(p_a[j] + f * (p_b[j] - p_c[j]), 0.0, 1.0) : p_curr[j];
    }
}


// =========================================================================
//                              SABRE SEED
// =========================================================================

void SabreDEHybridRouter::sabre_force_qubit_down(int p, int current_gate_id, SabreRunContext& ctx) {
    if (ctx.state.q_on_resonator[p] != -1) {
        int r = ctx.state.q_on_resonator[p];
        double step_start = fast_max_hyb(ctx.state.q_ready_time[p], ctx.state.r_ready_time[r]);

        ctx.out_moves++;
        ctx.total_routing_cost += MOVE_COST;
        ctx.state.q_ready_time[p] = step_start + MOVE_COST;
        ctx.state.r_ready_time[r] = step_start + MOVE_COST;
        ctx.q_logical_depth[p] += 1;
        ctx.state.q_on_resonator[p] = -1;
        ctx.state.res_occupant[r] = -1;
    }
}

void SabreDEHybridRouter::sabre_execute_gates(const HybridSabreDAG& current_dag, std::vector<int>& executable_gates, std::vector<int>& front_layer, std::vector<int>& deps, SabreRunContext& ctx, std::vector<int>* out_logical_order) {
    for (int gate_id : executable_gates) {
        if (out_logical_order) out_logical_order->push_back(gate_id);

        const auto& node = current_dag.nodes[gate_id];
        int p1 = ctx.state.logical_to_phys_qubit[node.lq1];

        if (node.lq2 == -1) {
            if (node.q1_1q_count > 0) sabre_force_qubit_down(p1, node.id, ctx);
            ctx.state.q_ready_time[p1] += node.q1_1q_count;
            ctx.q_logical_depth[p1] += node.q1_1q_count;
        } else {
            int p2 = ctx.state.logical_to_phys_qubit[node.lq2];
            if (node.q1_1q_count > 0) sabre_force_qubit_down(p1, node.id, ctx);
            if (node.q2_1q_count > 0) sabre_force_qubit_down(p2, node.id, ctx);

            ctx.state.q_ready_time[p1] += node.q1_1q_count;
            ctx.state.q_ready_time[p2] += node.q2_1q_count;
            ctx.q_logical_depth[p1] += node.q1_1q_count;
            ctx.q_logical_depth[p2] += node.q2_1q_count;

            int chosen_res = -1;
            const auto& p2_res = topo.qubits[p2].connected_resonators;
            for (int r : topo.qubits[p1].connected_resonators) {
                if (std::find(p2_res.begin(), p2_res.end(), r) != p2_res.end()) { chosen_res = r; break; }
            }

            if (chosen_res != -1) {
                double start_time = fast_max_hyb(ctx.state.q_ready_time[p1], ctx.state.q_ready_time[p2]);
                bool q1_goes_in = (ctx.state.q_ready_time[p1] <= ctx.state.q_ready_time[p2]);
                int p_in = q1_goes_in ? p1 : p2;
                int p_stay = q1_goes_in ? p2 : p1;

                sabre_force_qubit_down(p_stay, node.id, ctx);
                if (ctx.state.res_occupant[chosen_res] != -1 && ctx.state.res_occupant[chosen_res] != p_in) {
                    sabre_force_qubit_down(ctx.state.res_occupant[chosen_res], node.id, ctx);
                }

                if (ctx.state.q_on_resonator[p_in] != chosen_res) {
                    sabre_force_qubit_down(p_in, node.id, ctx);
                    start_time = fast_max_hyb(ctx.state.q_ready_time[p_in], fast_max_hyb(ctx.state.q_ready_time[p_stay], ctx.state.r_ready_time[chosen_res]));
                    ctx.out_moves++;
                    ctx.total_routing_cost += MOVE_COST;
                    start_time += MOVE_COST;
                    ctx.state.q_on_resonator[p_in] = chosen_res;
                    ctx.state.res_occupant[chosen_res] = p_in;
                    ctx.q_logical_depth[p_in] += 1;
                } else {
                    start_time = fast_max_hyb(ctx.state.q_ready_time[p_in], fast_max_hyb(ctx.state.q_ready_time[p_stay], ctx.state.r_ready_time[chosen_res]));
                }

                ctx.state.q_ready_time[p_in] = start_time + 1.0;
                ctx.state.q_ready_time[p_stay] = start_time + 1.0;
                ctx.state.r_ready_time[chosen_res] = start_time + 1.0;
                int gate_layer = fast_max_hyb(ctx.q_logical_depth[p_in], ctx.q_logical_depth[p_stay]);
                ctx.q_logical_depth[p_in] = gate_layer + 1;
                ctx.q_logical_depth[p_stay] = gate_layer + 1;
            }
        }
        for (int next_id : node.next_gates) {
            if (--deps[next_id] == 0) front_layer.push_back(next_id);
        }
    }
}

bool SabreDEHybridRouter::sabre_route_step(const HybridSabreDAG& current_dag, const std::vector<int>& front_layer, const std::vector<int>& extended_layer, SabreRunContext& ctx) {
    std::vector<bool> in_active(topo.qubits.size(), false);
    std::vector<int> active_phys_qubits;
    active_phys_qubits.reserve(num_qubits);

    for (int gate_id : front_layer) {
        const auto& node = current_dag.nodes[gate_id];
        if (node.lq1 != -1) {
            int p1 = ctx.state.logical_to_phys_qubit[node.lq1];
            if (!in_active[p1]) { in_active[p1] = true; active_phys_qubits.push_back(p1); }
        }
        if (node.lq2 != -1) {
            int p2 = ctx.state.logical_to_phys_qubit[node.lq2];
            if (!in_active[p2]) { in_active[p2] = true; active_phys_qubits.push_back(p2); }
        }
    }

    struct SwapCandidate { int p1, p2, r_shared; double score; };
    std::vector<SwapCandidate> candidates;
    candidates.reserve(topo.qubits.size() * 4);

    for (int p : active_phys_qubits) {
        for (int r : topo.qubits[p].connected_resonators) {
            for (int neighbor_p : topo.resonators[r].connected_qubits) {
                if (p == neighbor_p) continue;
                int u = std::min(p, neighbor_p);
                int v = std::max(p, neighbor_p);

                bool already_added = false;
                for (const auto& c : candidates) {
                    if (c.p1 == u && c.p2 == v) { already_added = true; break; }
                }
                if (already_added) continue;

                double score_F = 0.0, score_E = 0.0;
                for (int gate_id : front_layer) {
                    const auto& n = current_dag.nodes[gate_id];
                    if (n.lq2 != -1) {
                        int lq1_p = ctx.state.logical_to_phys_qubit[n.lq1], lq2_p = ctx.state.logical_to_phys_qubit[n.lq2];
                        if (lq1_p == p) lq1_p = neighbor_p; else if (lq1_p == neighbor_p) lq1_p = p;
                        if (lq2_p == p) lq2_p = neighbor_p; else if (lq2_p == neighbor_p) lq2_p = p;
                        score_F += get_distance(lq1_p, lq2_p);
                    }
                }
                for (int gate_id : extended_layer) {
                    const auto& n = current_dag.nodes[gate_id];
                    if (n.lq2 != -1) {
                        int lq1_p = ctx.state.logical_to_phys_qubit[n.lq1], lq2_p = ctx.state.logical_to_phys_qubit[n.lq2];
                        if (lq1_p == p) lq1_p = neighbor_p; else if (lq1_p == neighbor_p) lq1_p = p;
                        if (lq2_p == p) lq2_p = neighbor_p; else if (lq2_p == neighbor_p) lq2_p = p;
                        score_E += get_distance(lq1_p, lq2_p);
                    }
                }

                double f_norm = front_layer.empty() ? 0 : (score_F / front_layer.size());
                double e_norm = extended_layer.empty() ? 0 : (score_E / extended_layer.size());
                double score = (f_norm + EXTENDED_LAYER_WEIGHT * e_norm) * fast_max_hyb(ctx.decay[u], ctx.decay[v]);
                candidates.push_back({u, v, r, score});
            }
        }
    }

    SwapCandidate* best_candidate = nullptr;
    if (!candidates.empty()) {
        best_candidate = &candidates[0];
        for (size_t i = 1; i < candidates.size(); ++i) {
            if (candidates[i].score < best_candidate->score) best_candidate = &candidates[i];
        }
    }

    if (best_candidate != nullptr) {
        int p1 = best_candidate->p1, p2 = best_candidate->p2, r_step = best_candidate->r_shared;
        sabre_force_qubit_down(p1, -1, ctx);
        sabre_force_qubit_down(p2, -1, ctx);

        if (ctx.state.res_occupant[r_step] != -1) {
            sabre_force_qubit_down(ctx.state.res_occupant[r_step], -1, ctx);
        }

        double step_start = fast_max_hyb(ctx.state.q_ready_time[p1], fast_max_hyb(ctx.state.r_ready_time[r_step], ctx.state.q_ready_time[p2]));
        bool is_move = (ctx.state.phys_qubit_to_logical[p2] == -1) || (ctx.state.phys_qubit_to_logical[p1] == -1);
        double op_cost = is_move ? (2.0 * MOVE_COST) : SWAP_COST;

        int lq_1 = ctx.state.phys_qubit_to_logical[p1], lq_2 = ctx.state.phys_qubit_to_logical[p2];

        if (is_move) ctx.out_moves += 2; else ctx.out_swaps++;
        ctx.total_routing_cost += op_cost;

        std::swap(ctx.state.phys_qubit_to_logical[p1], ctx.state.phys_qubit_to_logical[p2]);
        if (lq_1 != -1) ctx.state.logical_to_phys_qubit[lq_1] = p2;
        if (lq_2 != -1) ctx.state.logical_to_phys_qubit[lq_2] = p1;

        int current_layer = fast_max_hyb(ctx.q_logical_depth[p1], ctx.q_logical_depth[p2]);
        ctx.q_logical_depth[p1] = current_layer + (is_move ? 2 : 1);
        ctx.q_logical_depth[p2] = current_layer + (is_move ? 2 : 1);

        ctx.state.q_ready_time[p1] = step_start + op_cost;
        ctx.state.r_ready_time[r_step] = step_start + op_cost;
        ctx.state.q_ready_time[p2] = step_start + op_cost;

        ctx.decay[p1] += 0.2; ctx.decay[p2] += 0.2;
        return true;
    }
    return false;
}

std::vector<int> SabreDEHybridRouter::run_sabre_pass(const HybridSabreDAG& current_dag, const std::vector<int>& initial_l2p, bool is_final, std::vector<int>* out_logical_order) {
    int n_gates = current_dag.nodes.size();
    SabreRunContext ctx;
    ctx.state.logical_to_phys_qubit = initial_l2p;
    ctx.state.phys_qubit_to_logical.assign(topo.qubits.size(), -1);

    for (int i = 0; i < num_qubits; ++i) {
        if (initial_l2p[i] != -1) ctx.state.phys_qubit_to_logical[initial_l2p[i]] = i;
    }

    ctx.state.q_on_resonator.assign(topo.qubits.size(), -1);
    ctx.state.res_occupant.assign(topo.resonators.size(), -1);
    ctx.state.q_ready_time.assign(topo.qubits.size(), 0.0);
    ctx.state.r_ready_time.assign(topo.resonators.size(), 0.0);
    ctx.q_logical_depth.assign(topo.qubits.size(), 0);
    ctx.decay.assign(topo.qubits.size(), 1.0);

    std::vector<int> front_layer;
    front_layer.reserve(num_qubits);
    std::vector<int> deps = current_dag.initial_deps;
    for (int i = 0; i < n_gates; ++i) {
        if (deps[i] == 0) front_layer.push_back(i);
    }

    int gates_routed = 0;
    std::vector<int> executable_gates, extended_layer;
    executable_gates.reserve(num_qubits);
    extended_layer.reserve(num_qubits * 4);
    std::vector<bool> in_extended(n_gates, false);

    int debug_iterations = 0;
    const int MAX_ITERATIONS = n_gates * 100;

    while (gates_routed < n_gates) {
        if (++debug_iterations > MAX_ITERATIONS) break;

        executable_gates.clear();
        for (size_t i = 0; i < front_layer.size(); ) {
            int gate_id = front_layer[i];
            const auto& node = current_dag.nodes[gate_id];
            bool is_executable = false;

            if (node.lq2 == -1) {
                is_executable = true;
            } else {
                int p1 = ctx.state.logical_to_phys_qubit[node.lq1];
                int p2 = ctx.state.logical_to_phys_qubit[node.lq2];
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
            sabre_execute_gates(current_dag, executable_gates, front_layer, deps, ctx, out_logical_order);
            gates_routed += executable_gates.size();
            std::fill(ctx.decay.begin(), ctx.decay.end(), 1.0);
            continue;
        }

        extended_layer.clear();
        std::fill(in_extended.begin(), in_extended.end(), false);
        for (int gate_id : front_layer) {
            for (int nxt : current_dag.nodes[gate_id].next_gates) {
                if (!in_extended[nxt]) { in_extended[nxt] = true; extended_layer.push_back(nxt); }
            }
        }

        sabre_route_step(current_dag, front_layer, extended_layer, ctx);
    }

    if (is_final) {
        for (size_t p = 0; p < topo.qubits.size(); ++p) sabre_force_qubit_down(p, -1, ctx);

        double max_physical_time = 0;
        int max_logical_depth = 0;
        for (double t : ctx.state.q_ready_time) max_physical_time = fast_max_hyb(max_physical_time, t);
        for (int d : ctx.q_logical_depth) max_logical_depth = fast_max_hyb(max_logical_depth, d);

        std::cout << "\n=== STATYSTYKI SAMEGO SABRE (Czyste wykonanie heurystyki) ===\n";
        std::cout << "Wykorzystane operacje SWAP: " << ctx.out_swaps << "\n";
        std::cout << "Wykorzystane operacje MOVE: " << ctx.out_moves << "\n";
        std::cout << "Czas fizyczny obwodu (Makespan): " << max_physical_time << " [jednostek]\n";
        std::cout << "Głębokość warstwowa (Logical Depth): " << max_logical_depth << " [warstw]\n";
        std::cout << "Koszt routingu (Total Routing Cost): " << ctx.total_routing_cost << "\n";
        std::cout << "==============================================================\n";
    }

    return ctx.state.logical_to_phys_qubit;
}

// =========================================================================
//                              DE SEED
// =========================================================================

// =========================================================================
//                          HYBRID ALGORITHMS
// =========================================================================

std::vector<int> SabreDEHybridRouter::prepare_sabre_seed() {
    std::vector<int> initial_mapping(num_qubits);
    saved_initial_l2p.assign(topo.qubits.size(), -1);
    for (int i = 0; i < num_qubits; ++i) { initial_mapping[i] = i; saved_initial_l2p[i] = i; }

    auto post_fwd = run_sabre_pass(dag, initial_mapping, false, nullptr);
    auto pre_final = run_sabre_pass(reversed_dag, post_fwd, false, nullptr);

    std::vector<int> logical_exec_order;
    run_sabre_pass(dag, pre_final, true, &logical_exec_order);
    return logical_exec_order;
}

// =========================================================================
//                               UTILS
// =========================================================================

bool SabreDEHybridRouter::validate_logical_circuit(const std::vector<RoutingRecord>& log) const {
    std::vector<std::vector<std::string>> expected(num_qubits);
    for (const auto& g : original_circuit) {
        if (g.type == "1Q") expected[g.target_q1].emplace_back("1Q");
        else {
            expected[g.target_q1].push_back("2Q_" + std::to_string(g.target_q2));
            expected[g.target_q2].push_back("2Q_" + std::to_string(g.target_q1));
        }
    }
    std::vector<std::vector<std::string>> actual(num_qubits);
    for (const auto& op : log) {
        if (op.op_type == OpType::GATE_1Q) actual[op.lq1].emplace_back("1Q");
        else if (op.op_type == OpType::GATE_2Q_LOCAL || op.op_type == OpType::GATE_2Q_GLOBAL) {
            actual[op.lq1].push_back("2Q_" + std::to_string(op.lq2));
            actual[op.lq2].push_back("2Q_" + std::to_string(op.lq1));
        }
    }
    bool is_valid = true;
    for (int i = 0; i < num_qubits; ++i) {
        if (expected[i] != actual[i]) {
            is_valid = false;
            std::cout << "[BŁĄD WALIDACJI] Niezgodność na kubicie logicznym lq" << i << "!\n";
        }
    }
    if (is_valid) std::cout << "[OK] Walidacja pomyślna. Czystość logiczna zachowana na 100% (DE Dynamic Routing)!\n";
    return is_valid;
}

void SabreDEHybridRouter::print_statistics(double best_score, double physical_time, int depth, int sw, int mv, double exec_time) {
    std::cout << "\n=== STATYSTYKI DE (HYBRYDA 3.0 - Zoptymalizowany wynik) ===\n";
    std::cout << "Wykorzystane operacje SWAP (Wygenerowane czysto przez DE): " << sw << "\n";
    std::cout << "Wykorzystane operacje MOVE: " << mv << "\n";
    std::cout << "Czas fizyczny: " << physical_time << " [jednostek] | Głębokość logiczna: " << depth << " [warstw]\n";
    std::cout << "Całkowity koszt funkcji celu: " << best_score << "\n";
    std::cout << "Czas przetwarzania algorytmu: " << exec_time << " sekund\n";
    std::cout << "============================================================\n\n";
}

void SabreDEHybridRouter::print_optimized_circuit(const std::vector<RoutingRecord>& log) const {
    if (log.empty()) { std::cout << "Zoptymalizowany obwód jest pusty!!\n"; return; }
    int n_q = topo.qubits.size();
    std::vector<std::stringstream> lines_q(n_q);
    for (int i = 0; i < n_q; ++i) { lines_q[i] << "q" << (i < 10 ? "0" : "") << i << (i < num_qubits ? " [L]: " : " [E]: "); }

    for (const auto& move : log) {
        int max_len = 0;
        for (const auto& ss : lines_q) if (ss.str().length() > max_len) max_len = ss.str().length();
        for (auto& ss : lines_q) for (int p = 0, pad = max_len - ss.str().length(); p < pad; ++p) ss << "-";

        int l1 = move.loc1, l2 = move.loc2;
        if (move.op_type == OpType::GATE_1Q) {
            if (l1 >= 0 && l1 < n_q) lines_q[l1] << "-[ U3 ]-";
            for (int i = 0; i < n_q; ++i) if (i != l1) lines_q[i] << "--------";
        } else if (move.op_type == OpType::SWAP) {
            int min_q = std::min(l1, l2), max_q = fast_max_hyb(l1, l2);
            for (int i = 0; i < n_q; ++i) lines_q[i] << ((i == l1 || i == l2) ? "--x--" : (i > min_q && i < max_q) ? "--|--" : "-----");
        } else if (move.op_type == OpType::MOVE) {
            if (l1 < RES_OFFSET && l2 < RES_OFFSET && l2 >= 0) {
                int min_q = std::min(l1, l2), max_q = fast_max_hyb(l1, l2);
                for (int i = 0; i < n_q; ++i) lines_q[i] << ((i == l1 || i == l2) ? "--m--" : (i > min_q && i < max_q) ? "--|--" : "-----");
            } else {
                int q_idx = (l1 < RES_OFFSET) ? l1 : l2;
                for (int i = 0; i < n_q; ++i) lines_q[i] << (i == q_idx ? "--m--" : "-----");
            }
        } else if (move.op_type == OpType::GATE_2Q_LOCAL || move.op_type == OpType::GATE_2Q_GLOBAL) {
            int q_stay = (l1 >= RES_OFFSET) ? l2 : l1, q_in = move.res_used, min_q = std::min(q_in, q_stay), max_q = fast_max_hyb(q_in, q_stay);
            bool q_in_is_control = false;
            for (const auto& g : original_circuit) if (g.target_q1 == move.lq1 || g.target_q2 == move.lq1) { q_in_is_control = (move.lq1 == g.target_q1); break; }
            std::string ctrl_char = (move.op_type == OpType::GATE_2Q_LOCAL) ? "--o--" : "--*--";
            std::string targ_char = (move.op_type == OpType::GATE_2Q_LOCAL) ? "-(+)-" : "-[+]-";
            for (int i = 0; i < n_q; ++i) lines_q[i] << (i == q_in ? (q_in_is_control ? ctrl_char : targ_char) : i == q_stay ? (q_in_is_control ? targ_char : ctrl_char) : (i > min_q && i < max_q && q_in != -1) ? "--|--" : "-----");
        }
    }
    std::cout << "\n=== WIZUALIZACJA ZOPTYMALIZOWANEGO OBWODU (HYBRYDA 3.0) ===\n";
    for (int i = 0; i < n_q; ++i) std::cout << lines_q[i].str() << "\n";
    std::cout << "===========================================================\n";
}

void SabreDEHybridRouter::export_moves_to_file(const std::vector<RoutingRecord>& moves, const std::string& filename) {
    if (moves.empty()) return;
    std::ofstream out(filename);
    out << "time_step,op_type,gate_id,lq1,lq2,loc1,loc2,res_used\n";
    int time_step = -1; double current_time = -1.0;
    for (const auto& move : moves) {
        if (std::abs(move.start_time - current_time) > 1e-5) { time_step++; current_time = move.start_time; }
        out << time_step << "," << to_string(move.op_type) << "," << move.gate_id << "," << move.lq1 << "," << move.lq2 << "," << move.loc1 << "," << move.loc2 << "," << move.res_used << "\n";
    }
    std::cout << ">>> Zapisano harmonogram (Hybryda 3.0) do: " << filename << "\n";
}