#include "SabreDEHybridRouter.h"
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include <fstream>
#include <numeric>
#include <sstream>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <unordered_map>
#include <unordered_set>

inline double fast_max_hyb(double a, double b) { return (a > b) ? a : b; }

SabreDEHybridRouter::SabreDEHybridRouter(const BaseTopology& topology, const std::vector<Gate>& circuit, int qubits_cnt, std::string name, bool slow, int max_look_ahead_depth)
    : topo(topology), original_circuit(circuit), num_qubits(qubits_cnt), name(std::move(name)), slow(slow), max_look_ahead_depth(max_look_ahead_depth)
{
    init_distance_cache();

    // Budowa grafu skierowanego (DAG) dla oryginalnego i odwróconego obwodu
    dag = build_dag(original_circuit);
    std::vector<Gate> reversed_circuit = original_circuit;
    std::reverse(reversed_circuit.begin(), reversed_circuit.end());
    reversed_dag = build_dag(reversed_circuit);
    p_cnt_cache = topo.qubits.size();
    precompute_physical_swaps();
    allocate_workspace(dag.nodes.size());
}

inline void SabreDEHybridRouter::init_distance_cache() {
    int p_cnt = topo.qubits.size();
    dist_cache.assign(p_cnt * p_cnt, -1);

    std::cout << "[Inicjalizacja] Pre-kalkulacja pełnej macierzy odległości...\n";
    for (int i = 0; i < p_cnt; ++i) {
        for (int j = i; j < p_cnt; ++j) {
            if (i == j) {
                dist_cache[i * p_cnt + j] = 0;
            } else {
                int dist = topo.get_routing_path(i, j).resonators.size();
                dist_cache[i * p_cnt + j] = dist;
                dist_cache[j * p_cnt + i] = dist; // Wykorzystanie symetrii w czasie O(1)
            }
        }
    }
}

inline void SabreDEHybridRouter::precompute_physical_swaps() {
    int p_cnt = topo.qubits.size();
    precomputed_swaps_offsets.clear();
    precomputed_swaps_offsets.reserve(p_cnt + 1);
    flat_precomputed_swaps.reserve(p_cnt * 6); // Szacunkowy bezpieczny zapas

    precomputed_swaps_offsets.push_back(0);
    for (int p = 0; p < p_cnt; ++p) {
        for (int r : topo.qubits[p].connected_resonators) {
            for (int neighbor_p : topo.resonators[r].connected_qubits) {
                if (p != neighbor_p) {
                    flat_precomputed_swaps.push_back({neighbor_p, r});
                }
            }
        }
        precomputed_swaps_offsets.push_back(flat_precomputed_swaps.size());
    }
}

inline void SabreDEHybridRouter::allocate_workspace(int g_cnt) {
    int p_cnt = topo.qubits.size();
    int r_cnt = topo.resonators.size();

    workspace.l2p.resize(p_cnt, -1);
    workspace.p2l.resize(p_cnt, -1);
    workspace.q_on_res.resize(p_cnt, -1);
    workspace.res_occ.resize(r_cnt, -1);
    workspace.q_depth.resize(p_cnt, 0);
    workspace.q_time.resize(p_cnt, 0.0);
    workspace.r_time.resize(r_cnt, 0.0);

    workspace.front_layer.reserve(g_cnt);
    workspace.exec_gates.reserve(p_cnt);
    workspace.ext_layer.reserve(g_cnt);

    workspace.visited_epoch.assign(g_cnt, 0);
    workspace.current_bfs_epoch = 0;

    workspace.decay.resize(p_cnt, 1.0);
    workspace.deps.resize(g_cnt, 0);

    workspace.shared_res.reserve(r_cnt);
    workspace.in_active.resize(p_cnt, 0);
    workspace.active_phys.reserve(p_cnt);

    workspace.interaction_counts.assign(p_cnt, 0);
    workspace.inter_targets.resize(p_cnt * HybridWorkspace::MAX_INTER_PER_Q, 0);
    workspace.inter_weights.resize(p_cnt * HybridWorkspace::MAX_INTER_PER_Q, 0.0);

    workspace.route_step_counter = 0;
    workspace.seen_edges.assign(p_cnt * p_cnt, 0);

    workspace.f_locs.reserve(g_cnt);
    workspace.lookahead_offsets.reserve(5);
    workspace.lookahead_locs_offsets.reserve(5);

    workspace.flat_lookahead_layers.reserve(g_cnt);
    workspace.flat_lookahead_locs.reserve(g_cnt);

    workspace.phys_gene_weights.resize(p_cnt, 0.0);
}

// ==============================================================================
// OBLICZANIE ODLEGŁOŚCI Z CACHE
// ==============================================================================
int SabreDEHybridRouter::get_distance_old(int p1, int p2) const {
    if (p1 == p2) return 0; // Odległość do samego siebie to zawsze 0
    int p_cnt = topo.qubits.size();

    // Używamy spłaszczonej tablicy 1D jako macierzy 2D.
    // Wzór na indeks to: wiersz * szerokość + kolumna
    int idx = p1 * p_cnt + p2;

    // Jeśli wartość nie jest -1, to znaczy, że już to kiedyś liczyliśmy.
    // O(1) zwrot z bufora - oszczędza mnóstwo czasu!
    if (dist_cache[idx] != -1) return dist_cache[idx];

    // Jeśli liczymy to pierwszy raz: odpytujemy topologię o ścieżkę
    // i liczymy ile po drodze jest rezonatorów (to nasza odległość).
    int dist = topo.get_routing_path(p1, p2).resonators.size();

    // Zapisujemy wynik do cache SYMETRYCZNIE.
    // Odległość z p1 do p2 jest taka sama jak z p2 do p1.
    dist_cache[idx] = dist;
    dist_cache[p2 * p_cnt + p1] = dist;

    return dist;
}

int SabreDEHybridRouter::get_distance(int p1, int p2) const {
    return dist_cache[p1 * p_cnt_cache + p2];
}

// ==============================================================================
// BUDOWANIE SKIEROWANEGO GRAFU ACYKLICZNEGO (DAG)
// ==============================================================================
HybridSabreDAG SabreDEHybridRouter::build_dag(const std::vector<Gate>& circuit) const {
    HybridSabreDAG local_dag;
    // [OPTYMALIZACJA] Pre-alokacja grafu, zakładając pesymistycznie 1 bramkę = 1 węzeł.
    // Zapobiega to realokacjom i kopiowaniu obiektów HybridGateNode podczas dodawania do wektora.
    local_dag.nodes.reserve(circuit.size());

    std::vector<int> last_gate(num_qubits, -1);
    std::vector<absl::InlinedVector<int, 4>> pending_1q(num_qubits);
    int group_id = 0;

    for (const auto& g : circuit) {
        if (g.type == "1Q") {
            pending_1q[g.target_q1].push_back(g.id);
        } else {
            add_2q_node_to_dag(local_dag, g, pending_1q, last_gate, group_id);
        }
    }

    // Zamknięcie grafu dla wiszących na końcu bramek 1-kubitowych
    flush_trailing_1q_nodes(local_dag, pending_1q, last_gate, group_id);

    // Finalizacja wag wejściowych
    compute_initial_deps(local_dag);

    return local_dag;
}

inline void SabreDEHybridRouter::add_2q_node_to_dag(
    HybridSabreDAG& local_dag,
    const Gate& g,
    std::vector<absl::InlinedVector<int, 4>>& pending_1q,
    std::vector<int>& last_gate,
    int& group_id) const
{
    auto& p1_pending = pending_1q[g.target_q1];
    auto& p2_pending = pending_1q[g.target_q2];

    HybridGateNode node(
        group_id++, g.target_q1, g.target_q2,
        p1_pending.size(), p2_pending.size(),
        p1_pending, p2_pending, // Kopiujemy wektory celowo...
        g.id
    );

    // ... aby zachować ich capacity i uniknąć ponownych alokacji w przyszłych operacjach push_back
    p1_pending.clear();
    p2_pending.clear();

    if (last_gate[node.lq1] != -1) {
        local_dag.nodes[last_gate[node.lq1]].next_gates.push_back(node.id);
        node.deps_count++;
    }

    if (last_gate[node.lq2] != -1 && last_gate[node.lq2] != last_gate[node.lq1]) {
        local_dag.nodes[last_gate[node.lq2]].next_gates.push_back(node.id);
        node.deps_count++;
    }

    local_dag.nodes.push_back(std::move(node)); // std::move wymusza użycie konstruktora przenoszącego (szybsze wstawianie)
    last_gate[g.target_q1] = local_dag.nodes.back().id;
    last_gate[g.target_q2] = local_dag.nodes.back().id;
}

inline void SabreDEHybridRouter::flush_trailing_1q_nodes(
    HybridSabreDAG& local_dag,
    std::vector<absl::InlinedVector<int, 4>>& pending_1q,
    std::vector<int>& last_gate,
    int& group_id) const
{
    for (int i = 0; i < num_qubits; ++i) {
        if (!pending_1q[i].empty()) {
            HybridGateNode node(
                group_id++, i, -1,
                pending_1q[i].size(), 0,
                pending_1q[i], {}, -1 // <--- ZMIANA: używamy {} zamiast std::vector<int>{}
            );

            if (last_gate[i] != -1) {
                local_dag.nodes[last_gate[i]].next_gates.push_back(node.id);
                node.deps_count++;
            }

            local_dag.nodes.push_back(std::move(node));
            last_gate[i] = local_dag.nodes.back().id;

            pending_1q[i].clear();
        }
    }
}

inline void SabreDEHybridRouter::compute_initial_deps(HybridSabreDAG& local_dag) const {
    local_dag.initial_deps.resize(local_dag.nodes.size());
    for (size_t i = 0; i < local_dag.nodes.size(); ++i) {
        local_dag.initial_deps[i] = local_dag.nodes[i].deps_count;
    }
}

// ==============================================================================
// Ściąganie kubitu z rezonatora z DE (podobnie jak w sabre - funkcja zoptymalizowana)
// ==============================================================================
inline int SabreDEHybridRouter::de_select_lazy_target(int p, int r, int gate_id, const DERunContext& ctx, const double* genes, int n_gates) const {
    int valid_count = 1;
    const auto& conn = topo.resonators[r].connected_qubits;

    for (int c : conn) {
        if (c != p && ctx.p2l[c] == -1 && ctx.q_on_res[c] == -1) valid_count++;
    }

    if (valid_count > 1) {
        int gene_idx = std::min(static_cast<int>(genes[n_gates * 3 + gate_id] * valid_count), valid_count - 1);
        if (gene_idx > 0) {
            int c_idx = 1;
            for (int c : conn) {
                if (c != p && ctx.p2l[c] == -1 && ctx.q_on_res[c] == -1) {
                    if (c_idx == gene_idx) return c;
                    c_idx++;
                }
            }
        }
    }
    return p;
}

template <bool LOGGING, bool CAN_LAZY>
void SabreDEHybridRouter::de_force_down(int& p, int gate_id, DERunContext& ctx, const double* genes, int n_gates, std::vector<RoutingRecord>* out_log) {
    int r = ctx.q_on_res[p];
    if (r == -1) return; // Szybkie wyjście

    int target_p = p;
    // Kompilator całkowicie usunie ten blok, jeśli CAN_LAZY == false (zero overhead)
    if constexpr (CAN_LAZY) {
        if (gate_id >= 0) {
            target_p = de_select_lazy_target(p, r, gate_id, ctx, genes, n_gates);
        }
    }

    double start = fast_max_hyb(ctx.q_time[p], ctx.r_time[r]);
    if (target_p != p) start = fast_max_hyb(start, ctx.q_time[target_p]);

    if constexpr (LOGGING) {
        out_log->push_back({OpType::MOVE, gate_id, ctx.p2l[p], -1, r + RES_OFFSET, target_p, -1, start});
    }

    ctx.out_moves++;
    ctx.total_routing_cost += MOVE_COST;

    ctx.q_time[target_p] = start + MOVE_COST;
    ctx.r_time[r] = start + MOVE_COST;

    if (target_p != p) {
        ctx.p2l[target_p] = ctx.p2l[p];
        ctx.l2p[ctx.p2l[p]] = target_p;
        ctx.p2l[p] = -1;
        ctx.q_depth[target_p] = ctx.q_depth[p] + 1;
        ctx.q_depth[p] = 0;
    } else {
        ctx.q_depth[p] += 1;
    }

    ctx.q_on_res[p] = -1;
    ctx.q_on_res[target_p] = -1;
    ctx.res_occ[r] = -1;
    p = target_p;
}

// ==============================================================================
// Wykonywanie bramek w zależności od piorytetu wygenerowanego przez DE
// ==============================================================================
template <bool LOGGING>
void SabreDEHybridRouter::de_execute_1q_logic(const HybridGateNode& node, int p1, DERunContext& ctx, const double* genes, int n_gates, std::vector<RoutingRecord>* out_log, bool can_lazy) {
    if (node.q1_1q_count > 0) {
        de_force_down<LOGGING, false>(p1, node.id, ctx, genes, n_gates, out_log);
        if constexpr (LOGGING) {
            for (int k = 0; k < node.q1_1q_count; ++k)
                out_log->push_back({OpType::GATE_1Q, node.q1_1q_ids[k], node.lq1, -1, p1, -1, -1, ctx.q_time[p1] + k});
        }
        ctx.q_time[p1] += node.q1_1q_count;
        ctx.q_depth[p1] += node.q1_1q_count;
    }
}

template <bool LOGGING>
bool SabreDEHybridRouter::de_execute_2q_logic(const HybridGateNode& node, int p1, int p2, DERunContext& ctx, const double* genes, int n_gates, std::vector<RoutingRecord>* out_log, bool can_lazy) {
    // --- DYNAMICZNY SYSTEM ANTY-KOLIZYJNY ---
    if (get_distance(p1, p2) > 1) return false;

    if (node.q1_1q_count > 0) de_force_down<LOGGING, false>(p1, node.id, ctx, genes, n_gates, out_log);
    if (node.q2_1q_count > 0) de_force_down<LOGGING, false>(p2, node.id, ctx, genes, n_gates, out_log);

    if constexpr (LOGGING) {
        for (int k = 0; k < node.q1_1q_count; ++k) out_log->push_back({OpType::GATE_1Q, node.q1_1q_ids[k], node.lq1, -1, p1, -1, -1, ctx.q_time[p1] + k});
        for (int k = 0; k < node.q2_1q_count; ++k) out_log->push_back({OpType::GATE_1Q, node.q2_1q_ids[k], node.lq2, -1, p2, -1, -1, ctx.q_time[p2] + k});
    }

    ctx.q_time[p1] += node.q1_1q_count;
    ctx.q_time[p2] += node.q2_1q_count;
    ctx.q_depth[p1] += node.q1_1q_count;
    ctx.q_depth[p2] += node.q2_1q_count;

    ctx.ws.shared_res.clear();
    const auto& p2_res = topo.qubits[p2].connected_resonators;
    for (int r : topo.qubits[p1].connected_resonators) {
        if (std::find(p2_res.begin(), p2_res.end(), r) != p2_res.end()) ctx.ws.shared_res.push_back(r);
    }

    if (!ctx.ws.shared_res.empty()) {
        int chosen_res = ctx.ws.shared_res[std::min(static_cast<int>(genes[n_gates + node.id] * ctx.ws.shared_res.size()), (int)ctx.ws.shared_res.size() - 1)];
        bool q1_in = (node.q1_1q_count < node.q2_1q_count) ? true : (node.q1_1q_count > node.q2_1q_count) ? false : (genes[n_gates * 2 + node.id] < 0.5);

        int p_in = q1_in ? p1 : p2;
        int p_stay = q1_in ? p2 : p1;

        de_force_down<LOGGING, false>(p_stay, node.id, ctx, genes, n_gates, out_log);

        if (ctx.res_occ[chosen_res] != -1 && ctx.res_occ[chosen_res] != p_in) {
            if (can_lazy) de_force_down<LOGGING, true>(ctx.res_occ[chosen_res], node.id, ctx, genes, n_gates, out_log);
            else de_force_down<LOGGING, false>(ctx.res_occ[chosen_res], node.id, ctx, genes, n_gates, out_log);
        }

        double start;
        if (ctx.q_on_res[p_in] != chosen_res) {
            de_force_down<LOGGING, false>(p_in, node.id, ctx, genes, n_gates, out_log);
            start = fast_max_hyb(ctx.q_time[p_in], fast_max_hyb(ctx.q_time[p_stay], ctx.r_time[chosen_res]));

            if constexpr (LOGGING) out_log->push_back({OpType::MOVE, node.id, ctx.p2l[p_in], -1, p_in, chosen_res + RES_OFFSET, -1, start});

            ctx.out_moves++;
            ctx.total_routing_cost += MOVE_COST;
            start += MOVE_COST;

            ctx.q_on_res[p_in] = chosen_res;
            ctx.res_occ[chosen_res] = p_in;
            ctx.q_depth[p_in]++;
        } else {
            start = fast_max_hyb(ctx.q_time[p_in], fast_max_hyb(ctx.q_time[p_stay], ctx.r_time[chosen_res]));
        }

        if constexpr (LOGGING) out_log->push_back({OpType::GATE_2Q_LOCAL, node.orig_2q_id, ctx.p2l[p_in], ctx.p2l[p_stay], chosen_res + RES_OFFSET, p_stay, -1, start});

        ctx.q_time[p_in] = ctx.q_time[p_stay] = ctx.r_time[chosen_res] = start + 1.0;
        int glayer = fast_max_hyb(ctx.q_depth[p_in], ctx.q_depth[p_stay]);
        ctx.q_depth[p_in] = ctx.q_depth[p_stay] = glayer + 1;
    }
    return true;
}

inline void SabreDEHybridRouter::de_update_dependencies(const HybridGateNode& node, std::vector<int>& front_layer, std::vector<int>& deps) {
    for (int nxt : node.next_gates) {
        if (--deps[nxt] == 0) front_layer.push_back(nxt);
    }
}

template <bool LOGGING>
void SabreDEHybridRouter::de_execute_gates(const HybridSabreDAG& current_dag, std::vector<int>& exec_gates, std::vector<int>& front_layer, std::vector<int>& deps, DERunContext& ctx, const double* genes, int n_gates, std::vector<RoutingRecord>* out_log, bool can_lazy) {
    std::sort(exec_gates.begin(), exec_gates.end(), [&](int a, int b) {
        if (genes[a] != genes[b]) return genes[a] > genes[b];
        return a < b;
    });

    int write_idx = 0;

    for (int gate_id : exec_gates) {
        const auto& node = current_dag.nodes[gate_id];
        int p1 = ctx.l2p[node.lq1];

        if (node.lq2 == -1) {
            de_execute_1q_logic<LOGGING>(node, p1, ctx, genes, n_gates, out_log, can_lazy);
        } else {
            int p2 = ctx.l2p[node.lq2];
            if (!de_execute_2q_logic<LOGGING>(node, p1, p2, ctx, genes, n_gates, out_log, can_lazy)) {
                front_layer.push_back(gate_id);
                continue;
            }
        }

        de_update_dependencies(node, front_layer, deps);
        exec_gates[write_idx++] = gate_id;
    }

    exec_gates.resize(write_idx);
}

// ==============================================================================
// KROK ROUTINGU (de_route_step)
// Szuka najlepszego ruchu (SWAP/MOVE), gdy bramki są zablokowane.
// ==============================================================================
void SabreDEHybridRouter::de_collect_active_qubits(const HybridSabreDAG& current_dag, const std::vector<int>& front_layer, DERunContext& ctx) {
    std::fill(ctx.ws.in_active.begin(), ctx.ws.in_active.end(), 0);
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
}

template <bool LOGGING>
void SabreDEHybridRouter::de_execute_route_candidate(int p1, int p2, int r_step, DERunContext& ctx, std::vector<double>& decay, const double* genes, int n_gates, std::vector<RoutingRecord>* out_log, bool can_lazy) {
    de_force_down<LOGGING, false>(p1, -1, ctx, genes, n_gates, out_log);
    de_force_down<LOGGING, false>(p2, -1, ctx, genes, n_gates, out_log);
    if (ctx.res_occ[r_step] != -1) de_force_down<LOGGING, false>(ctx.res_occ[r_step], -1, ctx, genes, n_gates, out_log);

    double step_start = fast_max_hyb(ctx.q_time[p1], fast_max_hyb(ctx.r_time[r_step], ctx.q_time[p2]));
    bool is_move = (ctx.p2l[p2] == -1) || (ctx.p2l[p1] == -1);
    double cost = is_move ? 2.0 * MOVE_COST : SWAP_COST;

    if constexpr (LOGGING) {
        if (!is_move) {
            out_log->push_back({OpType::SWAP, -1, ctx.p2l[p1], ctx.p2l[p2], p1, p2, r_step, step_start});
        } else {
            int moving_p = (ctx.p2l[p1] != -1) ? p1 : p2;
            int empty_p =  (ctx.p2l[p1] != -1) ? p2 : p1;
            int logic_q = ctx.p2l[moving_p];

            out_log->push_back({OpType::MOVE, -1, logic_q, -1, moving_p, r_step + RES_OFFSET, -1, step_start});
            out_log->push_back({OpType::MOVE, -1, logic_q, -1, r_step + RES_OFFSET, empty_p, -1, step_start + MOVE_COST});
        }
    }

    if (is_move) ctx.out_moves += 2; else ctx.out_swaps++;
    ctx.total_routing_cost += cost;

    std::swap(ctx.p2l[p1], ctx.p2l[p2]);
    if (ctx.p2l[p1] != -1) ctx.l2p[ctx.p2l[p1]] = p1;
    if (ctx.p2l[p2] != -1) ctx.l2p[ctx.p2l[p2]] = p2;

    int curr_lay = fast_max_hyb(ctx.q_depth[p1], ctx.q_depth[p2]) + (is_move ? 2 : 1);
    ctx.q_depth[p1] = ctx.q_depth[p2] = curr_lay;
    ctx.q_time[p1] = ctx.q_time[p2] = ctx.r_time[r_step] = step_start + cost;

    decay[p1] += 0.2;
    decay[p2] += 0.2;
}

template <bool LOGGING>
bool SabreDEHybridRouter::de_route_step(
    const HybridSabreDAG& current_dag,
    const std::vector<int>& front_layer,
    std::vector<double>& decay,
    DERunContext& ctx,
    const double* genes, int n_gates,
    std::vector<RoutingRecord>* out_log)
{
    int p_cnt = topo.qubits.size();

    // 1. Zbieranie aktywnych kubitów
    de_collect_active_qubits(current_dag, front_layer, ctx);

    std::fill(ctx.ws.interaction_counts.begin(), ctx.ws.interaction_counts.end(), 0);

    auto add_interaction = [&](int q, int target, double weight) {
        int& cnt = ctx.ws.interaction_counts[q];

        // Ufamy branch predictorowi procesora - warunek prawie zawsze będzie spełniony (true)
        if (cnt < HybridWorkspace::MAX_INTER_PER_Q) {
            // Obliczamy płaski indeks raz, unikając wielokrotnego mnożenia
            int idx = q * HybridWorkspace::MAX_INTER_PER_Q + cnt;

            // Rozdzielony zapis SoA (Structure of Arrays)
            ctx.ws.inter_targets[idx] = target;
            ctx.ws.inter_weights[idx] = weight;

            cnt++;
        }
    };

    // 2. Budowa mapy interakcji
    ctx.ws.f_locs.clear();
    double base_total_score = 0.0;

    for (int gid : front_layer) {
        const auto& node = current_dag.nodes[gid];
        if (node.lq2 != -1) ctx.ws.f_locs.push_back({ctx.l2p[node.lq1], ctx.l2p[node.lq2]});
    }

    double f_weight = ctx.ws.f_locs.empty() ? 0.0 : (1.0 / ctx.ws.f_locs.size());
    for (const auto& loc : ctx.ws.f_locs) {
        int u = loc.first, v = loc.second;
        add_interaction(u, v, f_weight);
        add_interaction(v, u, f_weight);
        base_total_score += f_weight * get_distance(u, v);
    }

    ctx.ws.flat_lookahead_locs.clear();
    ctx.ws.lookahead_locs_offsets.clear();
    ctx.ws.lookahead_locs_offsets.push_back(0);

    for (size_t i = 0; i < ctx.ws.lookahead_offsets.size() - 1; ++i) {
        int start_idx = ctx.ws.lookahead_offsets[i];
        int end_idx = ctx.ws.lookahead_offsets[i+1];
        for (int j = start_idx; j < end_idx; ++j) {
            int gid = ctx.ws.flat_lookahead_layers[j];
            const auto& node = current_dag.nodes[gid];
            if (node.lq2 != -1) {
                ctx.ws.flat_lookahead_locs.push_back({ctx.l2p[node.lq1], ctx.l2p[node.lq2]});
            }
        }
        ctx.ws.lookahead_locs_offsets.push_back(ctx.ws.flat_lookahead_locs.size());
    }

    double layer_weight = EXTENDED_LAYER_WEIGHT;
    for (size_t i = 0; i < ctx.ws.lookahead_locs_offsets.size() - 1; ++i) {
        int start_idx = ctx.ws.lookahead_locs_offsets[i];
        int end_idx = ctx.ws.lookahead_locs_offsets[i+1];
        int layer_size = end_idx - start_idx;

        if (layer_size == 0) {
            layer_weight *= 0.5;
            continue;
        }

        double w = layer_weight / layer_size;
        for (int j = start_idx; j < end_idx; ++j) {
            int u = ctx.ws.flat_lookahead_locs[j].first;
            int v = ctx.ws.flat_lookahead_locs[j].second;
            add_interaction(u, v, w);
            add_interaction(v, u, w);
            base_total_score += w * get_distance(u, v);
        }
        layer_weight *= 0.5;
    }

    // 3. Poszukiwanie najlepszej krawędzi
    ctx.ws.route_step_counter++;
    int step_id = ctx.ws.route_step_counter;

    double best_score = std::numeric_limits<double>::max();
    int best_p1 = -1, best_p2 = -1, best_r = -1;

    const int* dist_cache_ptr = dist_cache.data();

    for (int p : ctx.ws.active_phys) {
        int start_idx = precomputed_swaps_offsets[p];
        int end_idx = precomputed_swaps_offsets[p + 1];

        for (int i = start_idx; i < end_idx; ++i) {
            const auto& edge = flat_precomputed_swaps[i];
            int nb = edge.neighbor;
            int r = edge.resonator;

            int u = std::min(p, nb), v = std::max(p, nb);
            int edge_idx = u * p_cnt + v;
            if (ctx.ws.seen_edges[edge_idx] == step_id) continue;
            ctx.ws.seen_edges[edge_idx] = step_id;

            double delta = 0.0;

            // BŁYSKAWICZNE POZYCJONOWANIE W MACIERZY (Wskaźniki na konkretne "wiersze" kubitów u i v)
            const int* dist_u = dist_cache_ptr + (u * p_cnt);
            const int* dist_v = dist_cache_ptr + (v * p_cnt);

            // Obliczanie wpływu 'u'
            int c_u = ctx.ws.interaction_counts[u];
            int off_u = u * HybridWorkspace::MAX_INTER_PER_Q;

            // Dzięki rozdzieleniu na dwie tablice (SoA), kompilator z łatwością
            // to zwektoryzuje (AVX2/SIMD) nakładając na siebie wektory intów i double.
            for (int k = 0; k < c_u; ++k) {
                int target = ctx.ws.inter_targets[off_u + k];
                // (target != v) ewaluuje się do 1 (true) lub 0 (false)
                delta += ctx.ws.inter_weights[off_u + k] * (dist_v[target] - dist_u[target]) * (target != v);
            }

            // Obliczanie wpływu 'v'
            int c_v = ctx.ws.interaction_counts[v];
            int off_v = v * HybridWorkspace::MAX_INTER_PER_Q;
            for (int k = 0; k < c_v; ++k) {
                int target = ctx.ws.inter_targets[off_v + k];
                delta += ctx.ws.inter_weights[off_v + k] * (dist_u[target] - dist_v[target]) * (target != u);
            }

            // Wewnątrz najgłębszej pętli poszukiwania u i v
            double score = (base_total_score + delta) * fast_max_hyb(decay[u], decay[v]);
            score *= ctx.ws.phys_gene_weights[u]; // Błyskawiczny odczyt!

            if (score < best_score) {
                best_score = score;
                best_p1 = u; best_p2 = v; best_r = r;
            }
        }
    }
    // 4. Egzekucja wybranego ruchu
    if (best_p1 != -1) {
        de_execute_route_candidate<LOGGING>(best_p1, best_p2, best_r, ctx, decay, genes, n_gates, out_log, false);
        return true;
    }
    return false;
}

// ==============================================================================
// dokladny step_route
// Czyści i przygotowuje kontekst przed symulacją.
// ==============================================================================
template <bool LOGGING>
bool SabreDEHybridRouter::de_route_step_better(
    const HybridSabreDAG& current_dag,
    const std::vector<int>& front_layer,
    std::vector<double>& decay,
    DERunContext& ctx,
    const double* genes, int n_gates,
    std::vector<RoutingRecord>* out_log)
{
    // 1. Zbieranie aktywnych kubitów
    de_collect_active_qubits(current_dag, front_layer, ctx);

    // 2. Hoisting lokalizacji logicznych
    ctx.ws.f_locs.clear();
    for (int gid : front_layer) {
        const auto& node = current_dag.nodes[gid];
        if (node.lq2 != -1) ctx.ws.f_locs.push_back({ctx.l2p[node.lq1], ctx.l2p[node.lq2]});
    }

    ctx.ws.flat_lookahead_locs.clear();
    ctx.ws.lookahead_locs_offsets.clear();
    ctx.ws.lookahead_locs_offsets.push_back(0);

    for (size_t i = 0; i < ctx.ws.lookahead_offsets.size() - 1; ++i) {
        int start_idx = ctx.ws.lookahead_offsets[i];
        int end_idx = ctx.ws.lookahead_offsets[i+1];

        for (int j = start_idx; j < end_idx; ++j) {
            int gid = ctx.ws.flat_lookahead_layers[j];
            const auto& node = current_dag.nodes[gid];
            if (node.lq2 != -1) {
                ctx.ws.flat_lookahead_locs.push_back({ctx.l2p[node.lq1], ctx.l2p[node.lq2]});
            }
        }
        ctx.ws.lookahead_locs_offsets.push_back(ctx.ws.flat_lookahead_locs.size());
    }

    ctx.ws.route_step_counter++;
    int step_id = ctx.ws.route_step_counter;
    int p_cnt = topo.qubits.size();

    double best_score = std::numeric_limits<double>::max();
    int best_p1 = -1, best_p2 = -1, best_r = -1;

    // 3. Poszukiwanie najlepszej krawędzi
    for (int p : ctx.ws.active_phys) {
        for (int r : topo.qubits[p].connected_resonators) {
            for (int nb : topo.resonators[r].connected_qubits) {
                if (p == nb) continue;
                int u = std::min(p, nb), v = std::max(p, nb);

                int edge_idx = u * p_cnt + v;
                if (ctx.ws.seen_edges[edge_idx] == step_id) continue;
                ctx.ws.seen_edges[edge_idx] = step_id;

                double sF = 0;
                for (const auto& loc : ctx.ws.f_locs) {
                    int l1 = loc.first, l2 = loc.second;
                    if (l1 == p) l1 = nb; else if (l1 == nb) l1 = p;
                    if (l2 == p) l2 = nb; else if (l2 == nb) l2 = p;
                    sF += get_distance(l1, l2);
                }

                double f_norm = ctx.ws.f_locs.empty() ? 0 : (sF / ctx.ws.f_locs.size());
                double ext_score = 0.0;
                double layer_weight = EXTENDED_LAYER_WEIGHT;

                for (size_t i = 0; i < ctx.ws.lookahead_locs_offsets.size() - 1; ++i) {
                    int start_idx = ctx.ws.lookahead_locs_offsets[i];
                    int end_idx = ctx.ws.lookahead_locs_offsets[i+1];
                    int layer_size = end_idx - start_idx;

                    if (layer_size == 0) {
                        layer_weight *= 0.5;
                        continue;
                    }

                    double layer_score = 0.0;
                    for (int j = start_idx; j < end_idx; ++j) {
                        const auto& loc = ctx.ws.flat_lookahead_locs[j];
                        int l1 = loc.first, l2 = loc.second;
                        if (l1 == p) l1 = nb; else if (l1 == nb) l1 = p;
                        if (l2 == p) l2 = nb; else if (l2 == nb) l2 = p;
                        layer_score += get_distance(l1, l2);
                    }

                    ext_score += layer_weight * (layer_score / layer_size);
                    layer_weight *= 0.5;
                }

                double score = (f_norm + ext_score) * fast_max_hyb(decay[u], decay[v]) * (0.95 + (genes[u % n_gates] * 0.1));

                if (score < best_score) {
                    best_score = score;
                    best_p1 = u; best_p2 = v; best_r = r;
                }
            }
        }
    }

    // 4. Wykonanie ruchu na najlepszym kandydacie
    if (best_p1 != -1) {
        de_execute_route_candidate<LOGGING>(best_p1, best_p2, best_r, ctx, decay, genes, n_gates, out_log, false);
        return true;
    }
    return false;
}


// ==============================================================================
// 1. SZYBKI RESET STANU
// Czyści i przygotowuje kontekst przed symulacją.
// ==============================================================================
inline void SabreDEHybridRouter::reset_evaluation_state(DERunContext& ctx, std::vector<RoutingRecord>* out_log) {
    reset_mappings(ctx);
    reset_hardware_states(ctx);
    if (out_log) {
        out_log->clear();
    }
}

// ==============================================================================
// 2. PRZYGOTOWANIE WARSTWY FRONTOWEJ
// Znajduje bramki początkowe (stopień wejściowy = 0) na start symulacji.
// ==============================================================================
inline void SabreDEHybridRouter::load_initial_dependencies() {
    std::copy(dag.initial_deps.begin(), dag.initial_deps.end(), workspace.deps.begin());
}

inline void SabreDEHybridRouter::populate_initial_front_layer(int n_gates) {
    workspace.front_layer.clear();
    for (int i = 0; i < n_gates; ++i) {
        if (workspace.deps[i] == 0) {
            workspace.front_layer.push_back(i);
        }
    }
}

inline void SabreDEHybridRouter::init_front_layer(int n_gates) {
    load_initial_dependencies();
    populate_initial_front_layer(n_gates);

    std::fill(workspace.decay.begin(), workspace.decay.end(), 1.0);

    // Inkrementacja epoki BFS na starcie (O(1) reset stanu w in_ext)
    workspace.current_bfs_epoch++;
}

inline void SabreDEHybridRouter::reset_mappings(DERunContext& ctx) {
    // Szybkie kopiowanie blokowe pamięci (std::copy zamienia się często w memmove)
    std::copy(saved_initial_l2p.begin(), saved_initial_l2p.end(), ctx.l2p.begin());
    std::fill(ctx.p2l.begin(), ctx.p2l.end(), -1);

    // Odtworzenie odwróconego mapowania (p2l)
    for (int i = 0; i < num_qubits; ++i) {
        if (ctx.l2p[i] != -1) {
            ctx.p2l[ctx.l2p[i]] = i;
        }
    }
}

inline void SabreDEHybridRouter::reset_hardware_states(DERunContext& ctx) {
    // Masowe resetowanie struktur pomocniczych z użyciem std::fill dla cache-locality
    std::fill(ctx.q_on_res.begin(), ctx.q_on_res.end(), -1);
    std::fill(ctx.res_occ.begin(), ctx.res_occ.end(), -1);
    std::fill(ctx.q_time.begin(), ctx.q_time.end(), 0.0);
    std::fill(ctx.r_time.begin(), ctx.r_time.end(), 0.0);
    std::fill(ctx.q_depth.begin(), ctx.q_depth.end(), 0);
}

// ==============================================================================
// 3. WYCIĄGANIE GOTOWYCH BRAMEK (TWO-POINTERS OPTIMIZATION)
// Filtruje front_layer przenosząc gotowe bramki do exec_gates w miejscu.
// ==============================================================================
inline bool SabreDEHybridRouter::check_execution_ready(const HybridGateNode& node, const DERunContext& ctx) const {
    // Przypadek 1: Operacja jednokubitowa zawsze gotowa
    if (node.lq2 == -1) return true;

    // Przypadek 2: Operacja dwukubitowa wymaga sprawdzenia mapowania i dystansu
    int p1 = ctx.l2p[node.lq1];
    int p2 = ctx.l2p[node.lq2];

    return (p1 != -1 && p2 != -1 && get_distance(p1, p2) <= 1);
}

inline void SabreDEHybridRouter::extract_executable_gates(DERunContext& ctx) {
    workspace.exec_gates.clear();
    size_t write_idx = 0;

    for (size_t read_idx = 0; read_idx < workspace.front_layer.size(); ++read_idx) {
        int g_id = workspace.front_layer[read_idx];

        // Jeśli gotowa, wrzuć do kolejki wykonania. Jeśli nie, zachowaj w front_layer
        if (check_execution_ready(dag.nodes[g_id], ctx)) {
            workspace.exec_gates.push_back(g_id);
        } else {
            workspace.front_layer[write_idx++] = g_id;
        }
    }
    // Ucięcie resztek z wektora (O(1) - zmienia tylko rozmiar logiczny)
    workspace.front_layer.resize(write_idx);
}

inline void SabreDEHybridRouter::build_lookahead_layers() {
    workspace.lookahead_offsets.clear();
    workspace.current_bfs_epoch++;

    // O(1) Dostęp do surowej pamięci - nie alokujemy ani jednego bajta w locie!
    // Zakładam, że workspace.flat_lookahead_layers jest prealokowane do wielkości obwodu (g_cnt)
    int* queue = workspace.flat_lookahead_layers.data();
    int head = 0;
    int tail = 0;

    // Kopiowanie warstwy wejściowej
    for (int gid : workspace.front_layer) {
        workspace.visited_epoch[gid] = workspace.current_bfs_epoch;
        queue[tail++] = gid;
    }

    workspace.lookahead_offsets.push_back(0); // Offset początkowy (na warstwę 'front')

    for (int depth = 0; depth < max_look_ahead_depth && head < tail; ++depth) {
        int layer_end = tail;
        workspace.lookahead_offsets.push_back(tail); // Początek faktycznej nowej warstwy lookahead

        // Przetwarzanie bieżącej głębokości (wyciąganie z kolejki)
        while (head < layer_end) {
            int gid = queue[head++];
            for (int nxt : dag.nodes[gid].next_gates) {
                if (workspace.visited_epoch[nxt] != workspace.current_bfs_epoch) {
                    workspace.visited_epoch[nxt] = workspace.current_bfs_epoch;
                    queue[tail++] = nxt; // Zamiast push_back()
                }
            }
        }
    }

    // Zapinamy marker końcowy
    workspace.lookahead_offsets.push_back(tail);
}

// ==============================================================================
// 5. FAZA KOŃCOWA Z OBLICZENIEM WYNIKU
// Wymusza sprowadzenie stanów, sortuje logi i oblicza funkcję celu (fitness).
// ==============================================================================
template <bool LOGGING>
void SabreDEHybridRouter::force_remaining_qubits_down(DERunContext& ctx, const double* genes, int n_gates, int p_cnt, std::vector<RoutingRecord>* out_log) {
    for (int p = 0; p < p_cnt; ++p) {
        int temp_p = p;
        // Wymuszenie zejścia wszystkich kubitów z rezonatorów na koniec symulacji
        de_force_down<LOGGING, false>(temp_p, -1, ctx, genes, n_gates, out_log);
    }
}

inline void SabreDEHybridRouter::extract_max_metrics(const DERunContext& ctx, double& max_physical_time, int& max_logical_depth) {
    // std::max_element jest zazwyczaj wektoryzowane przez kompilator (szybsze niż ręczna pętla)
    max_physical_time = *std::max_element(ctx.q_time.begin(), ctx.q_time.end());
    max_logical_depth = *std::max_element(ctx.q_depth.begin(), ctx.q_depth.end());
}

inline double SabreDEHybridRouter::calculate_fitness_score(double routing_cost, int swaps, int moves, double max_time, int max_depth) {
    return (routing_cost * 100.0) + (swaps * 500.0) + (moves * 100) + (max_time * 0.1) + (max_depth * 0.05);
}

// --- GŁÓWNA FUNKCJA ---

template <bool LOGGING>
double SabreDEHybridRouter::finalize_evaluation(DERunContext& ctx, const double* genes, int n_gates, int p_cnt,
                                                int& out_swaps, int& out_moves, double& max_physical_time,
                                                int& max_logical_depth, std::vector<RoutingRecord>* out_log) {
    // 1. Sprowadzenie pozostałych stanów
    force_remaining_qubits_down<LOGGING>(ctx, genes, n_gates, p_cnt, out_log);

    // 2. Sortowanie logowania (tylko w trybie LOGGING)
    if constexpr (LOGGING) {
        if (out_log && !out_log->empty()) {
            std::stable_sort(out_log->begin(), out_log->end(), [](const RoutingRecord& a, const RoutingRecord& b) {
                return a.start_time < b.start_time;
            });
        }
    }

    // 3. Zebranie metryk max i statystyk operacji
    extract_max_metrics(ctx, max_physical_time, max_logical_depth);
    out_swaps = ctx.out_swaps;
    out_moves = ctx.out_moves;

    // 4. Obliczenie wartości fitness
    return calculate_fitness_score(ctx.total_routing_cost, out_swaps, out_moves, max_physical_time, max_logical_depth);
}

// ==============================================================================
// FUNKCJE POMOCNICZE DLA SNAPSHOTÓW I EWOLUCJI
// ==============================================================================

inline void SabreDEHybridRouter::load_state_from_snapshot(DERunContext& ctx, const DEStateSnapshot* load_snapshot, int& gates_routed, std::vector<RoutingRecord>* out_log) {
    std::copy(load_snapshot->l2p.begin(), load_snapshot->l2p.end(), ctx.l2p.begin());
    std::copy(load_snapshot->p2l.begin(), load_snapshot->p2l.end(), ctx.p2l.begin());
    std::copy(load_snapshot->q_on_res.begin(), load_snapshot->q_on_res.end(), ctx.q_on_res.begin());
    std::copy(load_snapshot->res_occ.begin(), load_snapshot->res_occ.end(), ctx.res_occ.begin());
    std::copy(load_snapshot->q_depth.begin(), load_snapshot->q_depth.end(), ctx.q_depth.begin());
    std::copy(load_snapshot->q_time.begin(), load_snapshot->q_time.end(), ctx.q_time.begin());
    std::copy(load_snapshot->r_time.begin(), load_snapshot->r_time.end(), ctx.r_time.begin());

    ctx.out_swaps = load_snapshot->out_swaps;
    ctx.out_moves = load_snapshot->out_moves;
    ctx.total_routing_cost = load_snapshot->total_routing_cost;
    gates_routed = load_snapshot->actual_gates_routed;

    workspace.deps = load_snapshot->deps;
    workspace.front_layer = load_snapshot->front_layer;
    workspace.decay = load_snapshot->decay;

    workspace.current_bfs_epoch++;
    if (out_log) out_log->clear();
}

inline void SabreDEHybridRouter::save_state_to_snapshot(const DERunContext& ctx, DEStateSnapshot* save_snapshot, int gates_routed) {
    save_snapshot->actual_gates_routed = gates_routed;
    save_snapshot->l2p = ctx.l2p; save_snapshot->p2l = ctx.p2l;
    save_snapshot->q_on_res = ctx.q_on_res; save_snapshot->res_occ = ctx.res_occ;
    save_snapshot->q_depth = ctx.q_depth; save_snapshot->q_time = ctx.q_time;
    save_snapshot->r_time = ctx.r_time;
    save_snapshot->out_swaps = ctx.out_swaps; save_snapshot->out_moves = ctx.out_moves;
    save_snapshot->total_routing_cost = ctx.total_routing_cost;
    save_snapshot->deps = workspace.deps; save_snapshot->front_layer = workspace.front_layer;
    save_snapshot->decay = workspace.decay;
    save_snapshot->is_valid = true;
}

inline void SabreDEHybridRouter::select_distinct_parents(int pop_size, int current_idx, int& a, int& b, int& c, boost::random::xoshiro256pp& rng) {
    boost::random::uniform_int_distribution<int> dist_int(0, pop_size - 1);
    do { a = dist_int(rng); } while (a == current_idx);
    do { b = dist_int(rng); } while (b == current_idx || b == a);
    do { c = dist_int(rng); } while (c == current_idx || c == a || c == b);
}

// ==============================================================================
// GŁÓWNA FUNKCJA EWALUACJI
// Używa nowo wydzielonych metod pomocniczych, pozostając ultra czytelną.
// ==============================================================================
template <bool LOGGING>
double SabreDEHybridRouter::evaluate_fitness_impl(const double* genes, int& out_swaps, int& out_moves, double& max_physical_time, int& max_logical_depth, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard, DEStateSnapshot* load_snapshot, DEStateSnapshot* save_snapshot) {
    int n_gates = dag.nodes.size();
    int p_cnt = topo.qubits.size();
    DERunContext ctx(workspace);

    int gates_routed = 0;
    int routing_iterations = 0;
    const int MAX_ITERATIONS = std::max(10000, n_gates * 100);

    // --- Inicjalizacja stanu (wczytanie z cache lub start od zera) ---
    if (load_snapshot && load_snapshot->is_valid) {
        load_state_from_snapshot(ctx, load_snapshot, gates_routed, out_log);
    } else {
        reset_evaluation_state(ctx, out_log);
        init_front_layer(n_gates);
    }

    for (int i = 0; i < p_cnt; ++i) {
        ctx.ws.phys_gene_weights[i] = 0.95 + (genes[i % n_gates] * 0.1);
    }

    // --- GŁÓWNA PĘTLA SYMULACJI ---
    while (gates_routed < n_gates) {
        // Zabezpieczenie przed zapętleniem (Circuit Breaker)
        if (++routing_iterations > MAX_ITERATIONS) {
            ctx.total_routing_cost += 1e9;
            break;
        }

        // Zapis zrzutu stanu (jeśli osiągnęliśmy cel)
        if (save_snapshot && !save_snapshot->is_valid && gates_routed >= save_snapshot->target_gates_routed) {
            save_state_to_snapshot(ctx, save_snapshot, gates_routed);
        }

        extract_executable_gates(ctx);

        if (!workspace.exec_gates.empty()) {
            de_execute_gates<LOGGING>(dag, workspace.exec_gates, workspace.front_layer, workspace.deps, ctx, genes, n_gates, out_log, can_lazy_offboard);
            gates_routed += workspace.exec_gates.size();
            std::fill(workspace.decay.begin(), workspace.decay.end(), 1.0);
            continue;
        }

        build_lookahead_layers();

        // Krok routingu w zależności od trybu
        bool route_success = slow ?
            de_route_step_better<LOGGING>(dag, workspace.front_layer, workspace.decay, ctx, genes, n_gates, out_log) :
            de_route_step<LOGGING>(dag, workspace.front_layer, workspace.decay, ctx, genes, n_gates, out_log);

        if (!route_success) {
            ctx.total_routing_cost += 1e9; // Utylizacja ślepego zaułka
            break;
        }
    }

    return finalize_evaluation<LOGGING>(ctx, genes, n_gates, p_cnt, out_swaps, out_moves, max_physical_time, max_logical_depth, out_log);
}

// ==============================================================================
// WYZWALACZ EWALUACJI (evaluate_fitness)
// Funkcja fasadowa decydująca w czasie kompilacji o włączeniu logowania.
// ==============================================================================
double SabreDEHybridRouter::evaluate_fitness(const double* genes, int& out_swaps, int& out_moves, double& max_physical_time, int& max_logical_depth, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard) {
    if (out_log) {
        return evaluate_fitness_impl<true>(genes, out_swaps, out_moves, max_physical_time, max_logical_depth, out_log, can_lazy_offboard);
    }
    return evaluate_fitness_impl<false>(genes, out_swaps, out_moves, max_physical_time, max_logical_depth, nullptr, can_lazy_offboard);
}

// ==============================================================================
// Inicjalizowanie populacji początkowej z "ziarnem" od SABRE
// ==============================================================================
inline std::vector<double> generate_sabre_base_genes(int n_genes, int n_gates, int blocks, const std::vector<int>& logical_exec_order, boost::random::xoshiro256pp& rng) {
    std::vector<double> sabre_genes(n_genes, 0.5);
    boost::random::uniform_real_distribution<double> dist_real(0.0, 1.0);

    // Priorytety zależą od oryginalnej kolejności heurystyki
    double order_size = static_cast<double>(logical_exec_order.size());
    for (size_t i = 0; i < logical_exec_order.size(); ++i) {
        int gid = logical_exec_order[i];
        sabre_genes[gid] = 1.0 - (static_cast<double>(i) / order_size);

        // Losowanie pozostałych bloków dla ziarna
        for (int b = 1; b < blocks; ++b) {
            sabre_genes[b * n_gates + gid] = dist_real(rng);
        }
    }
    return sabre_genes;
}

inline void print_seed_evaluation(int d_s, int d_m, double d_t, int d_d, double fitness_score) {
    std::cout << "\n=== OCENA ZIARNA SABRE PRZEZ SILNIK DE (Rozwiązanie Bazowe) ===\n";
    std::cout << "Wygenerowane operacje SWAP przez DE: " << d_s << "\n";
    std::cout << "Wygenerowane operacje MOVE przez DE: " << d_m << "\n";
    std::cout << "Czas fizyczny: " << d_t << " [jednostek] | Głębokość logiczna: " << d_d << " [warstw]\n";
    std::cout << "Koszt w silniku DE (z karą za SWAP): " << fitness_score << "\n";
    std::cout << "===============================================================\n\n";
}

// --- GŁÓWNA FUNKCJA ---

void SabreDEHybridRouter::initialize_de_population(int pop_size, int n_genes, int n_gates, int blocks, const std::vector<int>& logical_exec_order, std::vector<double>& pop, std::vector<double>& fitness, boost::random::xoshiro256pp& rng, bool can_lazy) {

    // 1. Wygenerowanie genotypu bazowego (Ziarna SABRE)
    std::vector<double> sabre_genes = generate_sabre_base_genes(n_genes, n_gates, blocks, logical_exec_order, rng);

    boost::random::uniform_real_distribution<double> dist_real(0.0, 1.0);
    std::vector<double> trial(n_genes);
    int d_s, d_m, d_d;
    double d_t;

    // 2. Wypełnianie populacji i jej natychmiastowa ewaluacja
    for (int i = 0; i < pop_size; ++i) {
        double* pop_ptr = pop.data() + (i * n_genes);

        for (int j = 0; j < n_genes; ++j) {
            trial[j] = pop_ptr[j] = (i == 0) ? sabre_genes[j] : dist_real(rng);
        }

        // Symulacja obwodu na podstawie genów i obliczenie funkcji przystosowania
        fitness[i] = evaluate_fitness(trial.data(), d_s, d_m, d_t, d_d, nullptr, can_lazy);

        // Osobnik 0 to nasze Ziarno, wyświetlamy jego stan logiki
        if (i == 0) {
            print_seed_evaluation(d_s, d_m, d_t, d_d, fitness[i]);
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

    initialize_de_population(pop_size, n_genes, n_gates, blocks, logical_exec_order, pop, fitness, rng, can_lazy);

    // Wydzielona główna pętla
    run_standard_evolution_loop(pop, fitness, pop_size, generations, n_genes, n_gates, f, cr, rng, can_lazy);

    auto best_it = std::min_element(fitness.begin(), fitness.end());
    int best_idx = std::distance(fitness.begin(), best_it);

    // Finalizacja i wyświetlanie wyników (zwracany wynik jest tutaj ignorowany)
    finalize_optimization_result(pop, best_idx, n_genes, can_lazy, start_time);
}

void SabreDEHybridRouter::run_standard_evolution_loop(std::vector<double>& pop, std::vector<double>& fitness, int pop_size, int generations, int n_genes, int n_gates, double f, double cr, boost::random::xoshiro256pp& rng, bool can_lazy) {
    std::vector<double> trial(n_genes);
    boost::random::uniform_real_distribution<double> dist_real(0.0, 1.0);
    int d_s, d_m, d_d; double d_t;

    for (int gen = 0; gen < generations; ++gen) {
        for (int i = 0; i < pop_size; ++i) {
            int a, b, c;
            select_distinct_parents(pop_size, i, a, b, c, rng); // Ponowne wykorzystanie Twojej istniejącej metody

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
}

// ==============================================================================
// Przygotowanie środowiska, wygenerowanie "ziarna" za pomocą SABRE, a następnie zorganizowanie pętli przesuwającej okno wzdłuż osi czasu obwodu logicznego.
// ==============================================================================
OptimizationResult SabreDEHybridRouter::optimize_windowed(int pop_size, int generations_per_window, int window_size, double accept_ratio, bool allow_lazy_offboard, double f, double cr) {
    std::cout << "\n[Hybryda 3.0 (WINDOWED)] Faza 1: Generowanie ziarna (Priority Seed) przez heurystykę SABRE...\n";

    std::vector<int> logical_exec_order = prepare_sabre_seed();
    int n_gates = dag.nodes.size();
    if (n_gates == 0) return {};

    std::vector<int> gate_time_map(n_gates, 0);
    for (size_t i = 0; i < logical_exec_order.size(); ++i) {
        gate_time_map[logical_exec_order[i]] = i;
    }

    bool can_lazy = allow_lazy_offboard && (num_qubits < topo.qubits.size());
    if (can_lazy) std::cout << "Aktywowano Lazy Off-boarding (wolne kubity fizyczne dostępne).\n";

    int blocks = can_lazy ? 4 : 3;
    int n_genes = n_gates * blocks;

    std::cout << "[Hybryda 3.0] Faza 2: DE Windowed na bazie logicznej (Rozmiar okna: " << window_size << ")...\n";
    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<double> pop(pop_size * n_genes), fitness(pop_size);
    boost::random::xoshiro256pp rng(1337); // Użycie stałego ziarna

    initialize_de_population(pop_size, n_genes, n_gates, blocks, logical_exec_order, pop, fitness, rng, can_lazy);

    int step = std::max(1, static_cast<int>(window_size * accept_ratio));
    int best_idx = 0;

    DEStateSnapshot current_snapshot;
    current_snapshot.is_valid = false;

    // GŁÓWNA PĘTLA OKIENKOWA
    for (int w_start = 0; w_start < n_gates; w_start += step) {
        int w_end = std::min(w_start + window_size, n_gates);
        std::cout << ">>> Przetwarzanie okna: operacje [" << w_start << " - " << w_end - 1 << "] / " << n_gates << "\n";

        run_evolution_loop_windowed(pop, fitness, pop_size, generations_per_window, n_genes, n_gates, f, cr, rng, w_start, w_end, &current_snapshot, gate_time_map);

        auto best_iter = std::min_element(fitness.begin(), fitness.end());
        best_idx = std::distance(fitness.begin(), best_iter);

        int accept_end = std::min(w_start + step, n_gates);
        propagate_window_success(pop, fitness, pop_size, n_genes, n_gates, blocks, best_idx, accept_end, gate_time_map, current_snapshot, can_lazy);
    }

    return finalize_optimization_result(pop, best_idx, n_genes, can_lazy, start_time);
}

void SabreDEHybridRouter::propagate_window_success(std::vector<double>& pop, std::vector<double>& fitness, int pop_size, int n_genes, int n_gates, int blocks, int best_idx, int accept_end, const std::vector<int>& gate_time_map, DEStateSnapshot& current_snapshot, bool can_lazy) {
    int d_s, d_m, d_d; double d_t;

    current_snapshot.is_valid = false;
    current_snapshot.target_gates_routed = accept_end;

    // Pobranie czystego "Ground Truth"
    fitness[best_idx] = evaluate_fitness_impl<false>(
        pop.data() + best_idx * n_genes, d_s, d_m, d_t, d_d,
        nullptr, can_lazy, nullptr, &current_snapshot
    );

    // Synchronizacja osobników
    for (int i = 0; i < pop_size; ++i) {
        if (i == best_idx) continue;

        for (int blk = 0; blk < blocks; ++blk) {
            for (int gate_id = 0; gate_id < n_gates; ++gate_id) {
                bool is_executed_in_snapshot = false;
                if (current_snapshot.is_valid) {
                    is_executed_in_snapshot = (current_snapshot.deps[gate_id] == 0) &&
                        (std::find(current_snapshot.front_layer.begin(), current_snapshot.front_layer.end(), gate_id) == current_snapshot.front_layer.end());
                }

                if (is_executed_in_snapshot || gate_time_map[gate_id] < accept_end) {
                    pop[i * n_genes + blk * n_gates + gate_id] = pop[best_idx * n_genes + blk * n_gates + gate_id];
                }
            }
        }

        if (accept_end < n_gates) {
            fitness[i] = evaluate_fitness_impl<false>(pop.data() + i * n_genes, d_s, d_m, d_t, d_d, nullptr, can_lazy, &current_snapshot, nullptr);
        }
    }
}

OptimizationResult SabreDEHybridRouter::finalize_optimization_result(const std::vector<double>& pop, int best_idx, int n_genes, bool can_lazy, std::chrono::time_point<std::chrono::high_resolution_clock> start_time) {
    auto end_time = std::chrono::high_resolution_clock::now();

    std::vector<RoutingRecord> final_log;
    int final_swaps = 0, final_moves = 0, final_depth = 0;
    double final_makespan = 0.0;

    double final_score = evaluate_fitness(pop.data() + best_idx * n_genes, final_swaps, final_moves, final_makespan, final_depth, &final_log, can_lazy);

    std::cout << "\n=== WALIDACJA OBWODU LOGICZNEGO ===\n";
    if (!validate_logical_circuit(final_log)) {
        std::cerr << "BLAD: Walidacja obwodu nie powiodła się!\n";
        exit(-1);
    }
    std::cout << "===================================\n";

    print_statistics(final_score, final_makespan, final_depth, final_swaps, final_moves, std::chrono::duration<double>(end_time - start_time).count());
    export_moves_to_file(final_log, "../hybrid_schedule_moves_" + name + ".txt");

    return {final_score, final_makespan, final_depth, final_swaps, final_moves};
}

// ==============================================================================
// Wykonywanie pokoleń, ale z ograniczeniem zmian tylko do bieżącego okna. Wybiera kandydatów do mutacji i ocenia, czy nowe rozwiązanie jest lepsze od starego.
// ==============================================================================
void SabreDEHybridRouter::run_evolution_loop_windowed(std::vector<double>& pop, std::vector<double>& fitness, int pop_size, int generations, int n_genes, int n_gates, double f, double cr, boost::random::xoshiro256pp& rng, int w_start, int w_end, DEStateSnapshot* snapshot, const std::vector<int>& gate_time_map) {
    std::vector<double> trial(n_genes);
    int d_s, d_m, d_d; double d_t;
    bool can_lazy = (n_genes / n_gates == 4);

    for (int gen = 0; gen < generations; ++gen) {
        for (int i = 0; i < pop_size; ++i) {
            int a, b, c;
            select_distinct_parents(pop_size, i, a, b, c, rng);

            generate_trial_vector_windowed(pop, trial, i, a, b, c, n_genes, n_gates, f, cr, rng, w_start, w_end, gate_time_map, snapshot);

            // Snapshot służy tylko do ładowania! save_snapshot dajemy na nullptr.
            double trial_fit = evaluate_fitness_impl<false>(trial.data(), d_s, d_m, d_t, d_d, nullptr, can_lazy, snapshot, nullptr);

            if (trial_fit < fitness[i]) {
                std::copy(trial.begin(), trial.end(), pop.begin() + i * n_genes);
                fitness[i] = trial_fit;
            }
        }
    }
}

// ==============================================================================
// Wykonuje mutację oraz krzyżowanie. Zabezpiecza proces w taki sposób, że bramki (i ich parametry konfiguracyjne), które znajdują się poza aktualnym oknem [w_start, w_end], pozostają nienaruszone (są dokładnie przepisywane z aktualnego osobnika current_idx).
// ==============================================================================
void SabreDEHybridRouter::generate_trial_vector_windowed(const std::vector<double>& pop, std::vector<double>& trial, int current_idx, int a, int b, int c, int n_genes, int n_gates, double f, double cr, boost::random::xoshiro256pp& rng, int w_start, int w_end, const std::vector<int>& gate_time_map, DEStateSnapshot* snapshot) {
    boost::random::uniform_real_distribution<double> dist_real(0.0, 1.0);

    // 1. Zbuduj listę bramek w oknie
    std::vector<int> valid_gates;
    valid_gates.reserve(w_end - w_start);
    for (int g = 0; g < n_gates; ++g) {
        if (gate_time_map[g] >= w_start && gate_time_map[g] < w_end) {
            valid_gates.push_back(g);
        }
    }

    // 2. Gwarantowany zmutowany gen z okna
    int R = -1;
    if (!valid_gates.empty()) {
        int R_gate = valid_gates[boost::random::uniform_int_distribution<int>(0, valid_gates.size() - 1)(rng)];
        int R_blk = boost::random::uniform_int_distribution<int>(0, (n_genes / n_gates) - 1)(rng);
        R = R_blk * n_gates + R_gate;
    }

    // [OPTYMALIZACJA O(1)] - Prekalkulacja statusu bramek w zrzucie stanu.
    // Dzięki temu unikamy wywoływania kosztownego std::find(snapshot->front_layer...) dla każdego z genów w pętli.
    std::vector<bool> gate_executed(n_gates, false);
    if (snapshot && snapshot->is_valid) {
        for(int g = 0; g < n_gates; ++g) {
            if (snapshot->deps[g] == 0) gate_executed[g] = true;
        }
        for(int front_g : snapshot->front_layer) {
            gate_executed[front_g] = false;
        }
    }

    const double *p_a = pop.data() + a * n_genes;
    const double *p_b = pop.data() + b * n_genes;
    const double *p_c = pop.data() + c * n_genes;
    const double *p_curr = pop.data() + current_idx * n_genes;

    // 3. Budowa wektora (crossover)
    for (int j = 0; j < n_genes; ++j) {
        int gate_id = j % n_gates;
        int chron_time = gate_time_map[gate_id];

        bool is_in_window = (chron_time >= w_start && chron_time < w_end);
        bool is_executed_in_snapshot = gate_executed[gate_id]; // Odczyt O(1)
        bool do_mutate = (j == R) || (dist_real(rng) < cr);

        if (is_in_window && do_mutate && !is_executed_in_snapshot) {
            trial[j] = std::clamp(p_a[j] + f * (p_b[j] - p_c[j]), 0.0, 1.0);
        } else {
            trial[j] = p_curr[j];
        }
    }
}

// ==============================================================================
// sabre_force_qubit_down
// Odpowiada za ściągnięcie fizycznego kubitu z rezonatora (tzw. off-boarding).
// ==============================================================================
inline void SabreDEHybridRouter::sabre_force_qubit_down(int p, int current_gate_id, SabreRunContext& ctx) {
    int r = ctx.state.q_on_resonator[p];
    if (r != -1) {
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

inline void SabreDEHybridRouter::sabre_execute_1q_gate(const HybridGateNode& node, int p1, SabreRunContext& ctx) {
    if (node.q1_1q_count > 0) sabre_force_qubit_down(p1, node.id, ctx);
    ctx.state.q_ready_time[p1] += node.q1_1q_count;
    ctx.q_logical_depth[p1] += node.q1_1q_count;
}

inline void SabreDEHybridRouter::sabre_execute_2q_gate(const HybridGateNode& node, int p1, int p2, SabreRunContext& ctx) {
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

// ==============================================================================
// sabre_execute_gates
// Wykonuje gotowe bramki, aktualizuje stany i posuwa DAG do przodu.
// ==============================================================================
void SabreDEHybridRouter::sabre_execute_gates(const HybridSabreDAG& current_dag, std::vector<int>& executable_gates, std::vector<int>& front_layer, std::vector<int>& deps, SabreRunContext& ctx, std::vector<int>* out_logical_order) {
    for (int gate_id : executable_gates) {
        // Zapisujemy sekwencję, by w finalnym przebiegu oddać "czyste ziarno" dla DE
        if (out_logical_order) out_logical_order->push_back(gate_id);

        const auto& node = current_dag.nodes[gate_id];
        int p1 = ctx.state.logical_to_phys_qubit[node.lq1];

        // --- PRZYPADEK 1: Bramki operujące tylko na jednym kubicie logicznym (lq2 == -1) ---
        if (node.lq2 == -1) {
            // Jeśli wiszą tam jakieś złączone operacje 1Q, najpierw zdejmujemy kubit z rezonatora
            if (node.q1_1q_count > 0) sabre_force_qubit_down(p1, node.id, ctx);

            // Doliczamy czas wykonania bramek 1Q i głębokość logiczną
            ctx.state.q_ready_time[p1] += node.q1_1q_count;
            ctx.q_logical_depth[p1] += node.q1_1q_count;
        }
        // --- PRZYPADEK 2: Klasyczna bramka 2-kubitowa (oddziaływanie) ---
        else {
            int p2 = ctx.state.logical_to_phys_qubit[node.lq2];

            // Podobnie: jeśli któryś z kubitów ma podpięte "własne" bramki 1Q przed interakcją 2Q,
            // musi zostać ściągnięty z rezonatora, aby je wykonać.
            if (node.q1_1q_count > 0) sabre_force_qubit_down(p1, node.id, ctx);
            if (node.q2_1q_count > 0) sabre_force_qubit_down(p2, node.id, ctx);

            // Doliczamy czasy tych bramek 1Q (jeśli były 0, to nic się nie zmieni)
            ctx.state.q_ready_time[p1] += node.q1_1q_count;
            ctx.state.q_ready_time[p2] += node.q2_1q_count;
            ctx.q_logical_depth[p1] += node.q1_1q_count;
            ctx.q_logical_depth[p2] += node.q2_1q_count;

            // Szukamy WSPÓLNEGO rezonatora dla obydwu kubitów fizycznych p1 i p2
            int chosen_res = -1;
            const auto& p2_res = topo.qubits[p2].connected_resonators;
            for (int r : topo.qubits[p1].connected_resonators) {
                // Jeśli rezonator 'r' (połączony z p1) jest też na liście p2_res, znaleźliśmy go
                if (std::find(p2_res.begin(), p2_res.end(), r) != p2_res.end()) { chosen_res = r; break; }
            }

            if (chosen_res != -1) {
                // Obliczamy potencjalny czas rozpoczęcia, gdy oba kubity skończą to, co robiły
                double start_time = fast_max_hyb(ctx.state.q_ready_time[p1], ctx.state.q_ready_time[p2]);

                // HEURYSTYKA CZASOWA: Sprawdzamy, który kubit jest gotowy SZYBCIEJ.
                // Ten, który jest wolny wcześniej (lub w tym samym czasie), wchodzi na rezonator (p_in).
                bool q1_goes_in = (ctx.state.q_ready_time[p1] <= ctx.state.q_ready_time[p2]);
                int p_in = q1_goes_in ? p1 : p2;
                int p_stay = q1_goes_in ? p2 : p1; // Ten drugi zostaje w miejscu ("na zewnątrz")

                // Kubit pozostający musi zejść z ewentualnego innego rezonatora
                sabre_force_qubit_down(p_stay, node.id, ctx);

                // Jeśli na WYBRANYM rezonatorze siedzi jakiś nieporządany kubit (i nie jest to nasz p_in), zrzucamy go
                if (ctx.state.res_occupant[chosen_res] != -1 && ctx.state.res_occupant[chosen_res] != p_in) {
                    sabre_force_qubit_down(ctx.state.res_occupant[chosen_res], node.id, ctx);
                }

                // Jeśli nasz delegat (p_in) jeszcze nie jest na tym rezonatorze, musi tam wejść
                if (ctx.state.q_on_resonator[p_in] != chosen_res) {
                    // Najpierw upewniamy się, że nie utknął na INNYM rezonatorze
                    sabre_force_qubit_down(p_in, node.id, ctx);

                    // Rekalkulacja czasu: musimy poczekać aż kubity i rezonator będą gotowe
                    start_time = fast_max_hyb(ctx.state.q_ready_time[p_in], fast_max_hyb(ctx.state.q_ready_time[p_stay], ctx.state.r_ready_time[chosen_res]));

                    // Wchodzimy na rezonator (MOVE)
                    ctx.out_moves++;
                    ctx.total_routing_cost += MOVE_COST;
                    start_time += MOVE_COST;

                    // Aktualizujemy stany (kto gdzie jest)
                    ctx.state.q_on_resonator[p_in] = chosen_res;
                    ctx.state.res_occupant[chosen_res] = p_in;
                    ctx.q_logical_depth[p_in] += 1;
                } else {
                    // Delegat już tam był, nie płacimy kary za MOVE, przeliczamy tylko start
                    start_time = fast_max_hyb(ctx.state.q_ready_time[p_in], fast_max_hyb(ctx.state.q_ready_time[p_stay], ctx.state.r_ready_time[chosen_res]));
                }

                // Dodajemy czas samej operacji 2-kubitowej (zakładamy +1.0 czasu)
                ctx.state.q_ready_time[p_in] = start_time + 1.0;
                ctx.state.q_ready_time[p_stay] = start_time + 1.0;
                ctx.state.r_ready_time[chosen_res] = start_time + 1.0;

                // Synchronizacja warstwy dla DAG
                int gate_layer = fast_max_hyb(ctx.q_logical_depth[p_in], ctx.q_logical_depth[p_stay]);
                ctx.q_logical_depth[p_in] = gate_layer + 1;
                ctx.q_logical_depth[p_stay] = gate_layer + 1;
            }
        }

        // PROPAGACJA DAG: Bramka została wykonana, więc przeglądamy jej następców
        for (int next_id : node.next_gates) {
            // Zmniejszamy licznik zależności. Jeśli spadnie do zera, następca jest gotowy!
            if (--deps[next_id] == 0) front_layer.push_back(next_id);
        }
    }
}


// ==============================================================================
// sabre_route_step
// Gdy utknęliśmy (brak bramek w odległości 1), szukamy najlepszego swapu/ruchu.
// ==============================================================================
bool SabreDEHybridRouter::sabre_route_step(
    const HybridSabreDAG& current_dag,
    const std::vector<int>& front_layer,
    const std::vector<std::vector<int>>& lookahead_layers,
    SabreRunContext& ctx)
{
    int p_cnt = topo.qubits.size();

    std::vector<uint8_t> in_active(p_cnt, 0);
    std::vector<int> active_phys_qubits;
    active_phys_qubits.reserve(num_qubits);

    double base_f_score = 0.0;

    for (int gate_id : front_layer) {
        const auto& node = current_dag.nodes[gate_id];
        if (node.lq1 != -1) {
            int p1 = ctx.state.logical_to_phys_qubit[node.lq1];
            if (!in_active[p1]) { in_active[p1] = true; active_phys_qubits.push_back(p1); }
        }
        if (node.lq2 != -1) {
            int p2 = ctx.state.logical_to_phys_qubit[node.lq2];
            if (!in_active[p2]) { in_active[p2] = true; active_phys_qubits.push_back(p2); }

            base_f_score += get_distance(ctx.state.logical_to_phys_qubit[node.lq1], ctx.state.logical_to_phys_qubit[node.lq2]);
        }
    }

    std::vector<double> base_ext_scores(lookahead_layers.size(), 0.0);
    for (size_t d = 0; d < lookahead_layers.size(); ++d) {
        for (int gate_id : lookahead_layers[d]) {
            const auto& n = current_dag.nodes[gate_id];
            if (n.lq2 != -1) {
                base_ext_scores[d] += get_distance(ctx.state.logical_to_phys_qubit[n.lq1], ctx.state.logical_to_phys_qubit[n.lq2]);
            }
        }
    }

    struct SwapCandidate { int p1, p2, r_shared; double score; };

    // [OPTYMALIZACJA] O(1) dla weryfikacji powielonych krawędzi
    std::vector<uint8_t> seen_edges(p_cnt * p_cnt, 0);

    // [OPTYMALIZACJA] Zamiast składować wszystko do wektora, wyznaczamy zwycięzcę w locie
    double best_score = std::numeric_limits<double>::max();
    SwapCandidate best_candidate = {-1, -1, -1, best_score};

    for (int p : active_phys_qubits) {
        for (int r : topo.qubits[p].connected_resonators) {
            for (int neighbor_p : topo.resonators[r].connected_qubits) {
                if (p == neighbor_p) continue;

                int u = std::min(p, neighbor_p);
                int v = std::max(p, neighbor_p);
                int edge_idx = u * p_cnt + v;

                if (seen_edges[edge_idx]) continue;
                seen_edges[edge_idx] = true;

                double delta_f = 0.0;
                for (int gate_id : front_layer) {
                    const auto& n = current_dag.nodes[gate_id];
                    if (n.lq2 != -1) {
                        int lq1_p = ctx.state.logical_to_phys_qubit[n.lq1];
                        int lq2_p = ctx.state.logical_to_phys_qubit[n.lq2];

                        if (lq1_p == u || lq1_p == v || lq2_p == u || lq2_p == v) {
                            int new_l1 = (lq1_p == u) ? v : ((lq1_p == v) ? u : lq1_p);
                            int new_l2 = (lq2_p == u) ? v : ((lq2_p == v) ? u : lq2_p);
                            delta_f += (get_distance(new_l1, new_l2) - get_distance(lq1_p, lq2_p));
                        }
                    }
                }

                double sF = base_f_score + delta_f;
                double f_norm = front_layer.empty() ? 0 : (sF / front_layer.size());

                double ext_score = 0.0;
                double layer_weight = EXTENDED_LAYER_WEIGHT;

                for (size_t d = 0; d < lookahead_layers.size(); ++d) {
                    if (lookahead_layers[d].empty()) {
                        layer_weight *= 0.5;
                        continue;
                    }

                    double delta_ext = 0.0;
                    for (int gate_id : lookahead_layers[d]) {
                        const auto& n = current_dag.nodes[gate_id];
                        if (n.lq2 != -1) {
                            int lq1_p = ctx.state.logical_to_phys_qubit[n.lq1];
                            int lq2_p = ctx.state.logical_to_phys_qubit[n.lq2];

                            if (lq1_p == u || lq1_p == v || lq2_p == u || lq2_p == v) {
                                int new_l1 = (lq1_p == u) ? v : ((lq1_p == v) ? u : lq1_p);
                                int new_l2 = (lq2_p == u) ? v : ((lq2_p == v) ? u : lq2_p);
                                delta_ext += (get_distance(new_l1, new_l2) - get_distance(lq1_p, lq2_p));
                            }
                        }
                    }

                    double layer_score = base_ext_scores[d] + delta_ext;
                    ext_score += layer_weight * (layer_score / lookahead_layers[d].size());
                    layer_weight *= 0.5;
                }

                double score = (f_norm + ext_score) * fast_max_hyb(ctx.decay[u], ctx.decay[v]);

                if (score < best_score) {
                    best_score = score;
                    best_candidate = {u, v, r, score};
                }
            }
        }
    }

    if (best_candidate.p1 != -1) {
        int p1 = best_candidate.p1, p2 = best_candidate.p2, r_step = best_candidate.r_shared;

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

inline void SabreDEHybridRouter::extract_sabre_executables(const HybridSabreDAG& current_dag, std::vector<int>& front_layer, std::vector<int>& executable_gates, SabreRunContext& ctx) {
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
            // Wyciągamy na wzór oryginału by zachować 100% zgodności kolejności:
            front_layer[i] = front_layer.back();
            front_layer.pop_back();
        } else {
            ++i;
        }
    }
}

inline void SabreDEHybridRouter::build_sabre_lookahead(const HybridSabreDAG& current_dag, const std::vector<int>& front_layer, std::vector<std::vector<int>>& lookahead_layers, std::vector<bool>& in_extended) {
    std::fill(in_extended.begin(), in_extended.end(), false);
    std::vector<int> current_level = front_layer;

    for (size_t d = 0; d < lookahead_layers.size(); ++d) {
        std::vector<int> next_level;
        for (int gate_id : current_level) {
            for (int nxt : current_dag.nodes[gate_id].next_gates) {
                if (!in_extended[nxt]) {
                    in_extended[nxt] = true;
                    next_level.push_back(nxt);
                }
            }
        }
        lookahead_layers[d] = next_level;
        current_level = next_level;

        if (current_level.empty()) break;
    }
}

// ==============================================================================
// GŁÓWNY PRZEBIEG HEURYSTYKI SABRE (DO PRZODU LUB WSTECZ)
// ==============================================================================
std::vector<int> SabreDEHybridRouter::run_sabre_pass(const HybridSabreDAG& current_dag, const std::vector<int>& initial_l2p, bool is_final, std::vector<int>* out_logical_order) {
    int n_gates = current_dag.nodes.size();
    int p_cnt = topo.qubits.size();

    SabreRunContext ctx;
    ctx.state.logical_to_phys_qubit = initial_l2p;
    ctx.state.phys_qubit_to_logical.assign(p_cnt, -1);

    for (int i = 0; i < num_qubits; ++i) {
        if (initial_l2p[i] != -1) ctx.state.phys_qubit_to_logical[initial_l2p[i]] = i;
    }

    ctx.state.q_on_resonator.assign(p_cnt, -1);
    ctx.state.res_occupant.assign(topo.resonators.size(), -1);
    ctx.state.q_ready_time.assign(p_cnt, 0.0);
    ctx.state.r_ready_time.assign(topo.resonators.size(), 0.0);
    ctx.q_logical_depth.assign(p_cnt, 0);
    ctx.decay.assign(p_cnt, 1.0);

    std::vector<int> front_layer;
    front_layer.reserve(num_qubits);
    std::vector<int> deps = current_dag.initial_deps;

    for (int i = 0; i < n_gates; ++i) {
        if (deps[i] == 0) front_layer.push_back(i);
    }

    int gates_routed = 0;
    std::vector<int> executable_gates;
    executable_gates.reserve(num_qubits);

    const int LOOKAHEAD_DEPTH = 20;
    std::vector<std::vector<int>> lookahead_layers(LOOKAHEAD_DEPTH);
    std::vector<bool> in_extended(n_gates, false);

    int debug_iterations = 0;
    const int MAX_ITERATIONS = n_gates * 100;

    while (gates_routed < n_gates) {
        if (++debug_iterations > MAX_ITERATIONS) break;

        extract_sabre_executables(current_dag, front_layer, executable_gates, ctx);

        if (!executable_gates.empty()) {
            sabre_execute_gates(current_dag, executable_gates, front_layer, deps, ctx, out_logical_order);
            gates_routed += executable_gates.size();
            std::fill(ctx.decay.begin(), ctx.decay.end(), 1.0);
            continue;
        }

        build_sabre_lookahead(current_dag, front_layer, lookahead_layers, in_extended);
        sabre_route_step(current_dag, front_layer, lookahead_layers, ctx);
    }

    if (is_final) {
        for (size_t p = 0; p < p_cnt; ++p) sabre_force_qubit_down(p, -1, ctx);

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

// ==============================================================================
// "Look-ahead/Look-back" pass
// ==============================================================================
inline std::vector<int> create_identity_mapping(int num_qubits) {
    std::vector<int> mapping(num_qubits);
    std::iota(mapping.begin(), mapping.end(), 0); // O(N) szybkie wypełnienie 0, 1, 2...
    return mapping;
}

inline void cache_initial_mapping(const std::vector<int>& pre_final_mapping, std::vector<int>& cache_dest) {
    std::fill(cache_dest.begin(), cache_dest.end(), -1);
    std::copy(pre_final_mapping.begin(), pre_final_mapping.end(), cache_dest.begin());
}

// --- GŁÓWNA FUNKCJA ---

std::vector<int> SabreDEHybridRouter::prepare_sabre_seed() {
    // 1. Inicjalizacja tożsamościowa
    std::vector<int> initial_mapping = create_identity_mapping(num_qubits);

    saved_initial_l2p.assign(topo.qubits.size(), -1);
    std::copy(initial_mapping.begin(), initial_mapping.end(), saved_initial_l2p.begin());

    // 2. Przebieg do przodu (Forward)
    auto post_fwd = run_sabre_pass(dag, initial_mapping, false, nullptr);

    // 3. Przebieg wsteczny (Backward)
    auto pre_final = run_sabre_pass(reversed_dag, post_fwd, false, nullptr);

    // 4. Buforowanie wyznaczonego mapowania
    cache_initial_mapping(pre_final, saved_initial_l2p);

    // 5. Przebieg finalny: Generuje faktyczną kolejność wykonania bramek
    std::vector<int> logical_exec_order;
    run_sabre_pass(dag, pre_final, true, &logical_exec_order);

    return logical_exec_order;
}

bool SabreDEHybridRouter::validate_logical_circuit(const std::vector<RoutingRecord>& log) const {
    int n_gates = dag.nodes.size();

    // 1. Mapowanie oryginalnych ID bramek na indeksy węzłów w DAG
    std::unordered_map<int, int> orig_to_dag;
    std::unordered_map<int, int> orig_1q_to_dag;
    for (int i = 0; i < n_gates; ++i) {
        if (dag.nodes[i].orig_2q_id != -1) orig_to_dag[dag.nodes[i].orig_2q_id] = i;
        for (int q1_id : dag.nodes[i].q1_1q_ids) orig_1q_to_dag[q1_id] = i;
        for (int q2_id : dag.nodes[i].q2_1q_ids) orig_1q_to_dag[q2_id] = i;
    }

    // LOKALNA PRZYCZYNOWOŚĆ DLA KUBITÓW
    std::vector<int> prev_node_lq1(n_gates, -1);
    std::vector<int> prev_node_lq2(n_gates, -1);
    std::vector<int> last_seen(num_qubits, -1);

    for (int i = 0; i < n_gates; ++i) {
        if (dag.nodes[i].lq1 != -1) {
            prev_node_lq1[i] = last_seen[dag.nodes[i].lq1];
            last_seen[dag.nodes[i].lq1] = i;
        }
        if (dag.nodes[i].lq2 != -1) {
            prev_node_lq2[i] = last_seen[dag.nodes[i].lq2];
            last_seen[dag.nodes[i].lq2] = i;
        }
    }

    std::vector<int> latest_completed_node(num_qubits, -1);
    std::vector<int> executed_count(n_gates, 0);
    std::unordered_set<int> executed_1q_gates;
    std::vector<int> current_deps = dag.initial_deps;
    std::vector<int> l2p = saved_initial_l2p;

    // Struktury śledzące czas na obu warstwach sprzętowych
    std::vector<double> qubit_ready_time(topo.qubits.size(), 0.0);
    std::vector<double> resonator_ready_time(topo.resonators.size(), 0.0);

    // --- NOWA STRUKTURA: Śledzenie obecności na rezonatorze ---
    std::vector<int> val_p_on_res(topo.qubits.size(), -1);

    const double EPS = 1e-5;
    int step = 0;

    for (const auto& record : log) {
        step++;
        double op_duration = 0.0;
        std::vector<int> active_qubits;
        std::vector<int> active_resonators;

        if (record.op_type == OpType::SWAP) {
            op_duration = SWAP_COST;
            int lq1 = record.lq1, lq2 = record.lq2;

            if (lq1 != -1 && l2p[lq1] != record.loc1) return false;
            if (lq2 != -1 && l2p[lq2] != record.loc2) return false;

            if (lq1 != -1) l2p[lq1] = record.loc2;
            if (lq2 != -1) l2p[lq2] = record.loc1;

            if (record.loc1 != -1) {
                active_qubits.push_back(record.loc1);
                val_p_on_res[record.loc1] = -1; // Gwarancja, że zeszły z rezonatorów
            }
            if (record.loc2 != -1) {
                active_qubits.push_back(record.loc2);
                val_p_on_res[record.loc2] = -1;
            }
            if (record.res_used != -1) active_resonators.push_back(record.res_used);
        }
        else if (record.op_type == OpType::MOVE) {
            op_duration = MOVE_COST;

            bool moving_from_res = (record.loc1 >= RES_OFFSET);
            bool moving_to_res = (record.loc2 >= RES_OFFSET);

            if (moving_to_res) {
                // Wejście na rezonator
                int phys_src = record.loc1;
                int res_dst = record.loc2 - RES_OFFSET;

                if (phys_src != -1) {
                    active_qubits.push_back(phys_src);
                    val_p_on_res[phys_src] = res_dst; // Zapamiętanie, gdzie wszedł
                }
                active_resonators.push_back(res_dst);
            }
            else if (moving_from_res) {
                // Zejście z rezonatora
                int res_src = record.loc1 - RES_OFFSET;
                int phys_dst = record.loc2;

                active_resonators.push_back(res_src);
                if (phys_dst != -1) {
                    active_qubits.push_back(phys_dst);
                    val_p_on_res[phys_dst] = -1; // Zwolnienie flagi obecności
                }
                if (record.lq1 != -1) l2p[record.lq1] = phys_dst;
            }
        }
        else if (record.op_type == OpType::GATE_1Q) {
            op_duration = 1.0;
            int lq1 = record.lq1;
            int p1 = record.loc1;
            int orig_id = record.gate_id;

            if (executed_1q_gates.count(orig_id)) {
                std::cerr << "Błąd 1Q [Krok " << step << "]: Bramka (ID: " << orig_id << ") wykonana ponownie!\n";
                return false;
            }

            auto it1q = orig_1q_to_dag.find(orig_id);
            if (it1q != orig_1q_to_dag.end()) {
                int parent_gid = it1q->second;

                if (executed_count[parent_gid] > 0) {
                    std::cerr << "Błąd Przyczynowości [Krok " << step << "]: Bramka 1Q wykonana za późno!\n";
                    return false;
                }

                int expected_prev = (dag.nodes[parent_gid].lq1 == lq1) ? prev_node_lq1[parent_gid] : prev_node_lq2[parent_gid];
                if (expected_prev != -1 && latest_completed_node[lq1] != expected_prev) {
                    std::cerr << "Błąd Przyczynowości [Krok " << step << "]: Naruszenie czasu lokalnego!\n";
                    return false;
                }
            }

            executed_1q_gates.insert(orig_id);

            // --- POPRAWKA: Blokada wykonywania 1Q na rezonatorze ---
            if (p1 >= RES_OFFSET) {
                std::cerr << "Błąd Sprzętowy 1Q [Krok " << step << "]: Próba wykonania na rezonatorze (" << p1 << ")!\n";
                return false;
            }
            if (lq1 != -1 && l2p[lq1] != p1) {
                std::cerr << "Błąd Mapowania 1Q [Krok " << step << "]: Błędna lokalizacja!\n";
                return false;
            }
            if (p1 != -1 && val_p_on_res[p1] != -1) {
                std::cerr << "Błąd Sprzętowy 1Q [Krok " << step << "]: Kubit fizyczny " << p1
                          << " musi zejść z rezonatora " << val_p_on_res[p1] << " przed wykonaniem 1Q!\n";
                return false;
            }

            if (p1 != -1) active_qubits.push_back(p1);

            // Zamykanie wiszącego węzła
            if (it1q != orig_1q_to_dag.end()) {
                int parent_gid = it1q->second;
                const auto& node = dag.nodes[parent_gid];
                if (node.lq2 == -1) {
                    bool all_1q_done = true;
                    for (int q1_id : node.q1_1q_ids) if (!executed_1q_gates.count(q1_id)) all_1q_done = false;

                    if (all_1q_done && executed_count[parent_gid] == 0) {
                        executed_count[parent_gid]++;
                        for (int nxt : node.next_gates) current_deps[nxt]--;
                        latest_completed_node[node.lq1] = parent_gid;
                    }
                }
            }
        }
        else if (record.op_type == OpType::GATE_2Q_LOCAL) {
            op_duration = 1.0;
            auto it = orig_to_dag.find(record.gate_id);
            if (it == orig_to_dag.end()) continue;

            int gid = it->second;

            if (executed_count[gid] > 0 || current_deps[gid] > 0) return false;

            const auto& node = dag.nodes[gid];
            bool lq_match = (node.lq1 == record.lq1 && node.lq2 == record.lq2) || (node.lq1 == record.lq2 && node.lq2 == record.lq1);
            if (!lq_match) return false;

            int p1 = l2p[record.lq1];
            int p2 = l2p[record.lq2];
            int res = record.loc1 - RES_OFFSET;

            if (p1 == -1 || p2 == -1) return false;
            if (record.loc2 != p1 && record.loc2 != p2) return false;

            // --- POPRAWKA: Weryfikacja Topologiczna wybranego Rezonatora ---
            const auto& p1_res = topo.qubits[p1].connected_resonators;
            const auto& p2_res = topo.qubits[p2].connected_resonators;
            bool valid_topology = (std::find(p1_res.begin(), p1_res.end(), res) != p1_res.end()) &&
                                  (std::find(p2_res.begin(), p2_res.end(), res) != p2_res.end());

            if (!valid_topology) {
                std::cerr << "Błąd Topologii 2Q [Krok " << step << "]: Rezonator " << res
                          << " nie łączy kubitów " << p1 << " i " << p2 << "!\n";
                return false;
            }

            // --- POPRAWKA: Weryfikacja faktycznego On-boarding'u ---
            if (val_p_on_res[p1] != res && val_p_on_res[p2] != res) {
                std::cerr << "Błąd Sprzętowy 2Q [Krok " << step << "]: Żaden z kubitów ("
                          << p1 << ", " << p2 << ") nie zszedł na rezonator " << res << " przed uderzeniem bramki!\n";
                return false;
            }
            if (val_p_on_res[p1] == res && val_p_on_res[p2] == res) {
                std::cerr << "Błąd Sprzętowy 2Q [Krok " << step << "]: Kolizja! Oba kubity sprzęgnięte do rezonatora " << res << " naraz!\n";
                return false;
            }

            executed_count[gid]++;
            for (int nxt : node.next_gates) current_deps[nxt]--;

            latest_completed_node[node.lq1] = gid;
            latest_completed_node[node.lq2] = gid;

            active_qubits.push_back(p1);
            active_qubits.push_back(p2);
            active_resonators.push_back(res);
        }

        // --- WALIDACJA CZASOWA MAKESPAN ---
        for (int q : active_qubits) {
            if (record.start_time < qubit_ready_time[q] - EPS) {
                std::cerr << "Błąd Czasowy [Krok " << step << "]: Kubit fizyczny " << q << " zderzył się w czasie!\n";
                return false;
            }
            qubit_ready_time[q] = record.start_time + op_duration;
        }

        for (int r : active_resonators) {
            if (record.start_time < resonator_ready_time[r] - EPS) {
                std::cerr << "Błąd Czasowy [Krok " << step << "]: Rezonator " << r << " zderzył się w czasie!\n";
                return false;
            }
            resonator_ready_time[r] = record.start_time + op_duration;
        }
    }

    // --- KOMPLETNOŚĆ ---
    for (int i = 0; i < n_gates; ++i) {
        const auto& node = dag.nodes[i];
        if (node.lq2 != -1 && executed_count[i] != 1) {
            std::cerr << "Błąd Kompletności: Węzeł 2Q (ID DAG: " << i << ") nie został wykonany!\n";
            return false;
        }
        for (int q1_id : node.q1_1q_ids) if (!executed_1q_gates.count(q1_id)) return false;
        for (int q2_id : node.q2_1q_ids) if (!executed_1q_gates.count(q2_id)) return false;
    }

    return true;
}

template <bool LOGGING>
bool SabreDEHybridRouter::de_naive_route_step(
    const HybridSabreDAG& current_dag,
    const std::vector<int>& front_layer,
    std::vector<double>& decay,
    DERunContext& ctx,
    const double* genes, int n_gates,
    std::vector<RoutingRecord>* out_log)
{
    // 1. Zbieranie aktywnych kubitów i lokalizacji
    std::fill(ctx.ws.in_active.begin(), ctx.ws.in_active.end(), 0);
    ctx.ws.active_phys.clear();
    ctx.ws.f_locs.clear();

    // --- [NOWE] Obliczanie KOSZTU BAZOWEGO ---
    double base_sF = 0.0;

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
        // Cache dla lokalizacji docelowych we front_layer + nabijanie BAZY
        if (node.lq2 != -1) {
            ctx.ws.f_locs.push_back({ctx.l2p[node.lq1], ctx.l2p[node.lq2]});
            base_sF += get_distance(ctx.l2p[node.lq1], ctx.l2p[node.lq2]);
        }
    }

    ctx.ws.route_step_counter++;
    int step_id = ctx.ws.route_step_counter;
    int p_cnt = topo.qubits.size();

    double best_score = std::numeric_limits<double>::max();
    int best_p1 = -1, best_p2 = -1, best_r = -1;

    // 2. Szukanie krawędzi przy użyciu DELTA SCORING O(k)
    for (int p : ctx.ws.active_phys) {
        // Wyciągamy zakres z płaskiego wektora w czasie O(1)
        int start_idx = precomputed_swaps_offsets[p];
        int end_idx = precomputed_swaps_offsets[p + 1];

        for (int i = start_idx; i < end_idx; ++i) {
            const auto& edge = flat_precomputed_swaps[i];
            int nb = edge.neighbor;
            int r = edge.resonator;

            int u = std::min(p, nb), v = std::max(p, nb);

            int edge_idx = u * p_cnt + v;
            if (ctx.ws.seen_edges[edge_idx] == step_id) continue;
            ctx.ws.seen_edges[edge_idx] = step_id;

            // Liczymy wyłącznie DELTĘ zmian
            double delta_sF = 0.0;
            for (const auto& loc : ctx.ws.f_locs) {
                int l1 = loc.first, l2 = loc.second;

                // Sprawdzamy tylko te bramki, których dotyka ten specyficzny SWAP
                if (l1 == u || l1 == v || l2 == u || l2 == v) {
                    int new_l1 = (l1 == u) ? v : ((l1 == v) ? u : l1);
                    int new_l2 = (l2 == u) ? v : ((l2 == v) ? u : l2);
                    delta_sF += (get_distance(new_l1, new_l2) - get_distance(l1, l2));
                }
            }

            // Ostateczny koszt to Baza + Delta
            double sF = base_sF + delta_sF;
            double f_norm = ctx.ws.f_locs.empty() ? 0 : (sF / ctx.ws.f_locs.size());

            // Mutacja wprowadzająca czynnik losowy z algorytmu DE
            double mutation = 0.95 + (genes[u % n_gates] * 0.1);

            // Wynik opiera się na f_norm i karach (decay)
            double score = f_norm * fast_max_hyb(decay[u], decay[v]) * mutation;

            if (score < best_score) {
                best_score = score;
                best_p1 = u; best_p2 = v; best_r = r;
            }
        }
    }

    // 3. Wykonanie ruchu na najlepszym kandydacie (część wykonawcza bez zmian)
    if (best_p1 != -1) {
        int p1 = best_p1, p2 = best_p2, r_step = best_r;

        de_force_down<LOGGING>(p1, -1, ctx, genes, n_gates, out_log, false);
        de_force_down<LOGGING>(p2, -1, ctx, genes, n_gates, out_log, false);
        if (ctx.res_occ[r_step] != -1) de_force_down<LOGGING>(ctx.res_occ[r_step], -1, ctx, genes, n_gates, out_log, false);

        double step_start = fast_max_hyb(ctx.q_time[p1], fast_max_hyb(ctx.r_time[r_step], ctx.q_time[p2]));
        bool is_move = (ctx.p2l[p2] == -1) || (ctx.p2l[p1] == -1);
        // Koszt podwójny, bo wykonujemy 2 skoki (Qubit -> Rezonator -> Qubit)
        double cost = is_move ? 2.0 * MOVE_COST : SWAP_COST;

        if constexpr (LOGGING) {
            if (!is_move) {
                // Dla klasycznego SWAP-a zakładamy, że to pojedyncza, współdzielona operacja
                out_log->push_back({OpType::SWAP, -1, ctx.p2l[p1], ctx.p2l[p2], p1, p2, r_step, step_start});
            } else {
                // POPRAWKA TOPOLOGII: Rozbijamy Lazy Off-boarding na dwa rzetelne ruchy fizyczne!
                int moving_p = (ctx.p2l[p1] != -1) ? p1 : p2;
                int empty_p =  (ctx.p2l[p1] != -1) ? p2 : p1;
                int logic_q = ctx.p2l[moving_p];

                // HOP 1: Z obecnego kubitu fizycznego wchodzimy NA rezonator
                out_log->push_back({OpType::MOVE, -1, logic_q, -1, moving_p, r_step + RES_OFFSET, -1, step_start});

                // HOP 2: Z rezonatora schodzimy na docelowy, wolny kubit (z uwzględnieniem opóźnienia MOVE_COST)
                out_log->push_back({OpType::MOVE, -1, logic_q, -1, r_step + RES_OFFSET, empty_p, -1, step_start + MOVE_COST});
            }
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

OptimizationResult SabreDEHybridRouter::evaluate_qiskit_routed_circuit(const std::vector<Gate>& qiskit_circuit) {
    int p_cnt = topo.qubits.size();
    int r_cnt = topo.resonators.size();

    // Rejestry stanu maszyny fizycznej (czyste, bo Qiskit operuje na ID fizycznych)
    std::vector<int> q_on_res(p_cnt, -1);
    std::vector<int> res_occ(r_cnt, -1);
    std::vector<double> q_time(p_cnt, 0.0);
    std::vector<double> r_time(r_cnt, 0.0);
    std::vector<int> q_depth(p_cnt, 0);

    int out_swaps = 0;
    int out_moves = 0;
    double total_routing_cost = 0.0;
    std::vector<RoutingRecord> final_log;

    // Lambda: Wymuszenie zejścia fizycznego kubitu z rezonatora (Automatyczny MOVE)
    auto force_down = [&](int p) {
        int r = q_on_res[p];
        if (r == -1) return;

        double start = fast_max_hyb(q_time[p], r_time[r]);
        final_log.push_back({OpType::MOVE, -1, -1, -1, r + RES_OFFSET, p, -1, start});

        out_moves++;
        total_routing_cost += MOVE_COST;

        q_time[p] = start + MOVE_COST;
        r_time[r] = start + MOVE_COST;
        q_depth[p]++;
        q_on_res[p] = -1;
        res_occ[r] = -1;
    };

    // PĘTLA SYMULACJI LINII CZASOWEJ
    for (const auto& g : qiskit_circuit) {
        // --- BRAMKA 1-KUBITOWA ---
        if (g.type == "1Q") {
            int p1 = g.target_q1;
            force_down(p1); // Musi zejść na węzeł, żeby wykonać 1Q

            final_log.push_back({OpType::GATE_1Q, g.id, -1, -1, p1, -1, -1, q_time[p1]});
            q_time[p1] += 1.0;
            q_depth[p1]++;
        }
        // --- BRAMKA 2-KUBITOWA (CX) ---
        else if (g.type == "2Q") {
            int p1 = g.target_q1;
            int p2 = g.target_q2;

            // Szukamy wspólnego rezonatora sprzętowego
            int chosen_res = -1;
            for (int r : topo.qubits[p1].connected_resonators) {
                const auto& p2_res = topo.qubits[p2].connected_resonators;
                if (std::find(p2_res.begin(), p2_res.end(), r) != p2_res.end()) {
                    chosen_res = r; break;
                }
            }

            if (chosen_res == -1) {
                std::cerr << "[BŁĄD WALIDACJI QISKITA] Próba wykonania 2Q na niesąsiadujących kubitach: " << p1 << " i " << p2 << "\n";
                continue;
            }

            // Decyzja czasowa - kto wchodzi, kto czeka na zewnątrz
            bool q1_goes_in = (q_time[p1] <= q_time[p2]);
            int p_in = q1_goes_in ? p1 : p2;
            int p_stay = q1_goes_in ? p2 : p1;

            force_down(p_stay);
            if (res_occ[chosen_res] != -1 && res_occ[chosen_res] != p_in) {
                force_down(res_occ[chosen_res]);
            }

            double start_time;
            if (q_on_res[p_in] != chosen_res) {
                force_down(p_in);
                start_time = fast_max_hyb(q_time[p_in], fast_max_hyb(q_time[p_stay], r_time[chosen_res]));

                final_log.push_back({OpType::MOVE, g.id, -1, -1, p_in, chosen_res + RES_OFFSET, -1, start_time});
                out_moves++;
                total_routing_cost += MOVE_COST;
                start_time += MOVE_COST;

                q_on_res[p_in] = chosen_res;
                res_occ[chosen_res] = p_in;
                q_depth[p_in]++;
            } else {
                start_time = fast_max_hyb(q_time[p_in], fast_max_hyb(q_time[p_stay], r_time[chosen_res]));
            }

            final_log.push_back({OpType::GATE_2Q_LOCAL, g.id, -1, -1, chosen_res + RES_OFFSET, p_stay, -1, start_time});
            q_time[p_in] = q_time[p_stay] = r_time[chosen_res] = start_time + 1.0;

            int max_d = fast_max_hyb(q_depth[p_in], q_depth[p_stay]);
            q_depth[p_in] = q_depth[p_stay] = max_d + 1;
        }
        // --- ABSTRACT SWAP (Wstawiony przez Qiskit SABRE) ---
        else if (g.type == "SWAP") {
            int p1 = g.target_q1;
            int p2 = g.target_q2;

            int chosen_res = -1;
            for (int r : topo.qubits[p1].connected_resonators) {
                const auto& p2_res = topo.qubits[p2].connected_resonators;
                if (std::find(p2_res.begin(), p2_res.end(), r) != p2_res.end()) {
                    chosen_res = r; break;
                }
            }

            // Wykonujemy fizyczny SWAP na architekturze magistrali rezonatorowej
            force_down(p1);
            force_down(p2);
            if (chosen_res != -1 && res_occ[chosen_res] != -1) force_down(res_occ[chosen_res]);

            // Dodatkowe zabezpieczenie: jeśli nie połączone rezonatorem, start time max dla kubitów
            double r_time_val = (chosen_res != -1) ? r_time[chosen_res] : 0.0;
            double start_time = fast_max_hyb(q_time[p1], fast_max_hyb(q_time[p2], r_time_val));

            final_log.push_back({OpType::SWAP, -1, -1, -1, p1, p2, chosen_res, start_time});
            out_swaps++;
            total_routing_cost += SWAP_COST;

            q_time[p1] = q_time[p2] = start_time + SWAP_COST;
            if (chosen_res != -1) r_time[chosen_res] = start_time + SWAP_COST;

            int max_d = fast_max_hyb(q_depth[p1], q_depth[p2]) + 1;
            q_depth[p1] = q_depth[p2] = max_d;
        }
    }

    // Czyszczenie końcowe - ściągamy wszystkie stany na węzły fizyczne
    for (int p = 0; p < p_cnt; ++p) force_down(p);

    double max_physical_time = *std::max_element(q_time.begin(), q_time.end());
    int max_logical_depth = *std::max_element(q_depth.begin(), q_depth.end());

    // Liczymy identyczny koszt fitness w celu bezpośredniego porównania
    double final_score = (total_routing_cost * 100.0) + (out_swaps * 500.0) + (max_physical_time * 0.1) + (max_logical_depth * 0.05);

    std::cout << "\n============================================================\n";
    std::cout << "=== FINALNA EWALUACJA PROFILU QISKIT SABRE NA MODELU IQM ===\n";

    // UWAGA: Upewnij się, że masz metodę print_statistics dostępną w swojej klasie.
    // print_statistics(final_score, max_physical_time, max_logical_depth, out_swaps, out_moves, 0.0);

    std::cout << "KOSZT: " << final_score << " | CZAS FIZYCZNY: " << max_physical_time << " | GŁĘBOKOŚĆ: " << max_logical_depth << "\n";
    std::cout << "SWAPY: " << out_swaps << " | MOVES: " << out_moves << "\n";
    std::cout << "============================================================\n";

    // export_moves_to_file(final_log, "../qiskit_faithful_schedule_" + name + ".txt");

    return {final_score, max_physical_time, max_logical_depth, out_swaps, out_moves};
}