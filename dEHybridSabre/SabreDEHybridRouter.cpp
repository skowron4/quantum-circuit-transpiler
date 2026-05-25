#include "SabreDEHybridRouter.h"
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include <fstream>
#include <sstream>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

inline double fast_max_hyb(double a, double b) { return (a > b) ? a : b; }

SabreDEHybridRouter::SabreDEHybridRouter(const BaseTopology& topology, const std::vector<Gate>& circuit, int qubits_cnt, std::string name, int max_look_ahead_depth)
    : topo(topology), original_circuit(circuit), num_qubits(qubits_cnt), name(std::move(name)), max_look_ahead_depth(max_look_ahead_depth)  {
    int p_cnt = topo.qubits.size();
    int r_cnt = topo.resonators.size();

    // Spłaszczony bufor cache dla odległości (p_cnt x p_cnt)
    dist_cache.assign(p_cnt * p_cnt, -1);

    dag = build_dag(original_circuit);
    int g_cnt = dag.nodes.size(); // <--- Tutaj definiujesz liczbę bramek jako g_cnt

    std::vector<Gate> reversed_circuit = original_circuit;
    std::reverse(reversed_circuit.begin(), reversed_circuit.end());
    reversed_dag = build_dag(reversed_circuit);

    // [OPTYMALIZACJA 4] Pełne resize() dla struktury workspace w celu wykorzystania std::fill w evaluate_fitness_impl
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

    // POPRAWKA: Użycie 'g_cnt' zamiast 'n_gates'
    workspace.visited_epoch.assign(g_cnt, 0);
    workspace.current_bfs_epoch = 0;

    workspace.decay.resize(p_cnt, 1.0);
    workspace.deps.resize(g_cnt, 0);

    workspace.shared_res.reserve(r_cnt);
    workspace.in_active.resize(p_cnt, 0);
    workspace.active_phys.reserve(p_cnt);

    // [OPTYMALIZACJA 1] Pre-alokacja dla seen_edges
    workspace.route_step_counter = 0;
    workspace.seen_edges.assign(p_cnt * p_cnt, 0);

    // [OPTYMALIZACJA 2] Pre-alokacja do hoistingu
    workspace.f_locs.reserve(g_cnt);

    // [SUPER-OPTYMALIZACJA] Płaskie struktury dla Cache Hit Rate (Flattening)
    // Zapas 5 offsetów bezpiecznie obsługuje MAX_LOOKAHEAD_DEPTH (np. 3 lub 4)
    workspace.lookahead_offsets.reserve(5);
    workspace.lookahead_locs_offsets.reserve(5);

    // Alokujemy flat-bufory na g_cnt (maksymalna możliwa liczba bramek w lookahead).
    // Daje to 100% pewność, że std::vector nigdy nie wykona malloc/free podczas oceny osobników DE.
    workspace.flat_lookahead_layers.reserve(g_cnt);
    workspace.flat_lookahead_locs.reserve(g_cnt);
}

// ==============================================================================
// OBLICZANIE ODLEGŁOŚCI Z CACHE
// ==============================================================================
int SabreDEHybridRouter::get_distance(int p1, int p2) const {
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

// ==============================================================================
// BUDOWANIE SKIEROWANEGO GRAFU ACYKLICZNEGO (DAG)
// ==============================================================================
HybridSabreDAG SabreDEHybridRouter::build_dag(const std::vector<Gate>& circuit) const {
    HybridSabreDAG local_dag;
    // last_gate - pamięta ID ostatniego węzła grafu, który operował na danym kubicie
    // pending_1q - zlicza "w locie" bramki 1-kubitowe, żeby nie robić dla nich osobnych węzłów
    std::vector<int> last_gate(num_qubits, -1), pending_1q(num_qubits, 0);
    int group_id = 0;

    for (const auto& g : circuit) {
        if (g.type == "1Q") {
            // Bramka 1-kubitowa? Super, po prostu podbij licznik dla tego kubitu,
            // zbijemy je w jeden blok przy następnej bramce 2-kubitowej.
            pending_1q[g.target_q1]++;
        }
        else {
            // Natrafiliśmy na bramkę 2-kubitową. Tworzymy nowy węzeł.
            // Dołączamy do niego od razu skumulowane zliczenia bramek 1Q z obu kubitów.
            HybridGateNode node{group_id++, g.target_q1, g.target_q2, pending_1q[g.target_q1], pending_1q[g.target_q2], 0, {}};

            // Zerujemy liczniki 1Q, bo właśnie je zużyliśmy
            pending_1q[g.target_q1] = 0; pending_1q[g.target_q2] = 0;

            // Sprawdzamy zależności czasowe (krawędzie DAGa).
            // Jeśli coś wcześniej operowało na kubicie nr 1 tej bramki (lq1):
            if (last_gate[node.lq1] != -1) {
                // Dodaj obecną bramkę do listy "następnych" dla tej starej bramki
                local_dag.nodes[last_gate[node.lq1]].next_gates.push_back(node.id);
                // Ta nowa bramka musi poczekać na tamtą (zwiększamy licznik zależności)
                node.deps_count++;
            }
            // To samo dla kubitu nr 2 (lq2), ale upewniamy się, że nie liczymy dwa razy
            // tego samego węzła, jeśli stara bramka też operowała na tych samych dwóch kubitach.
            if (last_gate[node.lq2] != -1 && last_gate[node.lq2] != last_gate[node.lq1]) {
                local_dag.nodes[last_gate[node.lq2]].next_gates.push_back(node.id);
                node.deps_count++;
            }

            // Wrzucamy gotowy węzeł do grafu i aktualizujemy "kto ostatni miał te kubity"
            local_dag.nodes.push_back(node);
            last_gate[node.lq1] = node.id; last_gate[node.lq2] = node.id;
        }
    }

    // Czyszczenie "resztek": jeśli obwód się skończył, a zostały nam na końcu
    // jakieś bramki 1Q, musimy utworzyć z nich sztuczne węzły na końcu DAGa.
    for (int i = 0; i < num_qubits; ++i) {
        if (pending_1q[i] > 0) {
            HybridGateNode node{group_id++, i, -1, pending_1q[i], 0, 0, {}};
            if (last_gate[i] != -1) {
                local_dag.nodes[last_gate[i]].next_gates.push_back(node.id);
                node.deps_count++;
            }
            local_dag.nodes.push_back(node);
            last_gate[i] = node.id;
            pending_1q[i] = 0;
        }
    }

    // Tworzymy kopię początkowych zależności (deps_count) dla każdego węzła,
    // żeby łatwo było resetować graf przy kolejnych ewaluacjach.
    local_dag.initial_deps.resize(local_dag.nodes.size());
    for (size_t i = 0; i < local_dag.nodes.size(); ++i) {
        local_dag.initial_deps[i] = local_dag.nodes[i].deps_count;
    }

    return local_dag;
}

// ==============================================================================
// Ściąganie kubitu z rezonatora z DE (podobnie jak w sabre - funkcja zoptymalizowana)
// ==============================================================================
template <bool LOGGING>
void SabreDEHybridRouter::de_force_down(int& p, int gate_id, DERunContext& ctx, const double* genes, int n_gates, std::vector<RoutingRecord>* out_log, bool can_lazy) {
    // Sprawdzamy, czy dany kubit fizyczny aktualnie "siedzi" na rezonatorze
    int r = ctx.q_on_res[p];
    if (r != -1) {
        int target_p = p; // Domyślnie zrzucamy go z powrotem tam, skąd przyszedł

        // Tryb "Lazy Off-boarding": Jeśli wolno nam zrzucić kubit w inne miejsce i nie jest to koniec obwodu
        if (can_lazy && gate_id >= 0) {
            int valid_count = 1;
            const auto& conn = topo.resonators[r].connected_qubits;

            // Szukamy sąsiadujących kubitów, które są aktualnie puste (nie mają mapowania logicznego)
            for (int c : conn) if (c != p && ctx.p2l[c] == -1 && ctx.q_on_res[c] == -1) valid_count++;

            // Jeśli jest więcej niż jedna opcja zrzutu, używamy genów ewolucji, aby podjąć decyzję
            if (valid_count > 1) {
                // Skalujemy gen (0.0-1.0) do indeksu dostępnego miejsca
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

        // Obliczamy czas startu tej operacji zrzutu (MOVE) biorąc pod uwagę zajętość zasobów
        double start = fast_max_hyb(ctx.q_time[p], ctx.r_time[r]);
        if (target_p != p) start = fast_max_hyb(start, ctx.q_time[target_p]);

        // Jeśli włączone jest logowanie, zapisujemy operację do historii
        if constexpr (LOGGING) out_log->push_back({OpType::MOVE, gate_id, ctx.p2l[p], -1, r + RES_OFFSET, target_p, -1, start});

        // Zwiększamy liczniki kar
        ctx.out_moves++;
        ctx.total_routing_cost += MOVE_COST;

        // Aktualizacja czasu i głębokości (dodanie kosztu MOVE)
        ctx.q_time[target_p] = start + MOVE_COST;
        ctx.r_time[r] = start + MOVE_COST;

        // Jeśli przenieśliśmy się na nowy fizyczny kubit, aktualizujemy mapowanie systemowe (l2p / p2l)
        if (target_p != p) {
            ctx.p2l[target_p] = ctx.p2l[p];
            ctx.l2p[ctx.p2l[p]] = target_p;
            ctx.p2l[p] = -1;
            ctx.q_depth[target_p] = ctx.q_depth[p] + 1;
            ctx.q_depth[p] = 0;
        } else {
            ctx.q_depth[p] += 1;
        }

        // Wyczyszczenie informacji o tym, że kubit był na rezonatorze i aktualizacja referencji
        ctx.q_on_res[p] = -1;
        ctx.q_on_res[target_p] = -1;
        ctx.res_occ[r] = -1;
        p = target_p; // Aktualizacja przekazanej przez referencję zmiennej
    }
}

// ==============================================================================
// Wykonywanie bramek w zależności od piorytetu wygenerowanego przez DE
// ==============================================================================
template <bool LOGGING>
void SabreDEHybridRouter::de_execute_gates(const HybridSabreDAG& current_dag, std::vector<int>& exec_gates, std::vector<int>& front_layer, std::vector<int>& deps, DERunContext& ctx, const double* genes, int n_gates, std::vector<RoutingRecord>* out_log, bool can_lazy) {
    // KLUCZOWY ELEMENT DE: Sortujemy bramki możliwe do wykonania według priorytetów zapisanych w genach.
    // Zamiast brać pierwszą z brzegu, algorytm uczy się optymalnej kolejności.
    std::sort(exec_gates.begin(), exec_gates.end(), [&](int a, int b) { return genes[a] > genes[b]; });

    for (int gate_id : exec_gates) {
        const auto& node = current_dag.nodes[gate_id];
        int p1 = ctx.l2p[node.lq1];

        // Przypadek 1: Bramka jednokubitowa (lq2 == -1)
        if (node.lq2 == -1) {
            if (node.q1_1q_count > 0) {
                // Jeśli kubit jest na rezonatorze, musimy go najpierw "ściągnąć" przed zaaplikowaniem bramki
                de_force_down<LOGGING>(p1, node.id, ctx, genes, n_gates, out_log, can_lazy);

                // Zapisujemy bramki do logów i zwiększamy czasy kubitu
                if constexpr (LOGGING) {
                    for (int k = 0; k < node.q1_1q_count; ++k) out_log->push_back({OpType::GATE_1Q, node.id, node.lq1, -1, p1, -1, -1, ctx.q_time[p1] + k});
                }
                ctx.q_time[p1] += node.q1_1q_count;
                ctx.q_depth[p1] += node.q1_1q_count;
            }
        }
        // Przypadek 2: Bramka dwukubitowa
        else {
            int p2 = ctx.l2p[node.lq2];

            // Ściągamy oba kubity z rezonatorów (jeśli tam były) do obsługi ewentualnych złączonych bramek 1Q
            if (node.q1_1q_count > 0) de_force_down<LOGGING>(p1, node.id, ctx, genes, n_gates, out_log, can_lazy);
            if (node.q2_1q_count > 0) de_force_down<LOGGING>(p2, node.id, ctx, genes, n_gates, out_log, can_lazy);

            // Logowanie i aktualizacja czasu dla wykonanych bramek 1-kubitowych
            if constexpr (LOGGING) {
                for (int k = 0; k < node.q1_1q_count; ++k) out_log->push_back({OpType::GATE_1Q, node.id, node.lq1, -1, p1, -1, -1, ctx.q_time[p1] + k});
                for (int k = 0; k < node.q2_1q_count; ++k) out_log->push_back({OpType::GATE_1Q, node.id, node.lq2, -1, p2, -1, -1, ctx.q_time[p2] + k});
            }
            ctx.q_time[p1] += node.q1_1q_count; ctx.q_time[p2] += node.q2_1q_count;
            ctx.q_depth[p1] += node.q1_1q_count; ctx.q_depth[p2] += node.q2_1q_count;

            // Szukanie wspólnego rezonatora łączącego p1 i p2
            ctx.ws.shared_res.clear();
            const auto& p2_res = topo.qubits[p2].connected_resonators;
            for (int r : topo.qubits[p1].connected_resonators) {
                if (std::find(p2_res.begin(), p2_res.end(), r) != p2_res.end()) ctx.ws.shared_res.push_back(r);
            }

            // Jeśli kubity mają wspólny rezonator, wykonujemy bramkę 2Q
            if (!ctx.ws.shared_res.empty()) {
                // Wykorzystujemy inny blok genów do wybrania konkretnego rezonatora, jeśli jest ich kilka
                int chosen_res = ctx.ws.shared_res[std::min(static_cast<int>(genes[n_gates + node.id] * ctx.ws.shared_res.size()), (int)ctx.ws.shared_res.size() - 1)];

                // Decydujemy, który kubit wchodzi na rezonator na podstawie zaległości czasowych lub genów
                bool q1_in = (node.q1_1q_count < node.q2_1q_count) ? true : (node.q1_1q_count > node.q2_1q_count) ? false : (genes[n_gates * 2 + node.id] < 0.5);

                int p_in = q1_in ? p1 : p2;
                int p_stay = q1_in ? p2 : p1; // Ten zostaje w miejscu, bramka uderza z rezonatora w niego

                // Upewniamy się, że kubit "stacjonarny" nie siedzi na rezonatorze
                de_force_down<LOGGING>(p_stay, node.id, ctx, genes, n_gates, out_log, can_lazy);

                // Zwalniamy rezonator docelowy, jeśli ktoś inny na nim siedzi
                if (ctx.res_occ[chosen_res] != -1 && ctx.res_occ[chosen_res] != p_in) {
                    int r_occ = ctx.res_occ[chosen_res];
                    de_force_down<LOGGING>(r_occ, node.id, ctx, genes, n_gates, out_log, can_lazy);
                }

                double start;

                // Wrzucamy kubit docelowy (p_in) na rezonator, o ile już tam nie jest
                if (ctx.q_on_res[p_in] != chosen_res) {
                    de_force_down<LOGGING>(p_in, node.id, ctx, genes, n_gates, out_log, can_lazy); // Zrzucamy go z ewentualnego innego rezonatora
                    start = fast_max_hyb(ctx.q_time[p_in], fast_max_hyb(ctx.q_time[p_stay], ctx.r_time[chosen_res]));

                    if constexpr (LOGGING) out_log->push_back({OpType::MOVE, node.id, ctx.p2l[p_in], -1, p_in, chosen_res + RES_OFFSET, -1, start});

                    ctx.out_moves++;
                    ctx.total_routing_cost += MOVE_COST;
                    start += MOVE_COST;

                    // Aktualizujemy status rezonatora
                    ctx.q_on_res[p_in] = chosen_res;
                    ctx.res_occ[chosen_res] = p_in;
                    ctx.q_depth[p_in]++;
                } else {
                    start = fast_max_hyb(ctx.q_time[p_in], fast_max_hyb(ctx.q_time[p_stay], ctx.r_time[chosen_res]));
                }

                // Logujemy wykonanie właściwej bramki 2-kubitowej (LOCAL - czyli w obrębie sąsiadów)
                if constexpr (LOGGING) out_log->push_back({OpType::GATE_2Q_LOCAL, node.id, ctx.p2l[p_in], ctx.p2l[p_stay], chosen_res + RES_OFFSET, p_stay, p_in, start});

                // Aktualizacja czasów gotowości zasobów o jednostkę czasu wykonania bramki (1.0)
                ctx.q_time[p_in] = ctx.q_time[p_stay] = ctx.r_time[chosen_res] = start + 1.0;

                // Synchronizacja warstw logicznych
                int glayer = fast_max_hyb(ctx.q_depth[p_in], ctx.q_depth[p_stay]);
                ctx.q_depth[p_in] = ctx.q_depth[p_stay] = glayer + 1;
            }
        }

        // Zwalniamy w DAGu zależności (dzieci tego węzła mogą teraz wejść do front_layer)
        for (int nxt : node.next_gates) {
            if (--deps[nxt] == 0) front_layer.push_back(nxt);
        }
    }
}

// ==============================================================================
// KROK ROUTINGU (de_route_step)
// Szuka najlepszego ruchu (SWAP/MOVE), gdy bramki są zablokowane.
// ==============================================================================
template <bool LOGGING>
bool SabreDEHybridRouter::de_route_step(
    const HybridSabreDAG& current_dag,
    const std::vector<int>& front_layer,
    std::vector<double>& decay,
    DERunContext& ctx,
    const double* genes, int n_gates,
    std::vector<RoutingRecord>* out_log)
{
    // 1. Zbieranie aktywnych kubitów
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

    // 2. [SUPER-OPTYMALIZACJA] Hoisting lokalizacji logicznych (Płaska Struktura)
    ctx.ws.f_locs.clear();
    for (int gid : front_layer) {
        const auto& node = current_dag.nodes[gid];
        if (node.lq2 != -1) ctx.ws.f_locs.push_back({ctx.l2p[node.lq1], ctx.l2p[node.lq2]});
    }

    // Przygotowanie płaskiego wektora lokacji z offsetami.
    // Działa tak samo jak BFS wyżej: zapisujemy wszystko w jednym ciągu, by maksymalizować
    // Cache Hit Rate i w ogóle nie alokować niczego w głębokich pętlach.
    ctx.ws.flat_lookahead_locs.clear();
    ctx.ws.lookahead_locs_offsets.clear();
    ctx.ws.lookahead_locs_offsets.push_back(0); // Punkt startowy dla lokalizacji

    // Iterujemy przez wyliczone w ewaluatorze płaskie warstwy przy użyciu offsetów
    for (size_t i = 0; i < ctx.ws.lookahead_offsets.size() - 1; ++i) {
        int start_idx = ctx.ws.lookahead_offsets[i];
        int end_idx = ctx.ws.lookahead_offsets[i+1];

        // Wyciągamy bramki należące tylko do aktualnie przetwarzanej warstwy (i)
        for (int j = start_idx; j < end_idx; ++j) {
            int gid = ctx.ws.flat_lookahead_layers[j];
            const auto& node = current_dag.nodes[gid];
            if (node.lq2 != -1) {
                ctx.ws.flat_lookahead_locs.push_back({ctx.l2p[node.lq1], ctx.l2p[node.lq2]});
            }
        }
        // Zatwierdzamy zakończenie warstwy, dodając na jej koniec nowy offset
        ctx.ws.lookahead_locs_offsets.push_back(ctx.ws.flat_lookahead_locs.size());
    }

    ctx.ws.route_step_counter++;
    int step_id = ctx.ws.route_step_counter;
    int p_cnt = topo.qubits.size();

    double best_score = std::numeric_limits<double>::max();
    int best_p1 = -1, best_p2 = -1, best_r = -1;

    // 3. Poszukiwanie najlepszej krawędzi (przeszukiwanie liniowe po Cache)
    for (int p : ctx.ws.active_phys) {
        for (int r : topo.qubits[p].connected_resonators) {
            for (int nb : topo.resonators[r].connected_qubits) {
                if (p == nb) continue;
                int u = std::min(p, nb), v = std::max(p, nb);

                int edge_idx = u * p_cnt + v;
                if (ctx.ws.seen_edges[edge_idx] == step_id) continue;
                ctx.ws.seen_edges[edge_idx] = step_id;

                double sF = 0;

                // Badanie natychmiastowych potrzeb (front_layer) z bufora
                for (const auto& loc : ctx.ws.f_locs) {
                    int l1 = loc.first, l2 = loc.second;
                    if (l1 == p) l1 = nb; else if (l1 == nb) l1 = p;
                    if (l2 == p) l2 = nb; else if (l2 == nb) l2 = p;
                    sF += get_distance(l1, l2);
                }

                double f_norm = ctx.ws.f_locs.empty() ? 0 : (sF / ctx.ws.f_locs.size());

                // Badanie przyszłości (Lookahead) z PŁASKIEGO bufora z offsetami
                double ext_score = 0.0;
                double layer_weight = EXTENDED_LAYER_WEIGHT;

                // Uwaga: lookahead_locs_offsets.size() przechowuje ilość offsetów, czyli (liczba warstw + 1)
                for (size_t i = 0; i < ctx.ws.lookahead_locs_offsets.size() - 1; ++i) {
                    int start_idx = ctx.ws.lookahead_locs_offsets[i];
                    int end_idx = ctx.ws.lookahead_locs_offsets[i+1];
                    int layer_size = end_idx - start_idx;

                    // Jeżeli warstwa nie miała żadnych bramek dwukubitowych
                    if (layer_size == 0) {
                        layer_weight *= 0.5;
                        continue;
                    }

                    double layer_score = 0.0;
                    // Iteracja odbywa się w jednym spójnym bloku pamięci Cache
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

                // Mutacja wprowadzająca czynnik losowy z algorytmu DE
                double mutation = 0.95 + (genes[u % n_gates] * 0.1);
                double score = (f_norm + ext_score) * fast_max_hyb(decay[u], decay[v]) * mutation;

                if (score < best_score) {
                    best_score = score;
                    best_p1 = u; best_p2 = v; best_r = r;
                }
            }
        }
    }

    // 4. Wykonanie ruchu na najlepszym kandydacie
    if (best_p1 != -1) {
        int p1 = best_p1, p2 = best_p2, r_step = best_r;

        de_force_down<LOGGING>(p1, -1, ctx, genes, n_gates, out_log, false);
        de_force_down<LOGGING>(p2, -1, ctx, genes, n_gates, out_log, false);
        if (ctx.res_occ[r_step] != -1) de_force_down<LOGGING>(ctx.res_occ[r_step], -1, ctx, genes, n_gates, out_log, false);

        double step_start = fast_max_hyb(ctx.q_time[p1], fast_max_hyb(ctx.r_time[r_step], ctx.q_time[p2]));
        bool is_move = (ctx.p2l[p2] == -1) || (ctx.p2l[p1] == -1);
        double cost = is_move ? 2.0 * MOVE_COST : SWAP_COST;

        if constexpr (LOGGING) {
            if (!is_move) out_log->push_back({OpType::SWAP, -1, ctx.p2l[p1], ctx.p2l[p2], p1, p2, r_step, step_start});
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

// ==============================================================================
// 1. SZYBKI RESET STANU
// Czyści i przygotowuje kontekst przed symulacją.
// ==============================================================================
inline void SabreDEHybridRouter::reset_evaluation_state(DERunContext& ctx, std::vector<RoutingRecord>* out_log) {
    // Kopiowanie początkowego mapowania logicznego na fizyczne (szybkie memmove)
    std::copy(saved_initial_l2p.begin(), saved_initial_l2p.end(), ctx.l2p.begin());
    std::fill(ctx.p2l.begin(), ctx.p2l.end(), -1);

    // Odwrotne mapowanie
    for (int i = 0; i < num_qubits; ++i) {
        if (ctx.l2p[i] != -1) ctx.p2l[ctx.l2p[i]] = i;
    }

    // Resetowanie tablic pomocniczych stanu kubitów i rezonatorów
    std::fill(ctx.q_on_res.begin(), ctx.q_on_res.end(), -1);
    std::fill(ctx.res_occ.begin(), ctx.res_occ.end(), -1);
    std::fill(ctx.q_time.begin(), ctx.q_time.end(), 0.0);
    std::fill(ctx.r_time.begin(), ctx.r_time.end(), 0.0);
    std::fill(ctx.q_depth.begin(), ctx.q_depth.end(), 0);

    if (out_log) { out_log->clear(); }
}

// ==============================================================================
// 2. PRZYGOTOWANIE WARSTWY FRONTOWEJ
// Znajduje bramki początkowe (stopień wejściowy = 0) na start symulacji.
// ==============================================================================
inline void SabreDEHybridRouter::init_front_layer(int n_gates) {
    workspace.front_layer.clear();
    std::copy(dag.initial_deps.begin(), dag.initial_deps.end(), workspace.deps.begin());

    // Zbieramy bramki bez zależności wjazdowych
    for (int i = 0; i < n_gates; ++i) {
        if (workspace.deps[i] == 0) workspace.front_layer.push_back(i);
    }

    std::fill(workspace.decay.begin(), workspace.decay.end(), 1.0);

    // Inkrementacja epoki BFS na starcie (O(1) reset stanu w in_ext)
    workspace.current_bfs_epoch++;
}

// ==============================================================================
// 3. WYCIĄGANIE GOTOWYCH BRAMEK (TWO-POINTERS OPTIMIZATION)
// Filtruje front_layer przenosząc gotowe bramki do exec_gates w miejscu.
// ==============================================================================
inline void SabreDEHybridRouter::extract_executable_gates(DERunContext& ctx) {
    workspace.exec_gates.clear();
    size_t write_idx = 0;

    for (size_t read_idx = 0; read_idx < workspace.front_layer.size(); ++read_idx) {
        int g_id = workspace.front_layer[read_idx];
        const auto& node = dag.nodes[g_id];

        // Sprawdzenie możliwości wykonania: operacja 1-kubitowa lq2 == -1
        // LUB sąsiedztwo dla operacji 2-kubitowych
        bool exec = (node.lq2 == -1);
        if (!exec) {
            int p1 = ctx.l2p[node.lq1];
            int p2 = ctx.l2p[node.lq2];
            exec = (p1 != -1 && p2 != -1 && get_distance(p1, p2) <= 1);
        }

        // Jeśli gotowa, wrzuć do kolejki wykonania. Jeśli nie, zachowaj we front layer
        if (exec) {
            workspace.exec_gates.push_back(g_id);
        } else {
            workspace.front_layer[write_idx++] = g_id;
        }
    }
    // Ucięcie resztek z wektora (O(1))
    workspace.front_layer.resize(write_idx);
}

// ==============================================================================
// 4. BUDOWANIE LOOKAHEAD (PŁASKIE BFS Z EPOKAMI)
// Poszukuje dalszych zależności bez alokacji nowych wektorów dla każdej głębokości.
// ==============================================================================
inline void SabreDEHybridRouter::build_lookahead_layers() {
    workspace.flat_lookahead_layers.clear();
    workspace.lookahead_offsets.clear();
    workspace.lookahead_offsets.push_back(0);

    // O(1) resetowanie odwiedzin
    workspace.current_bfs_epoch++;

    for (int gid : workspace.front_layer) {
        workspace.visited_epoch[gid] = workspace.current_bfs_epoch;
    }

    workspace.exec_gates = workspace.front_layer;

    // Iteracyjna budowa warstw
    for (size_t depth = 0; depth < max_look_ahead_depth && !workspace.exec_gates.empty(); ++depth) {
        workspace.ext_layer.clear();

        for (int gid : workspace.exec_gates) {
            for (int nxt : dag.nodes[gid].next_gates) {
                // Dodajemy bramkę tylko, jeśli nie była ruszana w obecnej epoce
                if (workspace.visited_epoch[nxt] != workspace.current_bfs_epoch) {
                    workspace.visited_epoch[nxt] = workspace.current_bfs_epoch;
                    workspace.ext_layer.push_back(nxt);
                    workspace.flat_lookahead_layers.push_back(nxt);
                }
            }
        }

        // Rejestracja offsetu płaskiego widoku
        if (!workspace.ext_layer.empty()) {
            workspace.lookahead_offsets.push_back(workspace.flat_lookahead_layers.size());
        }

        // Swap to najszybszy sposób na podmianę wektorów bez kopiowania
        std::swap(workspace.exec_gates, workspace.ext_layer);
    }

    workspace.exec_gates.clear();
    workspace.ext_layer.clear();
}

// ==============================================================================
// 5. FAZA KOŃCOWA Z OBLICZENIEM WYNIKU
// Wymusza sprowadzenie stanów, sortuje logi i oblicza funkcję celu (fitness).
// ==============================================================================
template <bool LOGGING>
double SabreDEHybridRouter::finalize_evaluation(DERunContext& ctx, const double* genes, int n_gates, int p_cnt,
                                                int& out_swaps, int& out_moves, double& max_physical_time,
                                                int& max_logical_depth, std::vector<RoutingRecord>* out_log) {
    // 1. Force down pozostałych stanów
    for (size_t p = 0; p < p_cnt; ++p) {
        int temp_p = p;
        de_force_down<LOGGING>(temp_p, -1, ctx, genes, n_gates, out_log, false);
    }

    // 2. Sortowanie logowania (tylko w trybie LOGGING)
    if constexpr (LOGGING) {
        if (out_log && !out_log->empty()) {
            std::stable_sort(out_log->begin(), out_log->end(), [](const RoutingRecord& a, const RoutingRecord& b) {
                return a.start_time < b.start_time;
            });
        }
    }

    // 3. Zebranie metryk max
    for (double t : ctx.q_time) max_physical_time = fast_max_hyb(max_physical_time, t);
    for (int d : ctx.q_depth) max_logical_depth = fast_max_hyb(max_logical_depth, d);

    out_swaps = ctx.out_swaps;
    out_moves = ctx.out_moves;

    // 4. Obliczenie wartości fitness
    return (ctx.total_routing_cost * 1000.0) + (out_swaps * 5000.0) + (max_physical_time * 1.0) + (max_logical_depth * 0.5);
}

// ==============================================================================
// GŁÓWNA FUNKCJA EWALUACJI
// Używa nowo wydzielonych metod pomocniczych, pozostając ultra czytelną.
// ==============================================================================
template <bool LOGGING>
double SabreDEHybridRouter::evaluate_fitness_impl(const double* genes, int& out_swaps, int& out_moves, double& max_physical_time, int& max_logical_depth, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard) {
    int n_gates = dag.nodes.size();
    int p_cnt = topo.qubits.size();

    out_swaps = 0;
    out_moves = 0;
    max_physical_time = 0.0;
    max_logical_depth = 0;

    DERunContext ctx(workspace);

    // 1. Reset i ustawienie warstwy początkowej
    reset_evaluation_state(ctx, out_log);
    init_front_layer(n_gates);

    int gates_routed = 0;

    // 2. Główna pętla symulacyjna
    while (gates_routed < n_gates) {

        // Oddzielenie bramek możliwych do wykonania (exec_gates) od front_layer
        extract_executable_gates(ctx);

        // Wykonaj bramki, jeśli cokolwiek znajduje się w exec_gates
        if (!workspace.exec_gates.empty()) {
            de_execute_gates<LOGGING>(dag, workspace.exec_gates, workspace.front_layer, workspace.deps, ctx, genes, n_gates, out_log, can_lazy_offboard);
            gates_routed += workspace.exec_gates.size();
            std::fill(workspace.decay.begin(), workspace.decay.end(), 1.0);
            continue;
        }

        // 3. Budowanie spłaszczonego Lookahead BFS
        build_lookahead_layers();

        // Krok decyzyjny DE (SWAPy/przesunięcia)
        if (!de_route_step<LOGGING>(dag, workspace.front_layer, workspace.decay, ctx, genes, n_gates, out_log)) {
            break; // Router zablokował się
        }
    }

    // 4. Pobranie i zwrócenie finalnych statystyk fitness
    return finalize_evaluation<LOGGING>(ctx, genes, n_gates, p_cnt, out_swaps, out_moves, max_physical_time, max_logical_depth, out_log);
}

// ==============================================================================
// WYZWALACZ EWALUACJI (evaluate_fitness)
// Funkcja fasadowa decydująca w czasie kompilacji o włączeniu logowania.
// ==============================================================================
double SabreDEHybridRouter::evaluate_fitness(const double* genes, int& out_swaps, int& out_moves, double& max_physical_time, int& max_logical_depth, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard) {
    // Template dispatching: Unikamy instrukcji warunkowych w gorącej pętli symulacji.
    // Jeśli wywołujemy tę funkcję prosząc o out_log (tylko na koniec dla zwycięzcy), odpala wersję <true>.
    // W trakcie milionów iteracji DE odpala wersję <false>, oszczędzając zasoby.
    if (out_log) return evaluate_fitness_impl<true>(genes, out_swaps, out_moves, max_physical_time, max_logical_depth, out_log, can_lazy_offboard);
    return evaluate_fitness_impl<false>(genes, out_swaps, out_moves, max_physical_time, max_logical_depth, nullptr, can_lazy_offboard);
}

// ==============================================================================
// Inicjalizowanie populacji początkowej z "ziarnem" od SABRE
// ==============================================================================
void SabreDEHybridRouter::initialize_de_population(int pop_size, int n_genes, int n_gates, int blocks, const std::vector<int>& logical_exec_order, std::vector<double>& pop, std::vector<double>& fitness, boost::random::xoshiro256pp& rng, bool can_lazy) {
    boost::random::uniform_real_distribution<double> dist_real(0.0, 1.0);

    // Tworzymy genotyp dla osobnika bazowego (ziarna SABRE). Domyślnie wartości 0.5.
    std::vector<double> sabre_genes(n_genes, 0.5);

    // Na podstawie kolejności wyznaczonej przez SABRE nadajemy genom odpowiednie priorytety.
    // Im wcześniej bramka została wykonana w SABRE, tym wyższą wartość genu (priorytet) otrzymuje.
    for (size_t i = 0; i < logical_exec_order.size(); ++i) {
        int gid = logical_exec_order[i];
        sabre_genes[gid] = 1.0 - (static_cast<double>(i) / logical_exec_order.size());

        // Pozostałe bloki genów (np. decydujące o wyborze rezonatora) są w ziarnie losowane.
        for (int b = 1; b < blocks; ++b) sabre_genes[b * n_gates + gid] = dist_real(rng);
    }

    int d_s, d_m, d_d; double d_t;
    std::vector<double> trial(n_genes);

    // Wypełnianie populacji
    for (int i = 0; i < pop_size; ++i) {
        for (int j = 0; j < n_genes; ++j) {
            // Osobnik 0 dostaje wyliczone ziarno SABRE, reszta jest kompletnie losowa.
            trial[j] = pop[i * n_genes + j] = (i == 0) ? sabre_genes[j] : dist_real(rng);
        }

        // Symulacja obwodu na podstawie genów i obliczenie funkcji przystosowania (fitness).
        fitness[i] = evaluate_fitness(trial.data(), d_s, d_m, d_t, d_d, nullptr, can_lazy);

        // Dodatkowe wypisanie logów dla pierwszego osobnika (ziarna).
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

    double final_score = evaluate_fitness(pop.data() + best_idx * n_genes, final_swaps, final_moves, final_makespan, final_depth, &final_log, can_lazy);

    std::cout << "\n=== WALIDACJA OBWODU LOGICZNEGO ===\n";
    if (!validate_logical_circuit(final_log)) exit(-1);
    std::cout << "===================================\n";

    // print_optimized_circuit(final_log);
    print_statistics(final_score, final_makespan, final_depth, final_swaps, final_moves, std::chrono::duration<double>(end_time - start_time).count());
    export_moves_to_file(final_log, "../hybrid_schedule_moves_" + name + ".txt");
}

// ==============================================================================
// Przygotowanie środowiska, wygenerowanie "ziarna" za pomocą SABRE, a następnie zorganizowanie pętli przesuwającej okno wzdłuż osi czasu obwodu logicznego.
// ==============================================================================
void SabreDEHybridRouter::optimize_windowed(int pop_size, int generations_per_window, int window_size, double accept_ratio, bool allow_lazy_offboard, double f, double cr) {
    std::cout << "\n[Hybryda 3.0 (WINDOWED)] Faza 1: Generowanie ziarna (Priority Seed) przez heurystykę SABRE...\n";

    // 1. Generujemy bazowy harmonogram logiczny przy pomocy czystej heurystyki SABRE.
    std::vector<int> logical_exec_order = prepare_sabre_seed();
    int n_gates = dag.nodes.size();
    if (n_gates == 0) return; // Zabezpieczenie przed pustym obwodem

    // 2. Sprawdzenie, czy możemy użyć "Lazy Off-boarding" (odkładanie kubitów na wolne miejsca).
    bool can_lazy = allow_lazy_offboard && (num_qubits < topo.qubits.size());
    if (can_lazy) std::cout << "Aktywowano Lazy Off-boarding (wolne kubity fizyczne dostępne).\n";

    // Wyliczamy liczbę bloków genów. Jeśli mamy wolne kubity, dodajemy dodatkowy blok dla logiki lazy.
    int blocks = can_lazy ? 4 : 3;
    int n_genes = n_gates * blocks;

    std::cout << "[Hybryda 3.0] Faza 2: DE Windowed na bazie logicznej (Rozmiar okna: " << window_size << ")...\n";
    auto start_time = std::chrono::high_resolution_clock::now();

    // 3. Inicjalizacja populacji i generatora liczb pseudolosowych
    std::vector<double> pop(pop_size * n_genes), fitness(pop_size);
    boost::random::xoshiro256pp rng(1337); // Użycie stałego ziarna (1337) zapewnia powtarzalność wyników

    // Populacja jest inicjalizowana m.in. ziarnem z SABRE, reszta jest losowa.
    initialize_de_population(pop_size, n_genes, n_gates, blocks, logical_exec_order, pop, fitness, rng, can_lazy);

    // Krok przesuwania okna zależy od accept_ratio (np. okno 100, ratio 0.5 -> przesuwamy o 50).
    int step = std::max(1, static_cast<int>(window_size * accept_ratio));
    int d_s, d_m, d_d; double d_t; // Zmienne pomocnicze dla statystyk ewaluacji

    // 4. GŁÓWNA PĘTLA OKIENKOWA
    for (int w_start = 0; w_start < n_gates; w_start += step) {
        // Zabezpieczenie, aby okno nie wyszło poza liczbę bramek
        int w_end = std::min(w_start + window_size, n_gates);
        std::cout << ">>> Przetwarzanie okna: operacje [" << w_start << " - " << w_end - 1 << "] / " << n_gates << "\n";

        // 5. Wywołanie Ewolucji Różnicowej TYLKO dla aktualnego okna
        run_evolution_loop_windowed(pop, fitness, pop_size, generations_per_window, n_genes, n_gates, f, cr, rng, w_start, w_end);

        // Znalezienie najlepszego osobnika w populacji po ewolucji tego okna
        auto best_iter = std::min_element(fitness.begin(), fitness.end());
        int best_idx = std::distance(fitness.begin(), best_iter);

        // Określenie, do którego miejsca "akceptujemy" zoptymalizowane rozwiązanie
        int accept_end = std::min(w_start + step, n_gates);

        // 6. ROZPROPAGOWANIE SUKCESU
        // Kopiujemy najlepsze geny z zaakceptowanej części okna do całej reszty populacji.
        // Dzięki temu w następnym kroku (dla kolejnego okna) cała populacja zaczyna
        // z wypracowaną do tej pory optymalną "przeszłością".
        for (int i = 0; i < pop_size; ++i) {
            if (i == best_idx) continue; // Najlepszy osobnik już to ma
            for (int blk = 0; blk < blocks; ++blk) {
                for (int g = w_start; g < accept_end; ++g) {
                    pop[i * n_genes + blk * n_gates + g] = pop[best_idx * n_genes + blk * n_gates + g];
                }
            }
            // Aktualizacja fitness po nadpisaniu genów z przeszłości
            fitness[i] = evaluate_fitness(pop.data() + i * n_genes, d_s, d_m, d_t, d_d, nullptr, can_lazy);
        }
    }

    // 7. Podsumowanie, walidacja i eksport wyników
    auto end_time = std::chrono::high_resolution_clock::now();
    auto best_it = std::min_element(fitness.begin(), fitness.end());
    int best_idx = std::distance(fitness.begin(), best_it);

    std::vector<RoutingRecord> final_log;
    int final_swaps = 0, final_moves = 0, final_depth = 0; double final_makespan = 0.0;
    double final_score = evaluate_fitness(pop.data() + best_idx * n_genes, final_swaps, final_moves, final_makespan, final_depth, &final_log, can_lazy);

    std::cout << "\n=== WALIDACJA OBWODU LOGICZNEGO ===\n";
    if (!validate_logical_circuit(final_log)) exit(-1);
    std::cout << "===================================\n";

    // print_optimized_circuit(final_log);
    print_statistics(final_score, final_makespan, final_depth, final_swaps, final_moves, std::chrono::duration<double>(end_time - start_time).count());
    export_moves_to_file(final_log, "../hybrid_schedule_moves_windowed_" + name + ".txt");
}

// ==============================================================================
// Wykonywanie pokoleń, ale z ograniczeniem zmian tylko do bieżącego okna. Wybiera kandydatów do mutacji i ocenia, czy nowe rozwiązanie jest lepsze od starego.
// ==============================================================================
void SabreDEHybridRouter::run_evolution_loop_windowed(std::vector<double>& pop, std::vector<double>& fitness, int pop_size, int generations, int n_genes, int n_gates, double f, double cr, boost::random::xoshiro256pp& rng, int w_start, int w_end) {
    boost::random::uniform_int_distribution<int> dist_int(0, pop_size - 1);
    std::vector<double> trial(n_genes);
    int d_s, d_m, d_d; double d_t;
    bool can_lazy = (n_genes / n_gates == 4); // Sprytne wyciągnięcie flagi z proporcji genów

    for (int gen = 0; gen < generations; ++gen) {
        for (int i = 0; i < pop_size; ++i) {
            // 1. Wybór 3 RÓŻNYCH losowych osobników (a, b, c) z populacji,
            // żaden z nich nie może być też obecnie ewaluowanym osobnikiem (i).
            int a, b, c;
            do { a = dist_int(rng); } while (a == i);
            do { b = dist_int(rng); } while (b == i || b == a);
            do { c = dist_int(rng); } while (c == i || c == a || c == b);

            // 2. Wygenerowanie nowego osobnika (wektora próbnego) - zmiany tylko w oknie!
            generate_trial_vector_windowed(pop, trial, i, a, b, c, n_genes, n_gates, f, cr, rng, w_start, w_end);

            // 3. Ocena funkcji przystosowania (fitness) dla nowego wektora próbnego
            double trial_fit = evaluate_fitness(trial.data(), d_s, d_m, d_t, d_d, nullptr, can_lazy);

            // 4. Selekcja: Jeśli nowy osobnik jest lepszy (ma mniejszy koszt),
            // zastępuje starego osobnika w bieżącej generacji populacji.
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
void SabreDEHybridRouter::generate_trial_vector_windowed(const std::vector<double>& pop, std::vector<double>& trial, int current_idx, int a, int b, int c, int n_genes, int n_gates, double f, double cr, boost::random::xoshiro256pp& rng, int w_start, int w_end) {
    boost::random::uniform_real_distribution<double> dist_real(0.0, 1.0);

    // Gwarancja mutacji przynajmniej jednego genu (wewnątrz obsługiwanego okna)
    // Aby to osiągnąć, losowany jest indeks R_idx w zmapowanej przestrzeni tylko dla okna,
    // a następnie jest przeliczany z powrotem na globalny indeks genu (R).
    int R_idx = boost::random::uniform_int_distribution<int>(0, (w_end - w_start) * (n_genes / n_gates) - 1)(rng);
    int R = (R_idx / (w_end - w_start)) * n_gates + w_start + (R_idx % (w_end - w_start));

    // Wskaźniki na początek bloków pamięci odpowiednich osobników
    const double *p_a = pop.data() + a * n_genes;
    const double *p_b = pop.data() + b * n_genes;
    const double *p_c = pop.data() + c * n_genes;
    const double *p_curr = pop.data() + current_idx * n_genes;

    for (int j = 0; j < n_genes; ++j) {
        int gate_id = j % n_gates; // Określenie, do której bramki należy ten konkretny gen

        // WARUNEK OKNA: Czy bramka mieści się w aktywnym oknie?
        bool is_in_window = (gate_id >= w_start && gate_id < w_end);

        // WARUNEK CROSSOVERU: Czy wylosowaliśmy szansę na mutację lub jesteśmy wymaganym genem R?
        bool do_mutate = (dist_real(rng) < cr || j == R);

        if (is_in_window && do_mutate) {
            // Ewolucja Różnicowa (DE/rand/1):
            // Mutacja wektora bazowego (a) różnicą dwóch pozostałych (b - c) przeskalowaną czynnikiem F
            double mutated_value = p_a[j] + f * (p_b[j] - p_c[j]);

            // Obcięcie (clamping) do dopuszczalnych granic <0.0, 1.0>
            trial[j] = std::clamp(mutated_value, 0.0, 1.0);
        } else {
            // Skopiowanie starego genu (zablokowanego przez bycie poza oknem lub brak mutacji)
            trial[j] = p_curr[j];
        }
    }
}

// ==============================================================================
// sabre_force_qubit_down
// Odpowiada za ściągnięcie fizycznego kubitu z rezonatora (tzw. off-boarding).
// ==============================================================================
void SabreDEHybridRouter::sabre_force_qubit_down(int p, int current_gate_id, SabreRunContext& ctx) {
    // Sprawdzamy, czy dany kubit fizyczny 'p' znajduje się obecnie na jakimkolwiek rezonatorze
    if (ctx.state.q_on_resonator[p] != -1) {
        int r = ctx.state.q_on_resonator[p];

        // Czas startu tej operacji: musimy poczekać, aż ZARÓWNO kubit, jak i rezonator
        // zakończą swoje dotychczasowe operacje. Wybieramy późniejszy czas.
        double step_start = fast_max_hyb(ctx.state.q_ready_time[p], ctx.state.r_ready_time[r]);

        // Rejestrujemy statystyki (off-boarding to zawsze jeden ruch/MOVE)
        ctx.out_moves++;
        ctx.total_routing_cost += MOVE_COST;

        // Aktualizujemy czas gotowości - kubit i rezonator będą zajęte przez czas trwania ruchu
        ctx.state.q_ready_time[p] = step_start + MOVE_COST;
        ctx.state.r_ready_time[r] = step_start + MOVE_COST;

        // Skoro kubit wykonał operację fizyczną, podbijamy jego głębokość w harmonogramie
        ctx.q_logical_depth[p] += 1;

        // Zwalniamy powiązania - kubit już nie jest na rezonatorze, a rezonator jest pusty
        ctx.state.q_on_resonator[p] = -1;
        ctx.state.res_occupant[r] = -1;
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
    // 1. Zbieramy "aktywne kubity" - tylko te fizyczne, które są zaangażowane
    // w bramki aktualnie czekające we front_layer, opłaca nam się przesuwać.
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

    // 2. Szukamy krawędzi (kandydatów na wymianę).
    // Dla każdego aktywnego kubitu sprawdzamy sąsiadów dostępnych poprzez jego rezonatory.
    for (int p : active_phys_qubits) {
        for (int r : topo.qubits[p].connected_resonators) {
            for (int neighbor_p : topo.resonators[r].connected_qubits) {
                if (p == neighbor_p) continue;

                // Standaryzacja identyfikatora krawędzi (mniejsze id najpierw), by uniknąć duplikatów
                int u = std::min(p, neighbor_p);
                int v = std::max(p, neighbor_p);

                // O(N) deduplikacja (można to przyspieszyć hashsetem lub flagami)
                bool already_added = false;
                for (const auto& c : candidates) {
                    if (c.p1 == u && c.p2 == v) { already_added = true; break; }
                }
                if (already_added) continue;

                // --- 3. Obliczanie punktacji (SCORE) dla ruchu ---
                double score_F = 0.0;

                // Badamy wpływ na front_layer (natychmiastowe potrzeby)
                for (int gate_id : front_layer) {
                    const auto& n = current_dag.nodes[gate_id];
                    if (n.lq2 != -1) {
                        int lq1_p = ctx.state.logical_to_phys_qubit[n.lq1];
                        int lq2_p = ctx.state.logical_to_phys_qubit[n.lq2];

                        // Symulacja zamiany
                        if (lq1_p == p) lq1_p = neighbor_p; else if (lq1_p == neighbor_p) lq1_p = p;
                        if (lq2_p == p) lq2_p = neighbor_p; else if (lq2_p == neighbor_p) lq2_p = p;

                        score_F += get_distance(lq1_p, lq2_p);
                    }
                }

                double f_norm = front_layer.empty() ? 0 : (score_F / front_layer.size());

                // --- Wpływ na warstwy w głębi (Lookahead) z zanikającą wagą ---
                double ext_score = 0.0;
                double layer_weight = EXTENDED_LAYER_WEIGHT; // Bazowo np. 0.5

                for (const auto& layer : lookahead_layers) {
                    if (layer.empty()) continue;

                    double layer_score = 0.0;
                    for (int gate_id : layer) {
                        const auto& n = current_dag.nodes[gate_id];
                        if (n.lq2 != -1) {
                            int lq1_p = ctx.state.logical_to_phys_qubit[n.lq1];
                            int lq2_p = ctx.state.logical_to_phys_qubit[n.lq2];

                            // Symulacja zamiany miejscami
                            if (lq1_p == p) lq1_p = neighbor_p; else if (lq1_p == neighbor_p) lq1_p = p;
                            if (lq2_p == p) lq2_p = neighbor_p; else if (lq2_p == neighbor_p) lq2_p = p;

                            layer_score += get_distance(lq1_p, lq2_p);
                        }
                    }
                    // Dodajemy znormalizowany wynik z odpowiednią wagą
                    ext_score += layer_weight * (layer_score / layer.size());

                    // Każda kolejna warstwa w głąb jest o połowę mniej ważna
                    layer_weight *= 0.5;
                }

                // OSTATECZNY WYNIK: łączymy obie metryki i uwzględniamy karę za ping-pong (decay)
                double score = (f_norm + ext_score) * fast_max_hyb(ctx.decay[u], ctx.decay[v]);
                candidates.push_back({u, v, r, score});
            }
        }
    }

    // 4. Wybór najlepszego ruchu (najniższy score)
    SwapCandidate* best_candidate = nullptr;
    if (!candidates.empty()) {
        best_candidate = &candidates[0];
        for (size_t i = 1; i < candidates.size(); ++i) {
            if (candidates[i].score < best_candidate->score) best_candidate = &candidates[i];
        }
    }

    // 5. Egzekucja wybranego ruchu
    if (best_candidate != nullptr) {
        int p1 = best_candidate->p1, p2 = best_candidate->p2, r_step = best_candidate->r_shared;

        // Zanim je zamienimy, zrzucamy oba z rezonatorów do stanu spoczynku
        sabre_force_qubit_down(p1, -1, ctx);
        sabre_force_qubit_down(p2, -1, ctx);

        // Zwalniamy sam rezonator "pośredniczący", przez który przeprowadzimy ruch
        if (ctx.state.res_occupant[r_step] != -1) {
            sabre_force_qubit_down(ctx.state.res_occupant[r_step], -1, ctx);
        }

        // Kiedy wszystko będzie gotowe na realizację ruchu?
        double step_start = fast_max_hyb(ctx.state.q_ready_time[p1], fast_max_hyb(ctx.state.r_ready_time[r_step], ctx.state.q_ready_time[p2]));

        // Weryfikacja: jeśli jedno z miejsc to wolny ("pusty") kubit fizyczny (Lazy Offboard), to nie jest to SWAP, tylko MOVE.
        bool is_move = (ctx.state.phys_qubit_to_logical[p2] == -1) || (ctx.state.phys_qubit_to_logical[p1] == -1);
        double op_cost = is_move ? (2.0 * MOVE_COST) : SWAP_COST;

        int lq_1 = ctx.state.phys_qubit_to_logical[p1], lq_2 = ctx.state.phys_qubit_to_logical[p2];

        // Zliczamy statystyki kosztów do raportowania
        if (is_move) ctx.out_moves += 2; else ctx.out_swaps++;
        ctx.total_routing_cost += op_cost;

        // FAKTYCZNA ZAMIANA MAPOWANIA W PAMIĘCI
        std::swap(ctx.state.phys_qubit_to_logical[p1], ctx.state.phys_qubit_to_logical[p2]);
        if (lq_1 != -1) ctx.state.logical_to_phys_qubit[lq_1] = p2;
        if (lq_2 != -1) ctx.state.logical_to_phys_qubit[lq_2] = p1;

        // Aktualizacja czasów warstw (MOVE zużywa zazwyczaj 2 kroki czasowe/zależności, SWAP 1 duży)
        int current_layer = fast_max_hyb(ctx.q_logical_depth[p1], ctx.q_logical_depth[p2]);
        ctx.q_logical_depth[p1] = current_layer + (is_move ? 2 : 1);
        ctx.q_logical_depth[p2] = current_layer + (is_move ? 2 : 1);

        // Zamrażamy czasowo oba kubity i rezonator na czas przeprowadzanej przeprowadzki
        ctx.state.q_ready_time[p1] = step_start + op_cost;
        ctx.state.r_ready_time[r_step] = step_start + op_cost;
        ctx.state.q_ready_time[p2] = step_start + op_cost;

        // PODBICIE KARY DECAY: zniechęca algorytm, by w następnym kroku znowu nie tknął tych samych kubitów
        ctx.decay[p1] += 0.2; ctx.decay[p2] += 0.2;

        return true; // Sukces: zrobiliśmy krok naprzód
    }
    return false; // Błąd: algorytm nie miał żadnego możliwego ruchu (w praktyce rzadkie)
}

// ==============================================================================
// GŁÓWNY PRZEBIEG HEURYSTYKI SABRE (DO PRZODU LUB WSTECZ)
// ==============================================================================
std::vector<int> SabreDEHybridRouter::run_sabre_pass(const HybridSabreDAG& current_dag, const std::vector<int>& initial_l2p, bool is_final, std::vector<int>* out_logical_order) {
    int n_gates = current_dag.nodes.size();

    // --- 3A. Przygotowanie czystego stanu ---
    SabreRunContext ctx;
    ctx.state.logical_to_phys_qubit = initial_l2p; // Od tego mapowania zaczynamy
    ctx.state.phys_qubit_to_logical.assign(topo.qubits.size(), -1);

    // Inicjalizujemy mapowanie odwrotne (fizyczny -> logiczny)
    for (int i = 0; i < num_qubits; ++i) {
        if (initial_l2p[i] != -1) ctx.state.phys_qubit_to_logical[initial_l2p[i]] = i;
    }

    // Resetujemy stany rezonatorów i czasy ułożenia kubitów
    ctx.state.q_on_resonator.assign(topo.qubits.size(), -1);
    ctx.state.res_occupant.assign(topo.resonators.size(), -1);
    ctx.state.q_ready_time.assign(topo.qubits.size(), 0.0);
    ctx.state.r_ready_time.assign(topo.resonators.size(), 0.0);
    ctx.q_logical_depth.assign(topo.qubits.size(), 0);
    ctx.decay.assign(topo.qubits.size(), 1.0);

    // --- 3B. Inicjalizacja Warstwy Frontowej (Bramki bez zależności) ---
    std::vector<int> front_layer;
    front_layer.reserve(num_qubits);
    std::vector<int> deps = current_dag.initial_deps; // Kopia zależności do odliczania

    // Wszystkie węzły, które nie mają bramki "rodzica" wchodzą od razu na start
    for (int i = 0; i < n_gates; ++i) {
        if (deps[i] == 0) front_layer.push_back(i);
    }

    int gates_routed = 0;
    std::vector<int> executable_gates, extended_layer;
    executable_gates.reserve(num_qubits);
    extended_layer.reserve(num_qubits * 4);
    std::vector<bool> in_extended(n_gates, false);

    // Zabezpieczenie przed nieskończoną pętlą, jakby coś w grafie się zapętliło
    int debug_iterations = 0;
    const int MAX_ITERATIONS = n_gates * 100;

    // --- 3C. Główna pętla routingu ---
    while (gates_routed < n_gates) {
        if (++debug_iterations > MAX_ITERATIONS) break;

        // KROK 1: Przeszukaj front_layer i znajdź to, co możemy wykonać natychmiast
        executable_gates.clear();
        for (size_t i = 0; i < front_layer.size(); ) {
            int gate_id = front_layer[i];
            const auto& node = current_dag.nodes[gate_id];
            bool is_executable = false;

            // Bramki 1Q zawieszone same ze sobą (lq2 == -1) zawsze można wykonać
            if (node.lq2 == -1) {
                is_executable = true;
            } else {
                // Dla 2Q sprawdzamy, czy ich fizyczne reprezentacje leżą na styk (odległość <= 1)
                int p1 = ctx.state.logical_to_phys_qubit[node.lq1];
                int p2 = ctx.state.logical_to_phys_qubit[node.lq2];
                if (get_distance(p1, p2) <= 1) is_executable = true;
            }

            // Jeśli można wykonać, wyciągamy ją z front_layer i wrzucamy do executable_gates
            if (is_executable) {
                executable_gates.push_back(gate_id);
                // Trikiem "pop_back" usuwamy element z wektora w czasie O(1) podmieniając go z ostatnim
                front_layer[i] = front_layer.back();
                front_layer.pop_back();
            } else {
                ++i; // Bramka wciąż czeka, bo kubity są za daleko od siebie
            }
        }

        // KROK 2: Jeśli znaleźliśmy coś do wykonania, "odpalamy" je
        if (!executable_gates.empty()) {
            // sabre_execute_gates zaktualizuje stany, odliczy zależności i ewentualnie
            // "wpuści" nowe bramki do front_layer (jeśli ich deps spadnie do 0)
            sabre_execute_gates(current_dag, executable_gates, front_layer, deps, ctx, out_logical_order);
            gates_routed += executable_gates.size();

            // Resetujemy wagi decay (zapobieganie ping-pongom w SWAPach), bo pchnęliśmy postęp logiczny
            std::fill(ctx.decay.begin(), ctx.decay.end(), 1.0);
            continue; // Wracamy na początek while'a, by znów szukać wykonywalnych bramek
        }

        // KROK 3: Jeśli tu dotarliśmy, oznacza to, że front_layer NIE JEST PUSTY,
        // ale ŻADNA bramka nie może być wykonana. Kubity są za daleko. Musimy zrobić SWAP/MOVE.

        // Zbieramy "warstwę rozszerzoną" - czyli patrzymy krok głębiej w DAG (na bezpośrednich
        // następców front_layer), żeby widzieć w którą stronę "wypada" nam ciągnąć kubity na przyszłość.
        // KROK 3: Zbierz bramki w wielowarstwowym "Lookaheadzie" (np. do 3 warstw w głąb)
        const int LOOKAHEAD_DEPTH = 20;
        std::vector<std::vector<int>> lookahead_layers(LOOKAHEAD_DEPTH);
        std::fill(in_extended.begin(), in_extended.end(), false);

        // Zaczynamy poszukiwania od bramek we front_layer
        std::vector<int> current_level = front_layer;

        for (int d = 0; d < LOOKAHEAD_DEPTH; ++d) {
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

            // Jeśli nie ma głębszych zależności, przerywamy wchodzenie głębiej
            if (current_level.empty()) break;
        }

        // Wywołujemy heurstykę decyzyjną z wielowarstwowym kontekstem
        sabre_route_step(current_dag, front_layer, lookahead_layers, ctx);
    } // koniec while

    // --- 3D. Logowanie (tylko na ostatnim przebiegu generującym ziarno) ---
    if (is_final) {
        // Sprzątanie na koniec: upewniamy się, że wszystkie kubity zeszły z rezonatorów
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

    // Zwracamy układ kubitów (mapowanie logiczne -> fizyczne) po wykonaniu całego obwodu
    return ctx.state.logical_to_phys_qubit;
}

// ==============================================================================
// "Look-ahead/Look-back" pass
// ==============================================================================
std::vector<int> SabreDEHybridRouter::prepare_sabre_seed() {
    std::vector<int> initial_mapping(num_qubits);
    saved_initial_l2p.assign(topo.qubits.size(), -1);
    // Inicjalizujemy mapowanie tożsamościowe (kubit logiczny i = fizyczny i)
    for (int i = 0; i < num_qubits; ++i) { initial_mapping[i] = i; saved_initial_l2p[i] = i; }

    // 1. Przebieg do przodu (Forward): Ustala wstępne rozproszenie kubitów
    auto post_fwd = run_sabre_pass(dag, initial_mapping, false, nullptr);

    // 2. Przebieg wsteczny (Backward): Wykonuje sabre na odwróconym grafie (reversed_dag).
    // To kluczowy trick – pozwala "skorygować" mapowanie z końca obwodu do początku.
    auto pre_final = run_sabre_pass(reversed_dag, post_fwd, false, nullptr);

    // 3. Przebieg finalny: Generuje faktyczną kolejność wykonania bramek (logical_exec_order)
    // bazując na wypracowanym przez dwa poprzednie kroki mapowaniu.
    std::vector<int> logical_exec_order;
    run_sabre_pass(dag, pre_final, true, &logical_exec_order);

    return logical_exec_order;
}

bool SabreDEHybridRouter::validate_logical_circuit(const std::vector<RoutingRecord>& log) const {
    // 1. Budujemy oczekiwaną sekwencję ID węzłów z DAG-a dla każdego kubitu logicznego
    std::vector<std::vector<int>> expected(num_qubits);
    for (const auto& node : dag.nodes) {
        if (node.lq1 != -1) expected[node.lq1].push_back(node.id);
        if (node.lq2 != -1) expected[node.lq2].push_back(node.id);
    }

    // 2. Budujemy faktyczną sekwencję na podstawie logu operacji
    std::vector<std::vector<int>> actual(num_qubits);
    for (const auto& op : log) {
        // Interesują nas tylko faktyczne bramki logiczne, ignorujemy SWAP/MOVE
        if (op.op_type == OpType::GATE_1Q || op.op_type == OpType::GATE_2Q_LOCAL || op.op_type == OpType::GATE_2Q_GLOBAL) {

            // Bramki 1Q z tego samego węzła DAG mogą być wbijane do logu w pętli (q1_1q_count),
            // więc dodajemy ID tylko, jeśli jest to nowa operacja dla tego kubitu.
            if (op.lq1 != -1 && (actual[op.lq1].empty() || actual[op.lq1].back() != op.gate_id)) {
                actual[op.lq1].push_back(op.gate_id);
            }
            if (op.lq2 != -1 && (actual[op.lq2].empty() || actual[op.lq2].back() != op.gate_id)) {
                actual[op.lq2].push_back(op.gate_id);
            }
        }
    }

    // 3. Porównujemy stan oczekiwany z faktycznym
    bool is_valid = true;
    for (int i = 0; i < num_qubits; ++i) {
        if (expected[i] != actual[i]) {
            is_valid = false;
            std::cout << "[BŁĄD WALIDACJI] Niezgodność kolejności bramek na kubicie logicznym lq" << i << "!\n";

            // Opcjonalny zrzut do łatwiejszego debugowania
            std::cout << "  Oczekiwano ID: ";
            for (int id : expected[i]) std::cout << id << " ";
            std::cout << "\n  Otrzymano ID:  ";
            for (int id : actual[i]) std::cout << id << " ";
            std::cout << "\n";
        }
    }

    if (is_valid) {
        std::cout << "[OK] Walidacja pomyślna. Czystość logiczna (ID bramek z DAG) zachowana na 100%!\n";
    }

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