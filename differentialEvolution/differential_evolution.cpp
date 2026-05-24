//
// Created by skowr on 18.05.2026.
//

#include "differential_evolution.h"
#include <iostream>
#include <algorithm>
#include <random>
#include "structures/Gate.h"
#include <chrono>
#include <fstream>
#include <queue>
#include <sstream>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/xoshiro.hpp>

class IQMConstellation;
using namespace std;

const int RES_OFFSET = 1000;

DifferentialEvolution::DifferentialEvolution(const BaseTopology& topology, const vector<Gate>& circuit, int qubits_cnt, string name)
    : topo(topology), original_circuit(circuit), num_qubits(qubits_cnt), name(name) {
    preprocess_circuit(circuit);
}

void DifferentialEvolution::preprocess_circuit(const vector<Gate>& circuit) {
    vector<int> last_gate_on_qubit(num_qubits, -1);
    vector<int> pending_1q(num_qubits, 0);
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
                int dep = last_gate_on_qubit[node.lq1];
                grouped_nodes[dep].next_gates.push_back(node.id);
                node.deps_count++;
            }
            if (last_gate_on_qubit[node.lq2] != -1 && last_gate_on_qubit[node.lq2] != last_gate_on_qubit[node.lq1]) {
                int dep = last_gate_on_qubit[node.lq2];
                grouped_nodes[dep].next_gates.push_back(node.id);
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

void DifferentialEvolution::optimize(int pop_size, int generations, int num_logical_qubits, double f, double cr) {
    std::cout << "Rozpoczynam optymalizacje DE... (proszę czekać)\n";

    int n_genes = grouped_nodes.size() * 3;

    // Opcjonalny 4. blok DODAJEMY TYLKO, jeśli są w ogóle jakieś wolne fizyczne kubity
    int logical_q = (num_logical_qubits <= 0) ? num_qubits : num_logical_qubits;

    if (logical_q < topo.qubits.size()) {
        n_genes = grouped_nodes.size() * 4; // lub n_gates * 4 w wersji windowed
        std::cout << "Aktywowano Lazy Off-boarding (wolne kubity fizyczne dostępne).\n";
    }

    if (n_genes == 0) {
        std::cout << "Obwód jest pusty lub preprocesing zawiódł!\n";
        return;
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    // 1. Alokacja pamięci i inicjalizacja generatora
    std::vector<double> pop(pop_size * n_genes);
    std::vector<double> fitness(pop_size);
    boost::random::xoshiro256pp rng(1337);

    // 2. Populacja początkowa
    initialize_population(pop, fitness, pop_size, n_genes, rng);

    // 3. Główna pętla ewolucji
    run_evolution_loop(pop, fitness, pop_size, generations, n_genes, f, cr, rng);

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> execution_time = end_time - start_time;

    // 4. Ewaluacja i eksport najlepszego wyniku
    process_and_export_best_solution(pop, fitness, n_genes, num_logical_qubits, execution_time.count());
}

void DifferentialEvolution::initialize_population(std::vector<double>& pop, std::vector<double>& fitness, int pop_size, int n_genes, boost::random::xoshiro256pp& rng) {
    boost::random::uniform_real_distribution<double> dist_real(0.0, 1.0);
    std::vector<double> trial(n_genes);

    int d_swaps, d_moves, d_depth;
    double d_time;

    for (int i = 0; i < pop_size; ++i) {
        double* __restrict p_pop = pop.data() + i * n_genes;
        double* __restrict p_trial = trial.data();

        for (int j = 0; j < n_genes; ++j) {
            double val = dist_real(rng);
            p_pop[j] = val;
            p_trial[j] = val;
        }

        fitness[i] = evaluate_fitness(trial, d_swaps, d_moves, d_time, d_depth, nullptr);
    }
}

void DifferentialEvolution::run_evolution_loop(std::vector<double>& pop, std::vector<double>& fitness, int pop_size, int generations, int n_genes, double f, double cr, boost::random::xoshiro256pp& rng) {
    boost::random::uniform_int_distribution<int> dist_int(0, pop_size - 1);
    std::vector<double> trial(n_genes);

    int d_swaps, d_moves, d_depth;
    double d_time;

    for (int gen = 0; gen < generations; ++gen) {
        for (int i = 0; i < pop_size; ++i) {
            int a, b, c;
            select_distinct_parents(i, pop_size, rng, dist_int, a, b, c);
            generate_trial_vector(pop, trial, i, a, b, c, n_genes, f, cr, rng);

            double trial_fitness = evaluate_fitness(trial, d_swaps, d_moves, d_time, d_depth, nullptr);
            if (trial_fitness < fitness[i]) {
                std::copy(trial.begin(), trial.end(), pop.begin() + i * n_genes);
                fitness[i] = trial_fitness;
            }
        }
    }
}

void DifferentialEvolution::select_distinct_parents(int current_idx, int pop_size, boost::random::xoshiro256pp& rng, boost::random::uniform_int_distribution<int>& dist_int, int& a, int& b, int& c) {
    do { a = dist_int(rng); } while (a == current_idx);
    do { b = dist_int(rng); } while (b == current_idx || b == a);
    do { c = dist_int(rng); } while (c == current_idx || c == a || c == b);
}

inline double fast_clamp(double val, double min_val, double max_val) {
    if (val < min_val) return min_val;
    if (val > max_val) return max_val;
    return val;
}

void DifferentialEvolution::generate_trial_vector(const std::vector<double>& pop, std::vector<double>& trial, int current_idx, int a, int b, int c, int n_genes, double f, double cr, boost::random::xoshiro256pp& rng) {
    boost::random::uniform_real_distribution<double> dist_real(0.0, 1.0);
    boost::random::uniform_int_distribution<int> dist_gene(0, n_genes - 1);
    int R = dist_gene(rng);

    // Surowe, zastrzeżone wskaźniki odblokowują dla kompilatora wektoryzację SIMD.
    const double* __restrict p_pop = pop.data();
    double* __restrict p_trial = trial.data();

    const double* __restrict p_a = p_pop + a * n_genes;
    const double* __restrict p_b = p_pop + b * n_genes;
    const double* __restrict p_c = p_pop + c * n_genes;
    const double* __restrict p_curr = p_pop + current_idx * n_genes;

    for (int j = 0; j < n_genes; ++j) {
        bool mutate = (dist_real(rng) < cr || j == R);
        double mutated_val = p_a[j] + f * (p_b[j] - p_c[j]);
        // Zastąpienie if/else potężnie wspiera CPU branch predictor (cmov w x86)
        p_trial[j] = mutate ? fast_clamp(mutated_val, 0.0, 1.0) : p_curr[j];
    }
}

void DifferentialEvolution::process_and_export_best_solution(const std::vector<double>& pop, const std::vector<double>& fitness, int n_genes, int num_logical_qubits, double exec_time_sec) {
    // Szybsze znalezienie minimum przy użyciu iteratorów
    auto best_iter = ranges::min_element(fitness);
    int best_idx = std::distance(fitness.begin(), best_iter);
    double final_best = *best_iter;

    // Kopiowanie najlepszego osobnika do osobnego wektora
    std::vector<double> best_individual(pop.begin() + best_idx * n_genes, pop.begin() + (best_idx + 1) * n_genes);

    int final_swaps = 0;
    int final_moves = 0;
    int final_logical_depth = 0;
    double final_physical_time = 0.0;
    std::vector<RoutingRecord> best_moves_log;

    // Ponowna ewaluacja w celu wyciągnięcia pełnych logów (best_moves_log)
    evaluate_fitness(best_individual, final_swaps, final_moves, final_physical_time, final_logical_depth, &best_moves_log);

    print_optimized_circuit(best_moves_log, num_logical_qubits);
    print_statistics(final_best, final_swaps, final_moves, final_physical_time, final_logical_depth, exec_time_sec);

    std::cout << "=== WALIDACJA OBWODU LOGICZNEGO ===\n";
    if (!validate_logical_circuit(original_circuit, best_moves_log)) {
        exit(-1);
    }
    std::cout << "===================================\n\n";

    export_moves_to_file(best_moves_log, "../best_schedule_moves" + name + ".txt");
}

void DifferentialEvolution::print_statistics(double best_score, int swaps, int moves, double physical_time, int depth, double exec_time) {
    std::cout << "\n=== STATYSTYKI PO ZAKONCZENIU OPTYMALIZACJI DE ===\n";
    std::cout << "Najlepszy znaleziony koszt (Funkcja celu z wagami): " << best_score << "\n";
    std::cout << "Wykorzystane operacje SWAP: " << swaps << "\n";
    std::cout << "Wykorzystane operacje MOVE: " << moves << "\n";
    std::cout << "Czas fizyczny obwodu (Makespan / Physical Time): " << physical_time << " [jednostek]\n";
    std::cout << "Głębokość warstwowa (Logical Depth / Layers): " << depth << " [warstw]\n";
    std::cout << "Czas wykonania algorytmu: " << exec_time << " sekund\n";
    std::cout << "==================================================\n\n";
}

inline double fast_max(double a, double b) {
    return (a > b) ? a : b;
}

double DifferentialEvolution::evaluate_fitness(const vector<double>& genes, int& out_swaps, int& out_moves, double &max_physical_time, int &max_logical_depth, std::vector<RoutingRecord>* out_log) {
    // Punkt 2: Obliczamy rozmiar raz i przekazujemy jako flagę
    bool can_lazy_offboard = genes.size() > grouped_nodes.size() * 3;

    // Punkt 3: Rozdzielacz szablonów. Wybiera w czasie kompilacji wersję z logowaniem lub bez.
    if (out_log) return evaluate_fitness_impl<true>(genes, out_swaps, out_moves, max_physical_time, max_logical_depth, out_log, can_lazy_offboard);
    return evaluate_fitness_impl<false>(genes, out_swaps, out_moves, max_physical_time, max_logical_depth, nullptr, can_lazy_offboard);
}

template <bool LOGGING>
double DifferentialEvolution::evaluate_fitness_impl(const vector<double>& genes, int& out_swaps, int& out_moves, double &max_physical_time, int &max_logical_depth, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard) {
    out_swaps = 0; out_moves = 0; max_physical_time = 0.0; max_logical_depth = 0;
    int n_gates = grouped_nodes.size();
    double total_routing_cost = 0.0;

    prepare_workspaces(n_gates, out_log);

    const double* genes_ptr = genes.data();

    auto cmp = [genes_ptr](int left, int right) {
        return genes_ptr[left] < genes_ptr[right];
    };

    ready_queue_workspace.clear();
    if (ready_queue_workspace.capacity() < n_gates) {
        ready_queue_workspace.reserve(n_gates);
    }
    for (int i = 0; i < n_gates; ++i) {
        deps_workspace[i] = grouped_nodes[i].deps_count;
        if (deps_workspace[i] == 0) ready_queue_workspace.push_back(i);
    }
    std::make_heap(ready_queue_workspace.begin(), ready_queue_workspace.end(), cmp);

    while (!ready_queue_workspace.empty()) {
        std::pop_heap(ready_queue_workspace.begin(), ready_queue_workspace.end(), cmp);
        int gate_id = ready_queue_workspace.back();
        ready_queue_workspace.pop_back();

        const GateNode& node = grouped_nodes[gate_id];
        int p1 = state_workspace.logical_to_phys_qubit[node.lq1];

        // Wywołania korzystają z dziedziczonego LOGGING
        apply_1q_gates<LOGGING>(node.id, node.lq1, p1, node.q1_1q_count, genes, n_gates, out_moves, total_routing_cost, out_log, can_lazy_offboard);

        if (node.lq2 == -1) {
            for (int next_id : node.next_gates) {
                deps_workspace[next_id]--;
                if (deps_workspace[next_id] == 0) {
                    ready_queue_workspace.push_back(next_id);
                    std::push_heap(ready_queue_workspace.begin(), ready_queue_workspace.end(), cmp);
                }
            }
            continue;
        }

        int p2 = state_workspace.logical_to_phys_qubit[node.lq2];
        apply_1q_gates<LOGGING>(node.id, node.lq2, p2, node.q2_1q_count, genes, n_gates, out_moves, total_routing_cost, out_log, can_lazy_offboard);

        shared_res_workspace.clear();
        const auto& p2_resonators = topo.qubits[p2].connected_resonators;
        for (int r : topo.qubits[p1].connected_resonators) {
            if (ranges::find(p2_resonators, r) != p2_resonators.end()) {
                shared_res_workspace.push_back(r);
            }
        }

        if (!shared_res_workspace.empty()) {
            execute_local_2q_gate<LOGGING>(node, p1, p2, gate_id, genes, n_gates, out_moves, total_routing_cost, out_log, can_lazy_offboard);
        } else {
            execute_routed_2q_gate<LOGGING>(node, p1, p2, gate_id, genes, n_gates, out_swaps, out_moves, total_routing_cost, out_log, can_lazy_offboard);
        }

        for (int next_id : node.next_gates) {
            deps_workspace[next_id]--;
            if (deps_workspace[next_id] == 0) {
                ready_queue_workspace.push_back(next_id);
                std::push_heap(ready_queue_workspace.begin(), ready_queue_workspace.end(), cmp);
            }
        }
    }

    return finalize_and_calculate_fitness<LOGGING>(genes, n_gates, out_moves, total_routing_cost, max_physical_time, max_logical_depth, out_log, can_lazy_offboard);
}

void DifferentialEvolution::prepare_workspaces(int n_gates, std::vector<RoutingRecord>* out_log) {
    if (out_log) {
        out_log->clear();
        out_log->reserve(n_gates * 3);
    }

    // Inicjalizacja rozmiarów TYLKO RAZ (alokacja stronnicowania)
    if (state_workspace.phys_qubit_to_logical.empty()) {
        state_workspace.phys_qubit_to_logical.resize(topo.qubits.size());
        state_workspace.logical_to_phys_qubit.resize(num_qubits);
        state_workspace.q_on_resonator.resize(topo.qubits.size());
        state_workspace.res_occupant.resize(topo.resonators.size());
        state_workspace.q_ready_time.resize(topo.qubits.size());
        state_workspace.r_ready_time.resize(topo.resonators.size());
        q_logical_depth_workspace.resize(topo.qubits.size());
        deps_workspace.resize(n_gates);
    }

    // Ekstremalnie szybki zapis blokowy bezpośrednio do pamięci
    // Wartość 0xFF daje dokładnie bitowe -1 (dla int w kodzie uzupełnień do 2)
    std::memset(state_workspace.phys_qubit_to_logical.data(), 0xFF, topo.qubits.size() * sizeof(int));
    std::memset(state_workspace.logical_to_phys_qubit.data(), 0xFF, num_qubits * sizeof(int));
    std::memset(state_workspace.q_on_resonator.data(), 0xFF, topo.qubits.size() * sizeof(int));
    std::memset(state_workspace.res_occupant.data(), 0xFF, topo.resonators.size() * sizeof(int));

    // W IEEE 754 zera to binarne same zera, więc możemy użyć memset do czyszczenia wektorów double
    std::memset(state_workspace.q_ready_time.data(), 0, topo.qubits.size() * sizeof(double));
    std::memset(state_workspace.r_ready_time.data(), 0, topo.resonators.size() * sizeof(double));
    std::memset(q_logical_depth_workspace.data(), 0, topo.qubits.size() * sizeof(int));

    // UWAGA: Usunąłem deps_workspace.assign(n_gates, 0)! Wyjaśnienie w punkcie 3.

    for (int i = 0; i < num_qubits; ++i) {
        state_workspace.logical_to_phys_qubit[i] = i;
        state_workspace.phys_qubit_to_logical[i] = i;
    }
}

template <bool LOGGING>
void DifferentialEvolution::force_qubit_down(int& p, int current_gate_id, const std::vector<double>& genes, int n_gates, int& out_moves, double& total_routing_cost, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard) {
    int r = state_workspace.q_on_resonator[p];
    if (r == -1) return;

    int target_p = p;

    // Punkt 2 w akcji - tanie sprawdzenie flagi zamiast .size()
    if (current_gate_id >= 0 && can_lazy_offboard) {
        int valid_cands[32];
        valid_cands[0] = p;
        int valid_count = 1;

        const auto& connected = topo.resonators[r].connected_qubits;

        for (int p_cand : connected) {
            if (p_cand != p && state_workspace.phys_qubit_to_logical[p_cand] == -1 && state_workspace.q_on_resonator[p_cand] == -1) {
                valid_cands[valid_count++] = p_cand;
            }
        }

        if (valid_count > 1) {
            double gene_val = genes[n_gates * 3 + current_gate_id];
            int target_idx = static_cast<int>(gene_val * valid_count);

            // Zabezpieczenie (tak jak miałeś)
            if (target_idx >= valid_count) target_idx = valid_count - 1;

            target_p = valid_cands[target_idx];
        }
    }

    double step_start = fast_max(state_workspace.q_ready_time[p], state_workspace.r_ready_time[r]);
    if (target_p != p) step_start = fast_max(step_start, state_workspace.q_ready_time[target_p]);

    int lq = state_workspace.phys_qubit_to_logical[p];

    // Punkt 3: Kompilator całkowicie USUNIE to sprawdzenie z głównej pętli!
    if constexpr (LOGGING) {
        out_log->push_back({OpType::MOVE, current_gate_id, lq, -1, r + RES_OFFSET, target_p, -1, step_start});
    }
    out_moves++;
    total_routing_cost += MOVE_COST;

    state_workspace.q_ready_time[target_p] = step_start + MOVE_COST;
    state_workspace.r_ready_time[r] = step_start + MOVE_COST;

    if (target_p != p) {
        state_workspace.phys_qubit_to_logical[target_p] = lq;
        state_workspace.logical_to_phys_qubit[lq] = target_p;
        state_workspace.phys_qubit_to_logical[p] = -1;
        q_logical_depth_workspace[target_p] = q_logical_depth_workspace[p] + 1;
        q_logical_depth_workspace[p] = 0;
    } else {
        q_logical_depth_workspace[p] += 1;
    }

    state_workspace.q_on_resonator[p] = -1;
    state_workspace.q_on_resonator[target_p] = -1;
    state_workspace.res_occupant[r] = -1;
    p = target_p;
}

template <bool LOGGING>
void DifferentialEvolution::apply_1q_gates(int gate_id, int lq, int& p, int q_1q_count, const std::vector<double>& genes, int n_gates, int& out_moves, double& total_routing_cost, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard) {
    if (q_1q_count > 0) {
        force_qubit_down<LOGGING>(p, gate_id, genes, n_gates, out_moves, total_routing_cost, out_log, can_lazy_offboard);

        if constexpr (LOGGING) {
            for (int i = 0; i < q_1q_count; ++i) {
                out_log->push_back({OpType::GATE_1Q, gate_id, lq, -1, p, -1, -1, state_workspace.q_ready_time[p] + i});
            }
        }
        state_workspace.q_ready_time[p] += q_1q_count;
        q_logical_depth_workspace[p] += q_1q_count;
    }
}

template <bool LOGGING>
void DifferentialEvolution::execute_local_2q_gate(const GateNode& node, int p1, int p2, int gate_id, const std::vector<double>& genes, int n_gates, int& out_moves, double& total_routing_cost, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard) {
    double start_time = fast_max(state_workspace.q_ready_time[p1], state_workspace.q_ready_time[p2]);

    bool q1_goes_in = (node.q1_1q_count < node.q2_1q_count) ? true :
                      (node.q1_1q_count > node.q2_1q_count) ? false :
                      (genes[n_gates * 2 + gate_id] < 0.5);

    int p_in = q1_goes_in ? p1 : p2;
    int p_stay = q1_goes_in ? p2 : p1;
    int lq_in = q1_goes_in ? node.lq1 : node.lq2;
    int lq_stay = q1_goes_in ? node.lq2 : node.lq1;

    int chosen_res = shared_res_workspace[0];
    if (shared_res_workspace.size() > 1) {
        int res_idx = (int)(genes[n_gates + gate_id] * shared_res_workspace.size());
        chosen_res = shared_res_workspace[std::min(res_idx, (int)shared_res_workspace.size() - 1)];
    }

    force_qubit_down<LOGGING>(p_stay, node.id, genes, n_gates, out_moves, total_routing_cost, out_log, can_lazy_offboard);

    if (state_workspace.res_occupant[chosen_res] != -1 && state_workspace.res_occupant[chosen_res] != p_in) {
        int occ = state_workspace.res_occupant[chosen_res];
        force_qubit_down<LOGGING>(occ, node.id, genes, n_gates, out_moves, total_routing_cost, out_log, can_lazy_offboard);
    }

    if (state_workspace.q_on_resonator[p_in] != chosen_res) {
        force_qubit_down<LOGGING>(p_in, node.id, genes, n_gates, out_moves, total_routing_cost, out_log, can_lazy_offboard);
        start_time = fast_max(state_workspace.q_ready_time[p_in], fast_max(state_workspace.q_ready_time[p_stay], state_workspace.r_ready_time[chosen_res]));

        if constexpr (LOGGING) {
            out_log->push_back({OpType::MOVE, node.id, lq_in, -1, p_in, chosen_res + RES_OFFSET, -1, start_time});
        }
        out_moves++;
        total_routing_cost += MOVE_COST;
        start_time += MOVE_COST;

        state_workspace.q_on_resonator[p_in] = chosen_res;
        state_workspace.res_occupant[chosen_res] = p_in;
        q_logical_depth_workspace[p_in] += 1;
    } else {
        start_time = fast_max(state_workspace.q_ready_time[p_in], fast_max(state_workspace.q_ready_time[p_stay], state_workspace.r_ready_time[chosen_res]));
    }

    if constexpr (LOGGING) {
        out_log->push_back({OpType::GATE_2Q_LOCAL, node.id, lq_in, lq_stay, chosen_res + RES_OFFSET, p_stay, p_in, start_time});
    }

    state_workspace.q_ready_time[p_in] = start_time + 1;
    state_workspace.q_ready_time[p_stay] = start_time + 1;
    state_workspace.r_ready_time[chosen_res] = start_time + 1;

    int gate_layer = fast_max(q_logical_depth_workspace[p_in], q_logical_depth_workspace[p_stay]);
    q_logical_depth_workspace[p_in] = gate_layer + 1;
    q_logical_depth_workspace[p_stay] = gate_layer + 1;
}

template <bool LOGGING>
void DifferentialEvolution::execute_routed_2q_gate(const GateNode& node, int p1, int p2, int gate_id, const std::vector<double>& genes, int n_gates, int& out_swaps, int& out_moves, double& total_routing_cost, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard) {
    force_qubit_down<LOGGING>(p1, node.id, genes, n_gates, out_moves, total_routing_cost, out_log, can_lazy_offboard);
    force_qubit_down<LOGGING>(p2, node.id, genes, n_gates, out_moves, total_routing_cost, out_log, can_lazy_offboard);

    const auto& path = topo.get_routing_path(p1, p2);
    if (path.resonators.empty()) return;

    int current_p = p1;
    process_routing_path<LOGGING>(node, path, current_p, genes, n_gates, out_swaps, out_moves, total_routing_cost, out_log, can_lazy_offboard);

    int final_res = path.resonators.back();
    int g_lq_in_cand = state_workspace.phys_qubit_to_logical[current_p];
    int g_lq_stay_cand = state_workspace.phys_qubit_to_logical[p2];

    int q_in_1q_count = (g_lq_in_cand == node.lq1) ? node.q1_1q_count : node.q2_1q_count;
    int q_stay_1q_count = (g_lq_stay_cand == node.lq1) ? node.q1_1q_count : node.q2_1q_count;

    bool g_q1_goes_in = (q_in_1q_count < q_stay_1q_count) ? true :
                        (q_in_1q_count > q_stay_1q_count) ? false :
                        (genes[n_gates * 2 + gate_id] < 0.5);

    int g_p_in = g_q1_goes_in ? current_p : p2;
    int g_p_stay = g_q1_goes_in ? p2 : current_p;
    int g_lq_in_final = g_q1_goes_in ? g_lq_in_cand : g_lq_stay_cand;
    int g_lq_stay_final = g_q1_goes_in ? g_lq_stay_cand : g_lq_in_cand;

    force_qubit_down<LOGGING>(g_p_stay, node.id, genes, n_gates, out_moves, total_routing_cost, out_log, can_lazy_offboard);

    if (state_workspace.res_occupant[final_res] != -1 && state_workspace.res_occupant[final_res] != g_p_in) {
        int occ = state_workspace.res_occupant[final_res];
        force_qubit_down<LOGGING>(occ, node.id, genes, n_gates, out_moves, total_routing_cost, out_log, can_lazy_offboard);
    }

    double exec_start = fast_max(state_workspace.q_ready_time[current_p], fast_max(state_workspace.q_ready_time[p2], state_workspace.r_ready_time[final_res]));

    if constexpr (LOGGING) {
        out_log->push_back({OpType::MOVE, node.id, g_lq_in_final, -1, g_p_in, final_res + RES_OFFSET, -1, exec_start});
    }

    out_moves++;
    total_routing_cost += MOVE_COST;
    exec_start += MOVE_COST;

    state_workspace.q_on_resonator[g_p_in] = final_res;
    state_workspace.res_occupant[final_res] = g_p_in;
    q_logical_depth_workspace[g_p_in] += 1;

    if constexpr (LOGGING) {
        out_log->push_back({OpType::GATE_2Q_GLOBAL, node.id, g_lq_in_final, g_lq_stay_final, final_res + RES_OFFSET, g_p_stay, g_p_in, exec_start});
    }

    state_workspace.q_ready_time[current_p] = exec_start + 1;
    state_workspace.q_ready_time[p2] = exec_start + 1;
    state_workspace.r_ready_time[final_res] = exec_start + 1;

    int global_gate_layer = fast_max(q_logical_depth_workspace[current_p], q_logical_depth_workspace[p2]);
    q_logical_depth_workspace[current_p] = global_gate_layer + 1;
    q_logical_depth_workspace[p2] = global_gate_layer + 1;
}

template <bool LOGGING>
void DifferentialEvolution::process_routing_path(const GateNode& node, const RoutingPath& path, int& current_p, const std::vector<double>& genes, int n_gates, int& out_swaps, int& out_moves, double& total_routing_cost, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard) {
    for (size_t i = 0; i < path.resonators.size() - 1; ++i) {
        int r_step = path.resonators[i];
        int next_p = path.qubits[i+1];

        int occ_next = next_p;
        force_qubit_down<LOGGING>(occ_next, node.id, genes, n_gates, out_moves, total_routing_cost, out_log, can_lazy_offboard);

        if (state_workspace.res_occupant[r_step] != -1) {
            int occ_res = state_workspace.res_occupant[r_step];
            force_qubit_down<LOGGING>(occ_res, node.id, genes, n_gates, out_moves, total_routing_cost, out_log, can_lazy_offboard);
        }
        double step_start = fast_max(state_workspace.q_ready_time[current_p], fast_max(state_workspace.r_ready_time[r_step], state_workspace.q_ready_time[next_p]));
        bool is_move = (state_workspace.phys_qubit_to_logical[next_p] == -1);
        double op_cost = is_move ? (2.0 * MOVE_COST) : SWAP_COST;
        total_routing_cost += op_cost;

        int lq_mov = state_workspace.phys_qubit_to_logical[current_p];
        int lq_target = state_workspace.phys_qubit_to_logical[next_p];

        if constexpr (LOGGING) {
            if (is_move) {
                out_log->push_back({OpType::MOVE, node.id, lq_mov, -1, current_p, r_step + RES_OFFSET, -1, step_start});
                out_log->push_back({OpType::MOVE, node.id, lq_mov, -1, r_step + RES_OFFSET, next_p, -1, step_start + MOVE_COST});
            } else {
                out_log->push_back({OpType::SWAP, node.id, lq_mov, lq_target, current_p, next_p, r_step, step_start});
            }
        }

        int current_layer = fast_max(q_logical_depth_workspace[current_p], q_logical_depth_workspace[next_p]);
        int layers_added = is_move ? 2 : 1;
        q_logical_depth_workspace[current_p] = current_layer + layers_added;
        q_logical_depth_workspace[next_p] = current_layer + layers_added;

        if (is_move) out_moves += 2;
        else out_swaps++;

        state_workspace.phys_qubit_to_logical[next_p] = lq_mov;
        state_workspace.logical_to_phys_qubit[lq_mov] = next_p;

        if (!is_move && lq_target != -1) {
            state_workspace.phys_qubit_to_logical[current_p] = lq_target;
            state_workspace.logical_to_phys_qubit[lq_target] = current_p;
        } else {
            state_workspace.phys_qubit_to_logical[current_p] = -1;
        }

        state_workspace.q_ready_time[current_p] = step_start + op_cost;
        state_workspace.r_ready_time[r_step] = step_start + op_cost;
        state_workspace.q_ready_time[next_p] = step_start + op_cost;

        current_p = next_p;
    }
}

template <bool LOGGING>
double DifferentialEvolution::finalize_and_calculate_fitness(const std::vector<double>& genes, int n_gates, int& out_moves, double& total_routing_cost, double& max_physical_time, int& max_logical_depth, std::vector<RoutingRecord>* out_log, bool can_lazy_offboard) {
    for (int p = 0; p < topo.qubits.size(); ++p) {
        int temp_p = p;
        force_qubit_down<LOGGING>(temp_p, -1, genes, n_gates, out_moves, total_routing_cost, out_log, can_lazy_offboard);
    }

    // Punkt 4: Sortujemy tylko wtedy, gdy faktycznie są jakieś logi
    if constexpr (LOGGING) {
        if (out_log && !out_log->empty()) {
            std::stable_sort(out_log->begin(), out_log->end(), [](const RoutingRecord& a, const RoutingRecord& b) {
                return a.start_time < b.start_time;
            });
        }
    }

    for (double t : state_workspace.q_ready_time) max_physical_time = fast_max(max_physical_time, t);
    for (int d : q_logical_depth_workspace) max_logical_depth = fast_max(max_logical_depth, d);

    double w_routing = 1000.0;
    double w_physical_time = 1.0;
    double w_logical_depth = 0.5;

    return (total_routing_cost * w_routing) + (max_physical_time * w_physical_time) + (max_logical_depth * w_logical_depth);
}

void DifferentialEvolution::print_optimized_circuit(const std::vector<RoutingRecord>& log, int num_logical_qubits) const {
    if (num_logical_qubits <= 0) num_logical_qubits = num_qubits;
    if (log.empty()) { cout << "Zoptymalizowany obwód jest pusty!!\n"; return; }

    int n_q = topo.qubits.size();
    vector<stringstream> lines_q(n_q);

    // Inicjalizacja nagłówków wyłącznie dla kubitów fizycznych
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
            int max_q = fast_max(l1, l2);
            for (int i = 0; i < n_q; ++i) {
                if (i == l1 || i == l2) lines_q[i] << "--x--";
                else if (i > min_q && i < max_q) lines_q[i] << "--|--";
                else lines_q[i] << "-----";
            }
        } else if (move.op_type == OpType::MOVE) {
            if (l1 < RES_OFFSET && l2 < RES_OFFSET && l2 >= 0) {
                // Routing MOVE między dwoma fizycznymi kubitami (jeden jest pusty)
                int min_q = std::min(l1, l2);
                int max_q = fast_max(l1, l2);
                for (int i = 0; i < n_q; ++i) {
                    if (i == l1 || i == l2) lines_q[i] << "--m--";
                    else if (i > min_q && i < max_q) lines_q[i] << "--|--";
                    else lines_q[i] << "-----";
                }
            } else {
                // Wejście lub Zejście z rezonatora przed bramką 2Q
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

            // ZAMIENIAMY CAŁĄ PĘTLE WYSZUKUJĄCĄ NA ODCZYT Z NOWEGO POLA
            int q_in = move.res_used; // <- Nowa niezawodna mechanika DE

            int min_q = std::min(q_in, q_stay);
            int max_q = fast_max(q_in, q_stay);

            bool q_in_is_control = (move.lq1 == grouped_nodes[move.gate_id].lq1);

            string ctrl_char = (move.op_type == OpType::GATE_2Q_LOCAL) ? "--o--" : "--*--";
            string targ_char = (move.op_type == OpType::GATE_2Q_LOCAL) ? "-(+)-" : "-[+]-";

            for (int i = 0; i < n_q; ++i) {
                if (i == q_in) {
                    // Jeśli kubit wchodzący był kontrolą, rysujemy mu o/*. Jeśli nie, rysujemy (+)/[+]
                    lines_q[i] << (q_in_is_control ? ctrl_char : targ_char);
                }
                else if (i == q_stay) {
                    // Analogicznie dla kubitu zostającego
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

    cout << "\n=== WIZUALIZACJA ZOPTYMALIZOWANEGO OBWODU (DE) ===\n";
    cout << "Legenda bramek: -[ U3 ]- [1Q] | --o--/-(+)- [2Q Lokalna] | --*--/-[+]- [2Q Globalna]\n";
    cout << "Operacje: --x-- [SWAP] | --m-- [MOVE / Wejscie na rezonator]\n\n";
    for (int i = 0; i < n_q; ++i) { cout << lines_q[i].str() << "\n"; }
    cout << "==================================================\n\n";
}

void DifferentialEvolution::export_moves_to_file(const std::vector<RoutingRecord>& moves, const std::string& filename) {
    if (moves.empty()) return;

    std::ofstream out_file(filename);
    if (!out_file.is_open()) return;

    // Dodano kolumny time_step oraz res_used
    out_file << "time_step,op_type,gate_id,lq1,lq2,loc1,loc2,res_used\n";

    double current_time = -1.0;
    int time_step = -1;

    for (const auto& move : moves) {
        // Grupowanie współbieżnych operacji (tolerancja na błędy zmiennoprzecinkowe)
        if (std::abs(move.start_time - current_time) > 1e-5) {
            time_step++;
            current_time = move.start_time;
        }

        out_file << time_step << ","
                 << to_string(move.op_type) << ","
                 << move.gate_id << ","
                 << move.lq1 << ","
                 << move.lq2 << ","
                 << move.loc1 << ","
                 << move.loc2 << ","
                 << move.res_used << "\n";
    }

    out_file.close();
    std::cout << ">>> Pomyślnie wyeksportowano współbieżny harmonogram do: " << filename << "\n";
}

bool DifferentialEvolution::validate_logical_circuit(const std::vector<Gate>& orig_circ, const std::vector<RoutingRecord>& log) const {
    std::vector<std::vector<std::string>> expected(num_qubits);
    for (const auto& g : orig_circ) {
        if (g.type == "1Q") {
            expected[g.target_q1].push_back("1Q");
        } else {
            expected[g.target_q1].push_back("2Q_" + std::to_string(g.target_q2));
            expected[g.target_q2].push_back("2Q_" + std::to_string(g.target_q1));
        }
    }

    std::vector<std::vector<std::string>> actual(num_qubits);
    for (const auto& op : log) {
        if (op.op_type == OpType::GATE_1Q) {
            actual[op.lq1].push_back("1Q");
        } else if (op.op_type == OpType::GATE_2Q_LOCAL || op.op_type == OpType::GATE_2Q_GLOBAL) {
            actual[op.lq1].push_back("2Q_" + std::to_string(op.lq2));
            actual[op.lq2].push_back("2Q_" + std::to_string(op.lq1));
        }
    }

    bool is_valid = true;
    for (int i = 0; i < num_qubits; ++i) {
        if (expected[i] != actual[i]) {
            is_valid = false;
            std::cout << "[BŁĄD WALIDACJI] Niezgodność na kubicie logicznym lq" << i << "!\n";
            std::cout << "  Oczekiwano: ";
            for (const auto& s : expected[i]) std::cout << s << " ";
            std::cout << "\n  Otrzymano:  ";
            for (const auto& s : actual[i]) std::cout << s << " ";
            std::cout << "\n";
        }
    }

    if (is_valid) {
        std::cout << "[OK] Walidacja pomyślna. Obwód wyjściowy zachowuje semantykę logiczną obwodu wejściowego.\n";
    }
    return is_valid;
}

void DifferentialEvolution::optimize_windowed(int pop_size, int generations_per_window, int window_size, double accept_ratio, int num_logical_qubits, double f, double cr) {
    std::cout << "Rozpoczynam optymalizacje DE (Oknem przesuwnym: " << window_size
              << " bramek, akceptacja: " << (accept_ratio * 100) << "%)... (proszę czekać)\n";

    int n_gates = grouped_nodes.size();
    int n_genes = n_gates * 3;
    int logical_q = (num_logical_qubits <= 0) ? num_qubits : num_logical_qubits;

    if (logical_q < topo.qubits.size()) {
        n_genes = grouped_nodes.size() * 4; // lub n_gates * 4 w wersji windowed
        std::cout << "Aktywowano Lazy Off-boarding (wolne kubity fizyczne dostępne).\n";
    }

    if (n_genes == 0) {
        std::cout << "Obwód jest pusty lub preprocesing zawiódł!\n";
        return;
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    // 1. Inicjalizacja populacji - wszystkie geny dostają randoma na start
    std::vector<double> pop(pop_size * n_genes);
    std::vector<double> fitness(pop_size);
    boost::random::xoshiro256pp rng(1337);

    initialize_population(pop, fitness, pop_size, n_genes, rng);

    // Wyliczamy ile bramek zamrażamy w każdym kroku (np. 50 z okna 100)
    int step = std::max(1, static_cast<int>(window_size * accept_ratio));

    // 2. Ślizgamy się po obwodzie
    for (int w_start = 0; w_start < n_gates; w_start += step) {
        int w_end = std::min(w_start + window_size, n_gates);
        std::cout << ">>> Przetwarzanie okna: bramki [" << w_start << " - " << w_end - 1 << "] / " << n_gates << "\n";

        // Odpalamy ewolucję tylko dla genów z aktualnego okna
        run_evolution_loop_windowed(pop, fitness, pop_size, generations_per_window, n_genes, n_gates, f, cr, rng, w_start, w_end);

        // Znajdujemy najlepszego mocarza z tego okna
        auto best_iter = ranges::min_element(fitness);
        int best_idx = std::distance(fitness.begin(), best_iter);

        // Kopiujemy jego ZAACKEPTOWANE decyzje (od w_start do w_start + step) do WSZYSTKICH osobników.
        // W ten sposób cała populacja ma wspólną zoptymalizowaną przeszłość i skupia się tylko na przyszłości.
        int accept_end = std::min(w_start + step, n_gates);
        int num_blocks = n_genes / n_gates;

        for (int i = 0; i < pop_size; ++i) {
            if (i == best_idx) continue;
            for (int blk = 0; blk < num_blocks; ++blk) {
                for (int g = w_start; g < accept_end; ++g) {
                    int gene_idx = blk * n_gates + g;
                    pop[i * n_genes + gene_idx] = pop[best_idx * n_genes + gene_idx];
                }
            }
            // Zmieniliśmy im geny (przeszłość), więc musimy przeliczyć ich fitness na nowy start
            std::vector<double> ind(pop.begin() + i * n_genes, pop.begin() + (i + 1) * n_genes);
            int ds, dm, dd; double dt;
            fitness[i] = evaluate_fitness(ind, ds, dm, dt, dd, nullptr);
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> execution_time = end_time - start_time;

    // 3. Eksport całości
    process_and_export_best_solution(pop, fitness, n_genes, num_logical_qubits, execution_time.count());
}

void DifferentialEvolution::run_evolution_loop_windowed(std::vector<double>& pop, std::vector<double>& fitness, int pop_size, int generations, int n_genes, int n_gates, double f, double cr, boost::random::xoshiro256pp& rng, int w_start, int w_end) {
    boost::random::uniform_int_distribution<int> dist_int(0, pop_size - 1);
    std::vector<double> trial(n_genes);

    int d_swaps, d_moves, d_depth;
    double d_time;

    for (int gen = 0; gen < generations; ++gen) {
        for (int i = 0; i < pop_size; ++i) {
            int a, b, c;
            select_distinct_parents(i, pop_size, rng, dist_int, a, b, c);

            // Tutaj magia - mutujemy tylko to co jest w oknie [w_start, w_end)
            generate_trial_vector_windowed(pop, trial, i, a, b, c, n_genes, n_gates, f, cr, rng, w_start, w_end);

            double trial_fitness = evaluate_fitness(trial, d_swaps, d_moves, d_time, d_depth, nullptr);
            if (trial_fitness < fitness[i]) {
                std::copy(trial.begin(), trial.end(), pop.begin() + i * n_genes);
                fitness[i] = trial_fitness;
            }
        }
    }
}

void DifferentialEvolution::generate_trial_vector_windowed(const std::vector<double>& pop, std::vector<double>& trial, int current_idx, int a, int b, int c, int n_genes, int n_gates, double f, double cr, boost::random::xoshiro256pp& rng, int w_start, int w_end) {
    boost::random::uniform_real_distribution<double> dist_real(0.0, 1.0);

    // Losujemy indeks 'R' wymuszający mutację, ale musi on spaść z puli aktywnych genów!
    int num_blocks = n_genes / n_gates;
    int active_size = (w_end - w_start) * num_blocks;
    boost::random::uniform_int_distribution<int> dist_active(0, active_size - 1);

    int R_idx = dist_active(rng);
    int R_block = R_idx / (w_end - w_start);
    int R_gate = w_start + (R_idx % (w_end - w_start));
    int R = R_block * n_gates + R_gate; // Rzeczywisty indeks upewniający, że wylosowaliśmy aktywny gen

    const double* __restrict p_pop = pop.data();
    double* __restrict p_trial = trial.data();
    const double* __restrict p_a = p_pop + a * n_genes;
    const double* __restrict p_b = p_pop + b * n_genes;
    const double* __restrict p_c = p_pop + c * n_genes;
    const double* __restrict p_curr = p_pop + current_idx * n_genes;

    for (int j = 0; j < n_genes; ++j) {
        int gate_id = j % n_gates; // Znajdujemy, do której bramki należy ten konkretny gen

        if (gate_id >= w_start && gate_id < w_end) {
            // Gen w oknie - lecimy normalną mutacją DE
            bool mutate = (dist_real(rng) < cr || j == R);
            double mutated_val = p_a[j] + f * (p_b[j] - p_c[j]);
            p_trial[j] = mutate ? fast_clamp(mutated_val, 0.0, 1.0) : p_curr[j];
        } else {
            // Gen POZA oknem (zamrożony z przeszłości lub czekający w przyszłości) - przepisujemy bez zmian
            p_trial[j] = p_curr[j];
        }
    }
}

template double DifferentialEvolution::evaluate_fitness_impl<true>(const std::vector<double>&, int&, int&, double&, int&, std::vector<RoutingRecord>*, bool);
template double DifferentialEvolution::evaluate_fitness_impl<false>(const std::vector<double>&, int&, int&, double&, int&, std::vector<RoutingRecord>*, bool);