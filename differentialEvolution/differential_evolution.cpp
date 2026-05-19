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

double DifferentialEvolution::evaluate_fitness(const vector<double>& genes, int& out_swaps, int& out_moves, double &max_physical_time, int &max_logical_depth, std::vector<RoutingRecord>* out_log) {
    out_swaps = 0;
    out_moves = 0;
    max_physical_time = 0.0;
    max_logical_depth = 0;
    if (out_log) out_log->clear();

    int n_gates = grouped_nodes.size();
    const double MOVE_COST = 1.0;
    const double SWAP_COST = 3.14;
    double total_routing_cost = 0.0;

    MachineState state;
    state.phys_qubit_to_logical.assign(topo.qubits.size(), -1);
    state.logical_to_phys_qubit.assign(num_qubits, -1);
    state.q_on_resonator.assign(topo.qubits.size(), -1);
    state.res_occupant.assign(topo.resonators.size(), -1);
    state.q_ready_time.assign(topo.qubits.size(), 0.0);
    state.r_ready_time.assign(topo.resonators.size(), 0.0);

    // NOWE: Tablica śledząca warstwy (głębokość logiczną) dla każdego kubitu fizycznego
    std::vector<int> q_logical_depth(topo.qubits.size(), 0);

    auto force_qubit_down = [&](int p, int current_gate_id) {
        if (state.q_on_resonator[p] != -1) {
            int r = state.q_on_resonator[p];
            double step_start = std::max(state.q_ready_time[p], state.r_ready_time[r]);

            if (out_log) {
                int lq = state.phys_qubit_to_logical[p];
                out_log->push_back({"MOVE", current_gate_id, lq, -1, r + RES_OFFSET, p, -1, step_start});
            }
            out_moves++;
            total_routing_cost += MOVE_COST;

            state.q_ready_time[p] = step_start + MOVE_COST;
            state.r_ready_time[r] = step_start + MOVE_COST;

            q_logical_depth[p] += 1; // Zrzucenie zajmuje warstwę operacyjną

            state.q_on_resonator[p] = -1;
            state.res_occupant[r] = -1;
        }
    };

    for (int i = 0; i < num_qubits; ++i) {
        state.logical_to_phys_qubit[i] = i;
        state.phys_qubit_to_logical[i] = i;
    }

    auto cmp = [&genes](int left, int right) { return genes[left] < genes[right]; };
    priority_queue<int, std::vector<int>, decltype(cmp)> ready_queue(cmp);

    std::vector<int> deps = std::vector<int>(n_gates);
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

            if (out_log) {
                for (int i = 0; i < node.q1_1q_count; ++i) {
                    out_log->push_back({"1Q", node.id, node.lq1, -1, p1, -1, -1, state.q_ready_time[p1] + i});
                }
            }
            state.q_ready_time[p1] += node.q1_1q_count;
            q_logical_depth[p1] += node.q1_1q_count; // Dodaj warstwy bramek 1Q

            for (int next_id : node.next_gates) {
                deps[next_id]--;
                if (deps[next_id] == 0) ready_queue.push(next_id);
            }
            continue;
        }

        int p2 = state.logical_to_phys_qubit[node.lq2];

        if (node.q1_1q_count > 0) force_qubit_down(p1, node.id);
        if (node.lq2 != -1 && node.q2_1q_count > 0) force_qubit_down(p2, node.id);

        if (out_log) {
            for (int i = 0; i < node.q1_1q_count; ++i)
                out_log->push_back({"1Q", node.id, node.lq1, -1, p1, -1, -1, state.q_ready_time[p1] + i});
            for (int i = 0; i < node.q2_1q_count; ++i)
                out_log->push_back({"1Q", node.id, node.lq2, -1, p2, -1, -1, state.q_ready_time[p2] + i});
        }

        state.q_ready_time[p1] += node.q1_1q_count;
        state.q_ready_time[p2] += node.q2_1q_count;
        q_logical_depth[p1] += node.q1_1q_count; // Aktualizacja głębokości logicznej 1Q
        q_logical_depth[p2] += node.q2_1q_count;

        std::vector<int> shared_res;
        for (int r : topo.qubits[p1].connected_resonators) {
            if (topo.qubits[p2].connected_resonators.count(r)) shared_res.push_back(r);
        }

        double start_time = std::max(state.q_ready_time[p1], state.q_ready_time[p2]);
        int gate_duration = 1;

        bool q1_goes_in = false;
        if (node.q1_1q_count < node.q2_1q_count) q1_goes_in = true;
        else if (node.q1_1q_count > node.q2_1q_count) q1_goes_in = false;
        else q1_goes_in = (genes[n_gates * 2 + gate_id] < 0.5);

        int p_in, p_stay, lq_in, lq_stay;
        if (q1_goes_in) {
            p_in = p1; p_stay = p2; lq_in = node.lq1; lq_stay = node.lq2;
        } else {
            p_in = p2; p_stay = p1; lq_in = node.lq2; lq_stay = node.lq1;
        }

        if (!shared_res.empty()) {
            int chosen_res = shared_res[0];
            if (shared_res.size() > 1) {
                int res_idx = (int)(genes[n_gates + gate_id] * shared_res.size());
                chosen_res = shared_res[std::min(res_idx, (int)shared_res.size() - 1)];
            }

            force_qubit_down(p_stay, node.id);

            if (state.res_occupant[chosen_res] != -1 && state.res_occupant[chosen_res] != p_in) {
                force_qubit_down(state.res_occupant[chosen_res], node.id);
            }

            if (state.q_on_resonator[p_in] != chosen_res) {
                force_qubit_down(p_in, node.id);

                start_time = std::max({state.q_ready_time[p_in], state.q_ready_time[p_stay], state.r_ready_time[chosen_res]});

                if (out_log) out_log->push_back({"MOVE", node.id, lq_in, -1, p_in, chosen_res + RES_OFFSET, -1, start_time});
                out_moves++;
                total_routing_cost += MOVE_COST;
                start_time += MOVE_COST;

                state.q_on_resonator[p_in] = chosen_res;
                state.res_occupant[chosen_res] = p_in;
                q_logical_depth[p_in] += 1; // Podniesienie zajmuje warstwę
            } else {
                start_time = std::max({state.q_ready_time[p_in], state.q_ready_time[p_stay], state.r_ready_time[chosen_res]});
            }

            if (out_log) {
                out_log->push_back({"2Q_LOCAL", node.id, lq_in, lq_stay, chosen_res + RES_OFFSET, p_stay, p_in, start_time});
            }

            state.q_ready_time[p_in] = start_time + gate_duration;
            state.q_ready_time[p_stay] = start_time + gate_duration;
            state.r_ready_time[chosen_res] = start_time + gate_duration;

            // Rejestracja warstwy bramki 2Q
            int gate_layer = std::max(q_logical_depth[p_in], q_logical_depth[p_stay]);
            q_logical_depth[p_in] = gate_layer + 1;
            q_logical_depth[p_stay] = gate_layer + 1;

        } else {
            force_qubit_down(p1, node.id);
            force_qubit_down(p2, node.id);

            auto path = topo.get_routing_path(p1, p2);
            int current_p = p1;

            if (path.resonators.empty()) continue;

            for (size_t i = 0; i < path.resonators.size() - 1; ++i) {
                int r_step = path.resonators[i];
                int next_p = path.qubits[i+1];

                force_qubit_down(next_p, node.id);

                double step_start = std::max({state.q_ready_time[current_p], state.r_ready_time[r_step], state.q_ready_time[next_p]});
                bool is_move = (state.phys_qubit_to_logical[next_p] == -1);

                double op_cost = is_move ? MOVE_COST : SWAP_COST;
                total_routing_cost += op_cost;

                int lq_mov = state.phys_qubit_to_logical[current_p];
                int lq_target = state.phys_qubit_to_logical[next_p];

                if (out_log) {
                    if (is_move) {
                        out_log->push_back({"MOVE", node.id, lq_mov, -1, current_p, r_step + RES_OFFSET, -1, step_start});
                        out_log->push_back({"MOVE", node.id, lq_mov, -1, r_step + RES_OFFSET, next_p, -1, step_start + (op_cost / 2.0)});
                    } else {
                        out_log->push_back({"SWAP", node.id, lq_mov, lq_target, current_p, next_p, r_step, step_start});
                    }
                }

                // Aktualizacja warstwy logicznej dla SWAP / MOVE
                int current_layer = std::max(q_logical_depth[current_p], q_logical_depth[next_p]);
                int layers_added = is_move ? 2 : 1; // Podwójny MOVE kosztuje 2 operacje
                q_logical_depth[current_p] = current_layer + layers_added;
                q_logical_depth[next_p] = current_layer + layers_added;

                if (is_move) out_moves += 2;
                else out_swaps++;

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

            bool g_q1_goes_in = false;
            int g_lq_in_cand = state.phys_qubit_to_logical[current_p];
            int g_lq_stay_cand = state.phys_qubit_to_logical[p2];

            int q_in_1q_count = (g_lq_in_cand == node.lq1) ? node.q1_1q_count : node.q2_1q_count;
            int q_stay_1q_count = (g_lq_stay_cand == node.lq1) ? node.q1_1q_count : node.q2_1q_count;

            if (q_in_1q_count < q_stay_1q_count) g_q1_goes_in = true;
            else if (q_in_1q_count > q_stay_1q_count) g_q1_goes_in = false;
            else g_q1_goes_in = (genes[n_gates * 2 + gate_id] < 0.5);

            int g_p_in, g_p_stay, g_lq_in_final, g_lq_stay_final;
            if (g_q1_goes_in) {
                g_p_in = current_p; g_p_stay = p2;
                g_lq_in_final = g_lq_in_cand; g_lq_stay_final = g_lq_stay_cand;
            } else {
                g_p_in = p2; g_p_stay = current_p;
                g_lq_in_final = g_lq_stay_cand; g_lq_stay_final = g_lq_in_cand;
            }
            force_qubit_down(g_p_stay, node.id);

            if (state.res_occupant[final_res] != -1 && state.res_occupant[final_res] != g_p_in) {
                force_qubit_down(state.res_occupant[final_res], node.id);
            }

            double exec_start = std::max({state.q_ready_time[current_p], state.q_ready_time[p2], state.r_ready_time[final_res]});

            if (out_log) out_log->push_back({"MOVE", node.id, g_lq_in_final, -1, g_p_in, final_res + RES_OFFSET, -1, exec_start});
            out_moves++;
            total_routing_cost += MOVE_COST;
            exec_start += MOVE_COST;

            state.q_on_resonator[g_p_in] = final_res;
            state.res_occupant[final_res] = g_p_in;
            q_logical_depth[g_p_in] += 1; // Podniesienie na finalny rezonator

            if (out_log) out_log->push_back({"2Q_GLOBAL", node.id, g_lq_in_final, g_lq_stay_final, final_res + RES_OFFSET, g_p_stay, g_p_in, exec_start});

            state.q_ready_time[current_p] = exec_start + gate_duration;
            state.q_ready_time[p2] = exec_start + gate_duration;
            state.r_ready_time[final_res] = exec_start + gate_duration;

            // Rejestracja globalnej bramki 2Q logicznie
            int global_gate_layer = std::max(q_logical_depth[current_p], q_logical_depth[p2]);
            q_logical_depth[current_p] = global_gate_layer + 1;
            q_logical_depth[p2] = global_gate_layer + 1;
        }

        for (int next_id : node.next_gates) {
            deps[next_id]--;
            if (deps[next_id] == 0) ready_queue.push(next_id);
        }
    }

    for (int p = 0; p < topo.qubits.size(); ++p) {
        force_qubit_down(p, -1);
    }

    if (out_log) {
        std::stable_sort(out_log->begin(), out_log->end(), [](const RoutingRecord& a, const RoutingRecord& b) {
            return a.start_time < b.start_time;
        });
    }

    for (double t : state.q_ready_time) max_physical_time = std::max(max_physical_time, t);
    for (int d : q_logical_depth) max_logical_depth = std::max(max_logical_depth, d);

    // Modyfikacja funkcji celu (Fitness): Suma ważona
    // Mozesz dostosowac wagi do swoich potrzeb badawczych
    double w_routing = 1000.0;     // Priorytet komunikacji (routing)
    double w_physical_time = 1.0;  // Priorytet Makespan (względem dekoherencji)
    double w_logical_depth = 0.5;  // Mniejsza waga, ale faworyzuje mniejszą liczbę warstw bramek

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

        if (move.op_type == "1Q") {
            if (l1 >= 0 && l1 < n_q) lines_q[l1] << "-[ U3 ]-";
            for (int i = 0; i < n_q; ++i) {
                if (i != l1) lines_q[i] << "--------";
            }
        } else if (move.op_type == "SWAP") {
            int min_q = std::min(l1, l2);
            int max_q = std::max(l1, l2);
            for (int i = 0; i < n_q; ++i) {
                if (i == l1 || i == l2) lines_q[i] << "--x--";
                else if (i > min_q && i < max_q) lines_q[i] << "--|--";
                else lines_q[i] << "-----";
            }
        } else if (move.op_type == "MOVE") {
            if (l1 < RES_OFFSET && l2 < RES_OFFSET && l2 >= 0) {
                // Routing MOVE między dwoma fizycznymi kubitami (jeden jest pusty)
                int min_q = std::min(l1, l2);
                int max_q = std::max(l1, l2);
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
        } else if (move.op_type == "2Q_LOCAL" || move.op_type == "2Q_GLOBAL") {
            int r_idx = -1, q_stay = -1;
            if (l1 >= RES_OFFSET) { r_idx = l1 - RES_OFFSET; q_stay = l2; }
            else if (l2 >= RES_OFFSET) { r_idx = l2 - RES_OFFSET; q_stay = l1; }

            // ZAMIENIAMY CAŁĄ PĘTLE WYSZUKUJĄCĄ NA ODCZYT Z NOWEGO POLA
            int q_in = move.res_used; // <- Nowa niezawodna mechanika DE

            int min_q = std::min(q_in, q_stay);
            int max_q = std::max(q_in, q_stay);

            bool q_in_is_control = (move.lq1 == grouped_nodes[move.gate_id].lq1);

            string ctrl_char = (move.op_type == "2Q_LOCAL") ? "--o--" : "--*--";
            string targ_char = (move.op_type == "2Q_LOCAL") ? "-(+)-" : "-[+]-";

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

void export_moves_to_file(const std::vector<RoutingRecord>& moves, const std::string& filename) {
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
                 << move.op_type << ","
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

void DifferentialEvolution::optimize(int pop_size, int generations, int num_logical_qubits) {
    using namespace std::chrono;

    cout << "Rozpoczynam optymalizacje DE... (proszę czekać)" << endl;

    auto start_time = high_resolution_clock::now();

    int n_genes = grouped_nodes.size() * 3;

    if (n_genes == 0) {
        cout << "Obwód jest pusty lub preprocesing zawiódł!" << endl;
        return;
    }

    vector<vector<double>> pop(pop_size, vector<double>(n_genes));
    vector<double> fitness(pop_size);
    mt19937 rng(1337);
    uniform_real_distribution<double> dist(0.0, 1.0);

    int dummy_swaps, dummy_moves, dummy_logical_depth;
    double dummy_physical_time = 0.0;

    for (int i = 0; i < pop_size; ++i) {
        for (int j = 0; j < n_genes; ++j) pop[i][j] = dist(rng);
        fitness[i] = evaluate_fitness(pop[i], dummy_swaps, dummy_moves, dummy_physical_time, dummy_logical_depth);
    }

    double f = 0.8;
    double cr = 0.9;

    for (int gen = 0; gen < generations; ++gen) {
        for (int i = 0; i < pop_size; ++i) {
            int a = i, b = i, c = i;
            while(a == i) a = rng() % pop_size;
            while(b == i || b == a) b = rng() % pop_size;
            while(c == i || c == a || c == b) c = rng() % pop_size;

            vector<double> trial = pop[i];
            int R = rng() % n_genes;

            for (int j = 0; j < n_genes; ++j) {
                if (dist(rng) < cr || j == R) {
                    trial[j] = pop[a][j] + f * (pop[b][j] - pop[c][j]);
                    if (trial[j] < 0.0) trial[j] = 0.0;
                    if (trial[j] > 1.0) trial[j] = 1.0;
                }
            }

            double trial_fitness = evaluate_fitness(trial, dummy_swaps, dummy_moves, dummy_physical_time, dummy_logical_depth);
            if (trial_fitness < fitness[i]) {
                pop[i] = trial;
                fitness[i] = trial_fitness;
            }
        }
    }

    auto end_time = high_resolution_clock::now();
    duration<double> execution_time = end_time - start_time;

    int best_idx = distance(fitness.begin(), min_element(fitness.begin(), fitness.end()));
    double final_best = fitness[best_idx];

    int final_swaps = 0;
    int final_moves = 0;
    int final_logical_depth = 0;
    double final_physical_time = 0;
    std::vector<RoutingRecord> best_moves_log;

    evaluate_fitness(pop[best_idx], final_swaps, final_moves, final_physical_time, final_logical_depth, &best_moves_log);
    print_optimized_circuit(best_moves_log, num_logical_qubits);

    cout << "\n=== STATYSTYKI PO ZAKONCZENIU OPTYMALIZACJI ===" << endl;
    cout << "Najlepszy znaleziony koszt (Funkcja celu z wagami): " << final_best << endl;
    cout << "Wykorzystane operacje SWAP: " << final_swaps << endl;
    cout << "Wykorzystane operacje MOVE: " << final_moves << endl;
    cout << "Czas fizyczny obwodu (Makespan / Physical Time): " << final_physical_time << " [jednostek]" << endl;
    cout << "Głębokość warstwowa (Logical Depth / Layers): " << final_logical_depth << " [warstw]" << endl;
    cout << "Czas wykonania algorytmu: " << execution_time.count() << " sekund" << endl;
    cout << "===============================================\n" << endl;

    cout << "=== WALIDACJA OBWODU LOGICZNEGO ===" << endl;
    if (!validate_logical_circuit(original_circuit, best_moves_log)) exit(-1);
    cout << "===================================\n" << endl;

    export_moves_to_file(best_moves_log, "../best_schedule_moves" + name + ".txt");
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
        if (op.op_type == "1Q") {
            actual[op.lq1].push_back("1Q");
        } else if (op.op_type == "2Q_LOCAL" || op.op_type == "2Q_GLOBAL") {
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