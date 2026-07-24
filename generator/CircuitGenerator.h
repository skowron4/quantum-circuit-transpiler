//
// Created by skowr on 6.05.2026.
//

#ifndef MAGISTERKA_CIRCUITGENERATOR_H
#define MAGISTERKA_CIRCUITGENERATOR_H


#include <random>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include "../architecture/IQMConstellation.h"
#include "../structures/Gate.h"
#include "structures/Resonator.h"

class CircuitGenerator {
private:
    const BaseTopology& topo;
    int num_physical_qubits;
    int num_logical_qubits;
    const std::vector<Resonator>& resonators;
    std::mt19937 rng;

    // ZOPTYMALIZOWANE: Sprawdza, czy dwa kubity są podpięte pod ten sam rezonator za pomocą std::ranges
    bool share_resonator(int q1, int q2) {
        for (const auto& res : resonators) {
            const auto& connected = res.connected_qubits;
            if (std::ranges::find(connected, q1) != connected.end() &&
                std::ranges::find(connected, q2) != connected.end()) {
                return true;
            }
        }
        return false;
    }

public:
    CircuitGenerator(const BaseTopology& t, int used_logical_qubits, int seed = 42)
    : topo(t),
      num_physical_qubits(t.get_num_qubits()),
      num_logical_qubits(used_logical_qubits),
      resonators(t.get_resonators()),
      rng(seed) {
        if (num_logical_qubits > num_physical_qubits) {
            num_logical_qubits = num_physical_qubits;
        }
    }

    std::vector<Gate> generate(int num_gates, double prob_2q, double prob_global) {
        std::vector<Gate> circuit;
        std::uniform_real_distribution<double> dist_prob(0.0, 1.0);
        std::uniform_int_distribution<int> dist_q(0, num_logical_qubits - 1);
        std::uniform_int_distribution<int> dist_res(0, resonators.size() - 1);

        for (int i = 0; i < num_gates; ++i) {
            double p1 = dist_prob(rng);

            if (p1 > prob_2q) {
                circuit.push_back({i, "1Q", dist_q(rng), -1});
            } else {
                double p2 = dist_prob(rng);

                if (p2 > prob_global) {
                    int r_id = dist_res(rng);
                    const auto& connected = resonators[r_id].connected_qubits;

                    if (connected.size() >= 2) {
                        std::vector<int> q_list(connected.begin(), connected.end());
                        std::vector<int> valid_q_list;
                        for (int q : q_list) {
                            if (q < num_logical_qubits) valid_q_list.push_back(q);
                        }

                        if (valid_q_list.size() >= 2) {
                            std::uniform_int_distribution<int> dist_local_q(0, valid_q_list.size() - 1);
                            int q1 = valid_q_list[dist_local_q(rng)];
                            int q2 = valid_q_list[dist_local_q(rng)];
                            while(q1 == q2) q2 = valid_q_list[dist_local_q(rng)];
                            circuit.push_back({i,"2Q_LOCAL", q1, q2});
                        } else {
                            circuit.push_back({i,"1Q", dist_q(rng), -1});
                        }
                    } else {
                        circuit.push_back({i,"1Q", dist_q(rng), -1});
                    }
                } else {
                    int q1 = dist_q(rng);
                    int q2 = dist_q(rng);
                    while(q1 == q2 || share_resonator(q1, q2)) {
                        q1 = dist_q(rng);
                        q2 = dist_q(rng);
                    }
                    circuit.push_back({i,"2Q_GLOBAL", q1, q2});
                }
            }
        }
        return circuit;
    }

    void print_circuit(const std::vector<Gate>& circuit) const {
        if (circuit.empty()) {
            std::cout << "Obwód jest pusty!!\n";
            return;
        }

        std::vector<std::stringstream> lines(num_physical_qubits);
        for (int i = 0; i < num_physical_qubits; ++i) {
            lines[i] << "q";
            if (i < 10) lines[i] << "0";
            lines[i] << i;

            if (i < num_logical_qubits) {
                lines[i] << " [L]: ";
            } else {
                lines[i] << " [E]: ";
            }
        }

        for (const auto& gate : circuit) {
            int max_len = 0;
            for (const auto& ss : lines) {
                int len = ss.str().length();
                if (len > max_len) max_len = len;
            }

            for (auto& ss : lines) {
                int pad = max_len - ss.str().length();
                for (int p = 0; p < pad; ++p) ss << "-";
            }

            if (gate.type == "1Q") {
                for (int i = 0; i < num_physical_qubits; ++i) {
                    if (i == gate.target_q1) lines[i] << "-[ U3 ]-";
                    else lines[i] << "--------";
                }
            }
            else if (gate.type == "2Q_LOCAL" || gate.type == "2Q_GLOBAL") {
                int c = gate.target_q1;
                int t = gate.target_q2;
                int min_q = std::min(c, t);
                int max_q = std::max(c, t);

                std::string ctrl_char = (gate.type == "2Q_GLOBAL") ? "--*--" : "--o--";
                std::string targ_char = (gate.type == "2Q_GLOBAL") ? "-[+]-" : "-(+)-";

                for (int i = 0; i < num_physical_qubits; ++i) {
                    if (i == c) lines[i] << ctrl_char;
                    else if (i == t) lines[i] << targ_char;
                    else if (i > min_q && i < max_q) lines[i] << "--|--";
                    else lines[i] << "-----";
                }
            }
        }

        // std::cout << "\n=== WIZUALIZACJA OBWODU KWANTOWEGO ===\n";
        // std::cout << "Oznaczenia kubitów: [L] - Logiczny (Zajęty) | [E] - Pusty (Dostępny dla MOVE)\n";
        // std::cout << "Legenda bramek:   --o--/-(+)- [Lokalna] | --*--/-[+]- [Globalna - wymaga routingu]\n\n";
        //
        // for (int i = 0; i < num_physical_qubits; ++i) {
        //     std::cout << lines[i].str() << "\n";
        // }
        // std::cout << "======================================\n\n";
    }

    void link_gates(std::vector<Gate>& circuit) {
        for (auto& g : circuit) {
            g.prev_gate_q1 = g.prev_gate_q2 = nullptr;
            g.next_gate_q1 = g.next_gate_q2 = nullptr;
        }

        std::vector<Gate*> last_gate_on_qubit(this->num_physical_qubits, nullptr);

        for (int i = 0; i < circuit.size(); ++i) {
            Gate* current = &circuit[i];

            if (last_gate_on_qubit[current->target_q1] != nullptr) {
                Gate* prev = last_gate_on_qubit[current->target_q1];
                current->prev_gate_q1 = prev;
                if (prev->target_q1 == current->target_q1) prev->next_gate_q1 = current;
                else prev->next_gate_q2 = current;
            }
            last_gate_on_qubit[current->target_q1] = current;

            if (current->type != "1Q" && current->target_q2 != -1) {
                if (last_gate_on_qubit[current->target_q2] != nullptr) {
                    Gate* prev = last_gate_on_qubit[current->target_q2];
                    current->prev_gate_q2 = prev;
                    if (prev->target_q1 == current->target_q2) prev->next_gate_q1 = current;
                    else prev->next_gate_q2 = current;
                }
                last_gate_on_qubit[current->target_q2] = current;
            }
        }
    }

    void print_gate_statistics(const std::vector<Gate>& circuit) const {
        int count_1q = 0, count_2q_local = 0, count_2q_global = 0;

        for (const auto& gate : circuit) {
            if (gate.type == "1Q") count_1q++;
            else if (gate.type == "2Q_LOCAL") count_2q_local++;
            else if (gate.type == "2Q_GLOBAL") count_2q_global++;
        }

        std::cout << "=== STATYSTYKI OBWODU ===\n";
        std::cout << "Całkowita liczba bramek: " << circuit.size() << "\n";
        std::cout << "Bramki 1-kubitowe: " << count_1q << "\n";
        std::cout << "Bramki 2-kubitowe (suma): " << (count_2q_local + count_2q_global) << "\n";
        std::cout << "  - Lokalne: " << count_2q_local << "\n";
        std::cout << "  - Globalne: " << count_2q_global << "\n";
        std::cout << "=========================\n\n";
    }
};


#endif //MAGISTERKA_CIRCUITGENERATOR_H