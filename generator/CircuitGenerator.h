//
// Created by skowr on 6.05.2026.
//

#ifndef MAGISTERKA_CIRCUITGENERATOR_H
#define MAGISTERKA_CIRCUITGENERATOR_H


#include <random> // Dodaj to na samej górze pliku dla RNG
#include <sstream> // Wymagane dla stringstream
#include <iomanip>
#include <iostream>
#include <vector>
#include <climits>
#include "../architecture/IQMConstellation.h"
#include "../structures/Gate.h"

// // Reprezentacja pojedynczej bramki w obwodzie
// struct Gate {
//     std::string type;     // "1Q", "2Q_LOCAL", "2Q_GLOBAL"
//     int target_q1;
//     int target_q2;   // -1 dla bramek 1-kubitowych
// };

class CircuitGenerator {
private:
    int num_qubits;
    const std::vector<Resonator>& resonators;
    std::mt19937 rng; // Generator liczb pseudolosowych (Mersenne Twister)

    // Funkcja pomocnicza: Sprawdza, czy dwa kubity są podpięte pod ten sam rezonator (gwiazdę)
    bool share_resonator(int q1, int q2) {
        for (const auto& res : resonators) {
            if (res.connected_qubits.count(q1) && res.connected_qubits.count(q2)) {
                return true;
            }
        }
        return false;
    }

public:
    // Konstruktor pobierający dane o topologii
    CircuitGenerator(int qubits_count, const vector<Resonator>& res, int seed = 42)
        : num_qubits(qubits_count), resonators(res), rng(seed) {}

    // Główna funkcja generująca obwód
    std::vector<Gate> generate(int num_gates, double prob_2q, double prob_global) {
        std::vector<Gate> circuit;
        std::uniform_real_distribution<double> dist_prob(0.0, 1.0);
        std::uniform_int_distribution<int> dist_q(0, num_qubits - 1);
        std::uniform_int_distribution<int> dist_res(0, resonators.size() - 1);

        for (int i = 0; i < num_gates; ++i) {
            double p1 = dist_prob(rng);

            if (p1 > prob_2q) {
                // Generuj bramkę 1-kubitową (bez udziału rezonatora, lokalna na kubicie)
                circuit.push_back({i, "1Q", dist_q(rng), -1});
            } else {
                double p2 = dist_prob(rng);

                if (p2 > prob_global) {
                    // BRAMKA 2Q LOKALNA (Wewnątrz jednej gwiazdy/rezonatora - brak SWAPów)
                    int r_id = dist_res(rng);
                    const auto& connected = resonators[r_id].connected_qubits;

                    if (connected.size() >= 2) {
                        std::vector<int> q_list(connected.begin(), connected.end());
                        std::uniform_int_distribution<int> dist_local_q(0, q_list.size() - 1);

                        int q1 = q_list[dist_local_q(rng)];
                        int q2 = q_list[dist_local_q(rng)];
                        while(q1 == q2) q2 = q_list[dist_local_q(rng)]; // Upewnij się, że są różne

                        circuit.push_back({i,"2Q_LOCAL", q1, q2});
                    } else {
                        // Fallback, jeśli rezonator ma tylko 1 kubit (nie powinno wystąpić w zdrowej siatce)
                        circuit.push_back({i,"1Q", dist_q(rng), -1});
                    }
                } else {
                    // BRAMKA 2Q GLOBALNA (Między-gwiazdowa - wymusi router i SWAPy)
                    int q1 = dist_q(rng);
                    int q2 = dist_q(rng);

                    // Szukaj dopóki kubity są takie same ALBO współdzielą ten sam rezonator
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

        std::vector<std::stringstream> lines(num_qubits);

        // Inicjalizacja nagłówków (q0, q1, ..., q10) z wyrównaniem
        for (int i = 0; i < num_qubits; ++i) {
            lines[i] << "q";
            if (i < 10) lines[i] << "0"; // Zero padding dla równych linii
            lines[i] << i << ": ";
        }

        for (const auto& gate : circuit) {
            // Szukamy najdłuższej linii w aktualnej "kolumnie" czasu
            int max_len = 0;
            for (const auto& ss : lines) {
                int len = ss.str().length();
                if (len > max_len) max_len = len;
            }

            // Dopełniamy wszystkie linie znakami '-', aby bramka weszła równo
            for (auto& ss : lines) {
                int pad = max_len - ss.str().length();
                for (int p = 0; p < pad; ++p) ss << "-";
            }

            // Rysowanie poszczególnych bramek
            if (gate.type == "1Q") {
                for (int i = 0; i < num_qubits; ++i) {
                    if (i == gate.target_q1) {
                        lines[i] << "-[ U3 ]-";
                    } else {
                        lines[i] << "--------";
                    }
                }
            }
            else if (gate.type == "2Q_LOCAL" || gate.type == "2Q_GLOBAL") {
                int c = gate.target_q1;
                int t = gate.target_q2;
                int min_q = std::min(c, t);
                int max_q = std::max(c, t);

                // Odroznienie wizualne trudnych bramek
                std::string ctrl_char = (gate.type == "2Q_GLOBAL") ? "--*--" : "--o--";
                std::string targ_char = (gate.type == "2Q_GLOBAL") ? "-[+]-" : "-(+)-";

                for (int i = 0; i < num_qubits; ++i) {
                    if (i == c) lines[i] << ctrl_char;
                    else if (i == t) lines[i] << targ_char;
                    else if (i > min_q && i < max_q) lines[i] << "--|--";
                    else lines[i] << "-----";
                }
            }
        }

        // Wydruk końcowy
        std::cout << "\n=== WIZUALIZACJA OBWODU KWANTOWEGO ===\n";
        std::cout << "Legenda: --o--/-(+)- [Lokalna]  |  --*--/-[+]- [Globalna - wymaga SWAP]\n\n";

        for (int i = 0; i < num_qubits; ++i) {
            std::cout << lines[i].str() << "\n";
        }
        std::cout << "======================================\n\n";
    }

    void link_gates(std::vector<Gate>& circuit) {
        // 1. Resetujemy wszystkie wskaźniki (na wypadek ponownego przeliczania)
        for (auto& g : circuit) {
            g.prev_gate_q1 = g.prev_gate_q2 = nullptr;
            g.next_gate_q1 = g.next_gate_q2 = nullptr;
        }

        // 2. Mapa: ID qubitu -> wskaźnik na ostatnią bramkę, która go użyła
        std::vector<Gate*> last_gate_on_qubit(this->num_qubits, nullptr);

        for (int i = 0; i < circuit.size(); ++i) {
            Gate* current = &circuit[i];

            // --- OBSŁUGA QUBITU 1 (zawsze obecny) ---
            if (last_gate_on_qubit[current->target_q1] != nullptr) {
                Gate* prev = last_gate_on_qubit[current->target_q1];

                // Ustawiamy wstecz w bieżącej
                current->prev_gate_q1 = prev;

                // Ustawiamy w przód w poprzedniej (sprawdzamy na którym qubicie była poprzednia)
                if (prev->target_q1 == current->target_q1) {
                    prev->next_gate_q1 = current;
                } else {
                    prev->next_gate_q2 = current;
                }
            }
            last_gate_on_qubit[current->target_q1] = current;

            // --- OBSŁUGA QUBITU 2 (tylko dla bramek 2Q) ---
            if (current->type != "1Q" && current->target_q2 != -1) {
                if (last_gate_on_qubit[current->target_q2] != nullptr) {
                    Gate* prev = last_gate_on_qubit[current->target_q2];

                    // Ustawiamy wstecz w bieżącej
                    current->prev_gate_q2 = prev;

                    // Ustawiamy w przód w poprzedniej
                    if (prev->target_q1 == current->target_q2) {
                        prev->next_gate_q1 = current;
                    } else {
                        prev->next_gate_q2 = current;
                    }
                }
                last_gate_on_qubit[current->target_q2] = current;
            }
        }
    }
};


#endif //MAGISTERKA_CIRCUITGENERATOR_H