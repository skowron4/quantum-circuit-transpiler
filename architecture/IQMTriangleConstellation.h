//
// Created by skowr on 6.05.2026.
//

#ifndef MAGISTERKA_IQMTRIANGLECONSTELLATION_H
#define MAGISTERKA_IQMTRIANGLECONSTELLATION_H

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include "../structures/Qubit.h"
#include "../structures/Resonator.h"

using namespace std;

class IQMTriangleConstellation {
private:
    vector<Resonator> resonators;
    vector<Qubit> qubits;

    map<set<int>, int> qubit_map;
    map<pair<int, int>, int> coord_to_id; // Mapa wspolrzednych do ID rezonatora

    // Macierze do prekomputacji (O(1) dla funkcji oceny)
    vector<vector<int>> distance_matrix;
    vector<vector<int>> next_node_matrix;
    const int INF = 1e9;

    // Helper dodajacy rezonator
    void add_resonator(int id, int y, int x) {
        resonators.push_back({id, y, x, {}});
        coord_to_id[{y, x}] = id;
    }

    // Helper pobierajacy ID na podstawie wspolrzednych
    int get_res_id(int y, int x) const {
        auto it = coord_to_id.find({y, x});
        if (it != coord_to_id.end()) return it->second;
        return -1;
    }

    void precompute_all_paths() {
        int n = resonators.size();
        distance_matrix.assign(n, vector<int>(n, INF));
        next_node_matrix.assign(n, vector<int>(n, -1));

        for (int i = 0; i < n; ++i) {
            distance_matrix[i][i] = 0;
            next_node_matrix[i][i] = i;

            for (int j = 0; j < n; ++j) {
                if (i != j && !get_shared_qubits(i, j).empty()) {
                    distance_matrix[i][j] = 1;
                    next_node_matrix[i][j] = j;
                }
            }
        }

        // Algorytm Floyda-Warshalla
        for (int k = 0; k < n; ++k) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (distance_matrix[i][k] < INF && distance_matrix[k][j] < INF) {
                        if (distance_matrix[i][k] + distance_matrix[k][j] < distance_matrix[i][j]) {
                            distance_matrix[i][j] = distance_matrix[i][k] + distance_matrix[k][j];
                            next_node_matrix[i][j] = next_node_matrix[i][k];
                        }
                    }
                }
            }
        }
    }

public:
    IQMTriangleConstellation() {
        build_topology();
        precompute_all_paths();
    }

    void build_topology() {
        // Ręczne definiowanie współrzędnych, które w siatce heksagonalnej tworzą ścisły trójkąt
        add_resonator(0, 0, 0);   // Rezonator 0
        add_resonator(1, 0, 1);   // Rezonator 1
        add_resonator(2, -1, 1);  // Rezonator 2

        // Wektory sąsiedztwa w siatce trójkątnej (6 kierunków)
        int dy[] = {0, -1, -1, 0, 1, 1};
        int dx[] = {1, 1, 0, -1, -1, 0};

        // Generowanie Kubitów (Wyszukiwanie części wspólnych)
        for (auto& res : resonators) {
            for (int i = 0; i < 6; ++i) {
                int n1_id = get_res_id(res.y + dy[i], res.x + dx[i]);
                int next_i = (i + 1) % 6;
                int n2_id = get_res_id(res.y + dy[next_i], res.x + dx[next_i]);

                set<int> res_group;
                res_group.insert(res.id);
                if (n1_id != -1) res_group.insert(n1_id);
                if (n2_id != -1) res_group.insert(n2_id);

                if (qubit_map.find(res_group) == qubit_map.end()) {
                    int new_q_id = qubits.size();
                    qubits.push_back({new_q_id, res_group});
                    qubit_map[res_group] = new_q_id;
                }

                // Przypinamy kubit do rezonatora
                res.connected_qubits.insert(qubit_map[res_group]);
            }
        }

        // --- DODANY KOD: Uzupełnianie do 6 kubitów (brzegowe) ---
        // Po wygenerowaniu współdzielonych połączeń, sprawdzamy każdy rezonator
        // i dodajemy unikalne kubity brzegowe (należące tylko do 1 rezonatora)
        // tak długo, aż dany rezonator będzie miał dokładnie 6 obsługiwanych kubitów.
        for (auto& res : resonators) {
            while (res.connected_qubits.size() < 6) {
                int new_q_id = qubits.size();
                set<int> res_group = {res.id}; // Kubit połączony tylko z tym konkretnym rezonatorem

                qubits.push_back({new_q_id, res_group});
                res.connected_qubits.insert(new_q_id);
            }
        }
    }

    vector<int> get_shared_qubits(int res1_id, int res2_id) const {
        vector<int> shared;
        for (int q : resonators[res1_id].connected_qubits) {
            if (resonators[res2_id].connected_qubits.count(q)) {
                shared.push_back(q);
            }
        }
        return shared;
    }

    int get_swap_cost(int start_res, int target_res) const {
        return distance_matrix[start_res][target_res];
    }

    vector<int> get_resonator_path(int start_res, int target_res) const {
        if (distance_matrix[start_res][target_res] == INF) return {};
        vector<int> path;
        int current = start_res;
        path.push_back(current);
        while (current != target_res) {
            current = next_node_matrix[current][target_res];
            path.push_back(current);
        }
        return path;
    }

    void print_stats() const {
        cout << "--- Statystyki QPU (Toy Model 3-Rezonatory) ---" << endl;
        cout << "Liczba Rezonatorów: " << resonators.size() << endl;
        cout << "Calkowita Liczba Kubitów: " << qubits.size() << endl;

        vector<int> shared_by_3, shared_by_2, shared_by_1;

        for (const auto& q : qubits) {
            if (q.connected_resonators.size() == 3) shared_by_3.push_back(q.id);
            else if (q.connected_resonators.size() == 2) shared_by_2.push_back(q.id);
            else if (q.connected_resonators.size() == 1) shared_by_1.push_back(q.id);
        }

        cout << "\nPodzial Kubitów ze wzgledu na ilosc polaczen:" << endl;
        cout << " - Polaczone z 3 rezonatorami (Centralne): " << shared_by_3.size() << " (ID: ";
        for(int id : shared_by_3) cout << id << " "; cout << ")" << endl;

        cout << " - Polaczone z 2 rezonatorami (Wewnetrzne Krawedzie): " << shared_by_2.size() << endl;
        cout << " - Polaczone z 1 rezonatorem (Brzegowe): " << shared_by_1.size() << endl;
    }

    int get_num_qubits() const { return qubits.size(); }
    const vector<Resonator>& get_resonators() const { return resonators; }
};

#endif //MAGISTERKA_IQMTRIANGLECONSTELLATION_H