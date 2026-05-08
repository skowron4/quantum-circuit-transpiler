//
// Created by skowr on 6.05.2026.
//

#ifndef MAGISTERKA_IQMCONSTELLATION_H
#define MAGISTERKA_IQMCONSTELLATION_H


#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <algorithm>

using namespace std;

// Struktura reprezentująca Kubit (połączony z 1, 2 lub 3 rezonatorami)
struct Qubit {
    int id;
    set<int> connected_resonators;
};

// Struktura reprezentująca Rezonator / Gwiazdę
struct Resonator {
    int id;
    int y, x;
    set<int> connected_qubits;
};

class IQMConstellation {
private:
    int Y_rows, X_cols;
    vector<Resonator> resonators;
    vector<Qubit> qubits;

    map<set<int>, int> qubit_map;

    // --- NOWE: Macierze do prekomputacji ---
    vector<vector<int>> distance_matrix;
    vector<vector<int>> next_node_matrix;
    const int INF = 1e9; // Nieskończoność dla początkowego stanu macierzy

    int get_res_id(int y, int x) const {
        if (y < 0 || y >= Y_rows || x < 0 || x >= X_cols) return -1;
        return y * X_cols + x;
    }

    // Algorytm Floyda-Warshalla wyliczający najkrótsze ścieżki i następne kroki
    void precompute_all_paths() {
        int n = resonators.size();
        distance_matrix.assign(n, vector<int>(n, INF));
        next_node_matrix.assign(n, vector<int>(n, -1));

        // 1. Inicjalizacja macierzy na podstawie bezpośrednich sąsiadów
        for (int i = 0; i < n; ++i) {
            distance_matrix[i][i] = 0;
            next_node_matrix[i][i] = i;

            for (int j = 0; j < n; ++j) {
                if (i != j && !get_shared_qubits(i, j).empty()) {
                    distance_matrix[i][j] = 1;      // 1 skok (SWAP)
                    next_node_matrix[i][j] = j;     // Aby dojść z 'i' do 'j', idź do 'j'
                }
            }
        }

        // 2. Właściwa pętla Floyda-Warshalla O(V^3)
        for (int k = 0; k < n; ++k) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (distance_matrix[i][k] < INF && distance_matrix[k][j] < INF) {
                        if (distance_matrix[i][k] + distance_matrix[k][j] < distance_matrix[i][j]) {
                            distance_matrix[i][j] = distance_matrix[i][k] + distance_matrix[k][j];
                            // Zapisujemy dokąd iść z 'i', żeby optymalnie dotrzeć do 'j' przez 'k'
                            next_node_matrix[i][j] = next_node_matrix[i][k];
                        }
                    }
                }
            }
        }
    }

public:
    IQMConstellation(int y, int x) : Y_rows(y), X_cols(x) {
        build_topology();
        precompute_all_paths(); // Wywoływane RAZ przy starcie
    }

    void build_topology() {
        for (int y = 0; y < Y_rows; ++y) {
            for (int x = 0; x < X_cols; ++x) {
                resonators.push_back({y * X_cols + x, y, x, {}});
            }
        }

        int dy[] = {0, -1, -1, 0, 1, 1};
        int dx[] = {1, 1, 0, -1, -1, 0};

        for (int y = 0; y < Y_rows; ++y) {
            for (int x = 0; x < X_cols; ++x) {
                int r_id = get_res_id(y, x);

                for (int i = 0; i < 6; ++i) {
                    int n1_id = get_res_id(y + dy[i], x + dx[i]);
                    int next_i = (i + 1) % 6;
                    int n2_id = get_res_id(y + dy[next_i], x + dx[next_i]);

                    set<int> res_group;
                    res_group.insert(r_id);
                    if (n1_id != -1) res_group.insert(n1_id);
                    if (n2_id != -1) res_group.insert(n2_id);

                    if (qubit_map.find(res_group) == qubit_map.end()) {
                        int new_q_id = qubits.size();
                        qubits.push_back({new_q_id, res_group});
                        qubit_map[res_group] = new_q_id;
                    }
                    resonators[r_id].connected_qubits.insert(qubit_map[res_group]);
                }
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

    // --- NOWE: O(1) funkcja oceny ---
    int get_swap_cost(int start_res, int target_res) const {
        return distance_matrix[start_res][target_res];
    }

    // --- NOWE: O(L) Błyskawiczne odtwarzanie trasy ---
    vector<int> get_resonator_path(int start_res, int target_res) const {
        if (distance_matrix[start_res][target_res] == INF) return {}; // Brak ścieżki

        vector<int> path;
        int current = start_res;
        path.push_back(current);

        // Dopóki nie dotrzemy do celu, odczytujemy następny węzeł z macierzy
        while (current != target_res) {
            current = next_node_matrix[current][target_res];
            path.push_back(current);
        }

        return path;
    }

    void print_stats() const {
        cout << "--- Statystyki QPU ---" << endl;
        cout << "Liczba Rezonatorów: " << resonators.size() << endl;
        cout << "Liczba Kubitów: " << qubits.size() << endl;
        cout << "Macierze odleglosci i routingu zainicjalizowane." << endl;
    }

    int get_num_qubits() const { return qubits.size(); }
    const vector<Resonator>& get_resonators() const { return resonators; }
};


#endif //MAGISTERKA_IQMCONSTELLATION_H