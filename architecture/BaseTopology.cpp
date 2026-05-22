#include "BaseTopology.h"
#include <iostream>

using namespace std;

void BaseTopology::precompute_all_paths() {
    int n = resonators.size();

    // Alokujemy pamięć jednorazowo w ciągłym bloku (n * n elementów)
    distance_matrix.assign(n * n, INF);
    next_node_matrix.assign(n * n, -1);

    // 1. Prekalkulacja ścieżek dla rezonatorów (Algorytm Floyda-Warshalla na płaskiej macierzy)
    for (int i = 0; i < n; ++i) {
        distance_matrix[i * n + i] = 0;
        next_node_matrix[i * n + i] = i;

        for (int j = 0; j < n; ++j) {
            if (i != j && !get_shared_qubits(i, j).empty()) {
                distance_matrix[i * n + j] = 1;
                next_node_matrix[i * n + j] = j;
            }
        }
    }

    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (distance_matrix[i * n + k] < INF && distance_matrix[k * n + j] < INF) {
                    if (distance_matrix[i * n + k] + distance_matrix[k * n + j] < distance_matrix[i * n + j]) {
                        distance_matrix[i * n + j] = distance_matrix[i * n + k] + distance_matrix[k * n + j];
                        next_node_matrix[i * n + j] = next_node_matrix[i * n + k];
                    }
                }
            }
        }
    }

    // 2. Prekalkulacja tras dla kubitów z użyciem płaskiej macierzy
    int num_qubits = qubits.size();
    routing_paths_matrix.assign(num_qubits * num_qubits, RoutingPath());

    for (int q1 = 0; q1 < num_qubits; ++q1) {
        for (int q2 = 0; q2 < num_qubits; ++q2) {
            routing_paths_matrix[q1 * num_qubits + q2] = calculate_routing_path(q1, q2);
        }
    }
}

vector<int> BaseTopology::get_shared_qubits(int res1_id, int res2_id) const {
    vector<int> shared;
    const auto& res2_qubits = resonators[res2_id].connected_qubits;
    for (int q : resonators[res1_id].connected_qubits) {
        if (ranges::find(res2_qubits, q) != res2_qubits.end()) {
            shared.push_back(q);
        }
    }
    return shared;
}

int BaseTopology::get_swap_cost(int start_res, int target_res) const {
    int n = resonators.size();
    return distance_matrix[start_res * n + target_res];
}

vector<int> BaseTopology::get_resonator_path(int start_res, int target_res) const {
    int n = resonators.size();
    if (distance_matrix[start_res * n + target_res] == INF) return {};

    vector<int> path;
    int current = start_res;
    path.push_back(current);

    while (current != target_res) {
        current = next_node_matrix[current * n + target_res];
        path.push_back(current);
    }
    return path;
}

void BaseTopology::print_stats() const {
    cout << "--- Statystyki QPU ---" << endl;
    cout << "Liczba Rezonatorow: " << resonators.size() << endl;
    cout << "Calkowita Liczba Kubitow: " << qubits.size() << endl;

    vector<int> shared_by_3, shared_by_2, shared_by_1;

    for (const auto& q : qubits) {
        if (q.connected_resonators.size() == 3) shared_by_3.push_back(q.id);
        else if (q.connected_resonators.size() == 2) shared_by_2.push_back(q.id);
        else if (q.connected_resonators.size() == 1) shared_by_1.push_back(q.id);
    }

    cout << "\nPodzial Kubitow ze wzgledu na ilosc polaczen:" << endl;
    cout << " - Polaczone z 3 rezonatorami (Centralne): " << shared_by_3.size() << " (ID: ";
    for(int id : shared_by_3) cout << id << " "; cout << ")" << endl;

    cout << " - Polaczone z 2 rezonatorami (Wewnetrzne Krawedzie): " << shared_by_2.size() << endl;
    cout << " - Polaczone z 1 rezonatorem (Brzegowe): " << shared_by_1.size() << endl;
}

RoutingPath BaseTopology::calculate_routing_path(int phys_q1, int phys_q2) const {
    RoutingPath path;
    int n = resonators.size();

    int min_dist = INF;
    int best_r1 = -1, best_r2 = -1;

    for (int r1 : qubits[phys_q1].connected_resonators) {
        for (int r2 : qubits[phys_q2].connected_resonators) {
            // Używamy odczytu 1D: [wiersz * szerokość + kolumna]
            if (distance_matrix[r1 * n + r2] < min_dist) {
                min_dist = distance_matrix[r1 * n + r2];
                best_r1 = r1;
                best_r2 = r2;
            }
        }
    }

    if (best_r1 == -1 || best_r2 == -1) return path;

    std::vector<int> res_path = get_resonator_path(best_r1, best_r2);

    path.qubits.push_back(phys_q1);

    for (size_t i = 0; i < res_path.size() - 1; ++i) {
        int r_curr = res_path[i];
        int r_next = res_path[i+1];

        path.resonators.push_back(r_curr);

        std::vector<int> shared = get_shared_qubits(r_curr, r_next);

        if (!shared.empty()) {
            path.qubits.push_back(shared.front());
        }
    }

    path.resonators.push_back(res_path.back());
    path.qubits.push_back(phys_q2);

    return path;
}

// Błyskawiczne pobranie gotowej ścieżki w czasie O(1)
const RoutingPath& BaseTopology::get_routing_path(int phys_q1, int phys_q2) const {
    int num_qubits = qubits.size();
    return routing_paths_matrix[phys_q1 * num_qubits + phys_q2];
}