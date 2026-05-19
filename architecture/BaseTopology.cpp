#include "BaseTopology.h"

#include <iostream>

using namespace std;

void BaseTopology::precompute_all_paths() {
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

vector<int> BaseTopology::get_shared_qubits(int res1_id, int res2_id) const {
    vector<int> shared;
    for (int q : resonators[res1_id].connected_qubits) {
        if (resonators[res2_id].connected_qubits.count(q)) {
            shared.push_back(q);
        }
    }
    return shared;
}

int BaseTopology::get_swap_cost(int start_res, int target_res) const {
    return distance_matrix[start_res][target_res];
}

vector<int> BaseTopology::get_resonator_path(int start_res, int target_res) const {
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

RoutingPath BaseTopology::get_routing_path(int phys_q1, int phys_q2) const {
    RoutingPath path;

    // 1. Szukamy najbliższej pary rezonatorów dla podanych fizycznych kubitów
    int min_dist = INF;
    int best_r1 = -1, best_r2 = -1;

    for (int r1 : qubits[phys_q1].connected_resonators) {
        for (int r2 : qubits[phys_q2].connected_resonators) {
            if (distance_matrix[r1][r2] < min_dist) {
                min_dist = distance_matrix[r1][r2];
                best_r1 = r1;
                best_r2 = r2;
            }
        }
    }

    // Zabezpieczenie na wypadek braku ścieżki (w spójnej topologii nie powinno wystąpić)
    if (best_r1 == -1 || best_r2 == -1) return path;

    // 2. Pobieramy optymalną sekwencję samych rezonatorów
    std::vector<int> res_path = get_resonator_path(best_r1, best_r2);

    // 3. Przeplatamy ścieżkę: Kubit -> Rezonator -> Kubit ...
    path.qubits.push_back(phys_q1);

    for (size_t i = 0; i < res_path.size() - 1; ++i) {
        int r_curr = res_path[i];
        int r_next = res_path[i+1];

        path.resonators.push_back(r_curr);

        // Pobieramy listę kubitów współdzielonych między obecnym a następnym rezonatorem
        std::vector<int> shared = get_shared_qubits(r_curr, r_next);

        // Zawsze bierzemy pierwszy wspólny kubit jako punkt przesiadkowy
        if (!shared.empty()) {
            path.qubits.push_back(shared.front());
        }
    }

    // Dodajemy ostatni rezonator na trasie i kubit docelowy
    path.resonators.push_back(res_path.back());
    path.qubits.push_back(phys_q2);

    return path;
}