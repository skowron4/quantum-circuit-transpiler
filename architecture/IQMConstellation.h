#pragma once
#include "BaseTopology.h"

class IQMConstellation : public BaseTopology {
private:
    int Y_rows, X_cols;

    int get_res_id(int y, int x) const {
        if (y < 0 || y >= Y_rows || x < 0 || x >= X_cols) return -1;
        return y * X_cols + x;
    }

protected:
    void build_topology() override {
        // 1. Inicjalizacja rezonatorów
        for (int y = 0; y < Y_rows; ++y) {
            for (int x = 0; x < X_cols; ++x) {
                resonators.push_back({y * X_cols + x, y, x, {}});
            }
        }

        int dy[] = {0, -1, -1, 0, 1, 1};
        int dx[] = {1, 1, 0, -1, -1, 0};

        // 2. Budowa połączeń
        for (int y = 0; y < Y_rows; ++y) {
            for (int x = 0; x < X_cols; ++x) {
                int r_id = get_res_id(y, x);

                for (int i = 0; i < 6; ++i) {
                    int n1_id = get_res_id(y + dy[i], x + dx[i]);
                    int next_i = (i + 1) % 6;
                    int n2_id = get_res_id(y + dy[next_i], x + dx[next_i]);
                    std::set<int> res_group = {r_id};
                    if (n1_id != -1) res_group.insert(n1_id);
                    if (n2_id != -1) res_group.insert(n2_id);

                    if (qubit_map.find(res_group) == qubit_map.end()) {
                        int new_q_id = static_cast<int>(qubits.size());
                        // Konwersja set -> InlinedVector
                        absl::InlinedVector<int, 5> vec_group(res_group.begin(), res_group.end());
                        qubits.push_back({new_q_id, vec_group});
                        qubit_map[res_group] = new_q_id;
                    }

                    int q_id = qubit_map[res_group];

                    auto& connected = resonators[r_id].connected_qubits;
                    if (std::ranges::find(connected, q_id) == connected.end()) {
                        connected.push_back(q_id);
                    }
                }
            }
        }

        for (auto& res : resonators) {
            while (res.connected_qubits.size() < 6) {
                int new_q_id = static_cast<int>(qubits.size());
                absl::InlinedVector<int, 5> res_group = {res.id};
                qubits.push_back({new_q_id, res_group});
                res.connected_qubits.push_back(new_q_id);
            }
        }
    }

public:
    IQMConstellation(int y, int x) : Y_rows(y), X_cols(x) {
        build_topology();
        precompute_all_paths();
    }
};