#pragma once
#include "BaseTopology.h"
#include <map>
#include <utility>

class IQMTriangleConstellation : public BaseTopology {
private:
    std::map<std::pair<int, int>, int> coord_to_id;

    void add_resonator(int id, int y, int x) {
        resonators.push_back({id, y, x, {}});
        coord_to_id[{y, x}] = id;
    }

    int get_res_id(int y, int x) const {
        auto it = coord_to_id.find({y, x});
        if (it != coord_to_id.end()) return it->second;
        return -1;
    }

protected:
    void build_topology() override {
        add_resonator(0, 0, 0);
        add_resonator(1, 0, 1);
        add_resonator(2, -1, 1);

        int dy[] = {0, -1, -1, 0, 1, 1};
        int dx[] = {1, 1, 0, -1, -1, 0};

        for (auto& res : resonators) {
            for (int i = 0; i < 6; ++i) {
                int n1_id = get_res_id(res.y + dy[i], res.x + dx[i]);
                int next_i = (i + 1) % 6;
                int n2_id = get_res_id(res.y + dy[next_i], res.x + dx[next_i]);

                std::set<int> res_group_set;
                res_group_set.insert(res.id);
                if (n1_id != -1) res_group_set.insert(n1_id);
                if (n2_id != -1) res_group_set.insert(n2_id);

                if (qubit_map.find(res_group_set) == qubit_map.end()) {
                    int new_q_id = qubits.size();
                    absl::InlinedVector<int, 5> res_group_vec(res_group_set.begin(), res_group_set.end());
                    qubits.push_back({new_q_id, res_group_vec});
                    qubit_map[res_group_set] = new_q_id;
                }

                int mapped_q_id = qubit_map[res_group_set];
                if (std::ranges::find(res.connected_qubits, mapped_q_id) == res.connected_qubits.end()) {
                    res.connected_qubits.push_back(mapped_q_id);
                }
            }
        }

        for (auto& res : resonators) {
            while (res.connected_qubits.size() < 6) {
                int new_q_id = qubits.size();
                absl::InlinedVector<int, 5> res_group = {res.id};
                qubits.push_back({new_q_id, res_group});
                res.connected_qubits.push_back(new_q_id);
            }
        }
    }

public:
    IQMTriangleConstellation() {
        build_topology();
        precompute_all_paths();
    }
};