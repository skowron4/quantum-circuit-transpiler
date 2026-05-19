//
// Created by skowr on 18.05.2026.
//

#ifndef MAGISTERKA_BASETOPOLOGY_H
#define MAGISTERKA_BASETOPOLOGY_H


#pragma once
#include <vector>
#include <map>
#include <set>

#include "structures/Qubit.h"
#include "structures/Resonator.h"

struct RoutingPath {
    std::vector<int> qubits;
    std::vector<int> resonators;
};

class BaseTopology {
public:
    std::vector<Resonator> resonators;
    std::vector<Qubit> qubits;
    const int INF = 1e9;

protected:
    std::map<std::set<int>, int> qubit_map;
    std::vector<std::vector<int>> distance_matrix;
    std::vector<std::vector<int>> next_node_matrix;

    // Czysto wirtualna funkcja - każda topologia buduje się inaczej
    virtual void build_topology() = 0;

    // Wspólny algorytm Floyda-Warshalla
    void precompute_all_paths();

public:
    virtual ~BaseTopology() = default;

    std::vector<int> get_shared_qubits(int res1_id, int res2_id) const;
    int get_swap_cost(int start_res, int target_res) const;
    std::vector<int> get_resonator_path(int start_res, int target_res) const;
    RoutingPath get_routing_path(int phys_q1, int phys_q2) const;
    void print_stats() const;

    int get_num_qubits() const { return qubits.size(); }
    const std::vector<Resonator>& get_resonators() const { return resonators; }
};

#endif //MAGISTERKA_BASETOPOLOGY_H