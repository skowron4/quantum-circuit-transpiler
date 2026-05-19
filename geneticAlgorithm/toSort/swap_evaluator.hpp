#ifndef SWAP_EVALUATOR_HPP
#define SWAP_EVALUATOR_HPP

#include <cstdint>
#include <numeric>
#include <vector>
#include <algorithm>
#include <map>
#include <queue>
#include <unordered_map>
#include <unordered_set>

#include "../templates/problem.hpp"
#include "../rng.hpp"
#include "structures/Qubit.h"
#include "structures/Resonator.h"
#include "structures/Gate.h"
#include "structures/physical_qubit.hpp"
#include "structures/routing_record.hpp"


namespace metah {
    class SwapEvaluator final : public metah::Problem<std::vector<std::pair<int, int>>, uint64_t> {
    public:
        std::vector<Gate> gates;
        std::vector<Resonator> resonators;

        std::unordered_map<int, Gate> gate_map;
        std::unordered_map<int, PhysicalQubit> locations; // Wszystkie unikalne lokacje fizyczne (na bazie kombinacji rezonatorów)
        std::unordered_map<int, std::vector<int>> loc_adj; // Graf sąsiedztwa lokacji fizycznych

        // Cache ścieżek: klucz to spakowane (start_loc_id, target_res_id, avoid_loc_id)
        std::unordered_map<uint64_t, std::vector<int>> path_cache;

        // Początkowe mapowanie kubitów logicznych na ID lokacji fizycznych
        std::unordered_map<int, int> initial_lq_to_loc;

        [[nodiscard]] static uint64_t make_cache_key(int start_loc, int target_res, int avoid_loc) {
            return (static_cast<uint64_t>(start_loc) << 42) |
                   (static_cast<uint64_t>(target_res) << 21) |
                   static_cast<uint64_t>(avoid_loc);
        }

    public:
        explicit SwapEvaluator(const std::vector<Gate>& gates_in,
                               const std::vector<Resonator>& resonators_in,
                               const std::vector<Qubit>& qubits_in)
            : gates(gates_in), resonators(resonators_in)
        {
            for (const auto& g : gates) gate_map[g.id] = g;

            // 1. Inicjalizacja fizycznych lokalizacji na bazie obiektów Qubit
            for (const auto& lq : qubits_in) {
                int loc_id = lq.id; // Zakładamy, że ID startowego kubita to ID jego gniazda fizycznego
                locations[loc_id] = {loc_id, lq.connected_resonators};
                initial_lq_to_loc[lq.id] = loc_id;
            }

            // 2. Budowa grafu sąsiedztwa lokacji fizycznych
            for (const auto& [loc1_id, loc1] : locations) {
                for (const auto& [loc2_id, loc2] : locations) {
                    if (loc1_id == loc2_id) continue;

                    // Lokacje są sąsiadami, jeśli mają chociaż jeden wspólny rezonator
                    bool share_resonator = false;
                    for (int res_id : loc1.connected_resonators) {
                        if (loc2.connected_resonators.contains(res_id)) {
                            share_resonator = true;
                            break;
                        }
                    }
                    if (share_resonator) {
                        loc_adj[loc1_id].push_back(loc2_id);
                    }
                }
            }

            // 3. PREKALKULACJA ŚCIEŻEK (tylko dla istniejących, prawidłowych lokacji fizycznych!)
            for (const auto& [start_loc_id, _] : locations) {
                for (const auto& res : resonators) {
                    for (const auto& [avoid_loc_id, __] : locations) {
                        if (start_loc_id == avoid_loc_id) continue;

                        uint64_t key = make_cache_key(start_loc_id, res.id, avoid_loc_id);
                        path_cache[key] = compute_shortest_path(start_loc_id, res.id, avoid_loc_id);
                    }
                }
            }
        }

        [[nodiscard]] uint64_t evaluate(const std::vector<std::pair<int, int>>& repr) const override {
            uint64_t total_swaps = 0;

            // Kopiujemy stan początkowy mapowania dla każdego wywołania evaluate
            std::unordered_map<int, int> lq_to_loc = initial_lq_to_loc;

            // Odwrotne mapowanie: loc_id -> lq_id (potrzebne do swapów)
            std::unordered_map<int, int> loc_to_lq;
            for (const auto& [lq, loc] : lq_to_loc) {
                loc_to_lq[loc] = lq;
            }

            auto apply_swaps = [&](const std::vector<int>& path, uint64_t& total_cost) {
                if (path.empty() || path.size() == 1) return;

                if (path[0] == -1) {
                    total_cost += 10000; // Kara za brak możliwości rutingu
                    printf("KURWAAAAAAAAAAAAAAA");
                    return;
                }

                // Przechodzimy po fizycznej ścieżce i zamieniamy miejscami kubity logiczne
                for (size_t i = 0; i < path.size() - 1; ++i) {
                    int loc1 = path[i];
                    int loc2 = path[i+1];

                    int lq1 = loc_to_lq[loc1];
                    int lq2 = loc_to_lq[loc2];

                    // Zamiana mapowań
                    lq_to_loc[lq1] = loc2;
                    lq_to_loc[lq2] = loc1;
                    loc_to_lq[loc1] = lq2;
                    loc_to_lq[loc2] = lq1;

                    total_cost++;
                }
            };

            for (const auto& [gate_id, res_id] : repr) {
                const auto& gate = gate_map.at(gate_id);

                if (gate.type == "2Q_LOCAL" || gate.type == "2Q_GLOBAL") {
                    int lq1 = gate.target_q1;
                    int lq2 = gate.target_q2;

                    // Gdzie fizycznie znajdują się teraz te kubity?
                    int loc1 = lq_to_loc[lq1];
                    int loc2 = lq_to_loc[lq2];

                    // Ścieżka dla pierwszego kubita
                    uint64_t key1 = make_cache_key(loc1, res_id, loc2);
                    const auto& path1 = path_cache.at(key1);
                    apply_swaps(path1, total_swaps);

                    // Zaktualizowana lokalizacja pierwszego kubita (bo mógł wykonać swapy)
                    loc1 = lq_to_loc[lq1];

                    // Ścieżka dla drugiego kubita
                    uint64_t key2 = make_cache_key(loc2, res_id, loc1);
                    const auto& path2 = path_cache.at(key2);
                    apply_swaps(path2, total_swaps);
                }
            }

            return total_swaps;
        }

        [[nodiscard]] std::vector<std::pair<int, int>> random_repr(metah::ThreadLocalRandom& rng) const override {
            std::vector<std::pair<int, int>> repr;
            repr.reserve(gates.size());
            std::unordered_set<int> completed_gates;
            std::vector<const Gate*> available_gates;

            for (const auto& gate : gates) {
                if (gate.prev_gate_q1 == nullptr && gate.prev_gate_q2 == nullptr) {
                    available_gates.push_back(&gate);
                }
            }

            auto is_ready = [&](const Gate* g) {
                if (!g) return false;
                if (g->prev_gate_q1 && completed_gates.find(g->prev_gate_q1->id) == completed_gates.end()) return false;
                if (g->prev_gate_q2 && completed_gates.find(g->prev_gate_q2->id) == completed_gates.end()) return false;
                return true;
            };

            while (!available_gates.empty()) {
                size_t idx = rng.random_idx(available_gates.size() - 1);
                const Gate* current = available_gates[idx];

                std::swap(available_gates[idx], available_gates.back());
                available_gates.pop_back();

                repr.emplace_back(current->id, rng.random_idx(resonators.size() - 1));
                completed_gates.insert(current->id);

                if (current->next_gate_q1 && is_ready(current->next_gate_q1)) {
                    available_gates.push_back(current->next_gate_q1);
                }

                if (current->next_gate_q2 && current->next_gate_q2 != current->next_gate_q1) {
                    if (is_ready(current->next_gate_q2)) {
                        available_gates.push_back(current->next_gate_q2);
                    }
                }
            }

            return repr;
        }

        //TODO
        [[nodiscard]] std::vector<RoutingRecord> evaluate_and_track_swaps(const std::vector<std::pair<int, int>>& repr) const /* override - jeśli tworzysz nową wirtualną */ {
            std::vector<RoutingRecord> swap_history;

            // Kopiujemy stan początkowy mapowania dla każdego wywołania
            std::unordered_map<int, int> lq_to_loc = initial_lq_to_loc;

            // Odwrotne mapowanie: loc_id -> lq_id (potrzebne do swapów)
            std::unordered_map<int, int> loc_to_lq;
            for (const auto& [lq, loc] : lq_to_loc) {
                loc_to_lq[loc] = lq;
            }

            // Zmodyfikowana lambda - teraz przyjmuje ID aktualnej bramki zamiast referencji do kosztu
            auto apply_swaps = [&](const std::vector<int>& path, int current_gate_id) {
                if (path.empty() || path.size() == 1) return;

                if (path[0] == -1) {
                    // Tutaj w zależności od potrzeb możesz wrzucić do wektora jakiś błąd,
                    // rzucić wyjątek lub po prostu zignorować.
                    printf("KURWAAAAAAAAAAAAAAA\n");
                    return;
                }

                // Przechodzimy po fizycznej ścieżce i zamieniamy miejscami kubity logiczne
                for (size_t i = 0; i < path.size() - 1; ++i) {
                    int loc1 = path[i];
                    int loc2 = path[i+1];

                    int lq1 = loc_to_lq[loc1];
                    int lq2 = loc_to_lq[loc2];

                    // Rejestrujemy wykonany swap
                    // swap_history.push_back({current_gate_id, lq1, lq2, loc1, loc2});

                    // Zamiana mapowań
                    lq_to_loc[lq1] = loc2;
                    lq_to_loc[lq2] = loc1;
                    loc_to_lq[loc1] = lq2;
                    loc_to_lq[loc2] = lq1;
                }
            };

            for (const auto& [gate_id, res_id] : repr) {
                const auto& gate = gate_map.at(gate_id);

                if (gate.type == "2Q_LOCAL" || gate.type == "2Q_GLOBAL") {
                    int lq1 = gate.target_q1;
                    int lq2 = gate.target_q2;

                    // Gdzie fizycznie znajdują się teraz te kubity?
                    int loc1 = lq_to_loc[lq1];
                    int loc2 = lq_to_loc[lq2];

                    // Ścieżka dla pierwszego kubita
                    uint64_t key1 = make_cache_key(loc1, res_id, loc2);
                    const auto& path1 = path_cache.at(key1);
                    apply_swaps(path1, gate_id);

                    // Zaktualizowana lokalizacja pierwszego kubita (bo mógł wykonać swapy)
                    loc1 = lq_to_loc[lq1];

                    // Ścieżka dla drugiego kubita
                    uint64_t key2 = make_cache_key(loc2, res_id, loc1);
                    const auto& path2 = path_cache.at(key2);
                    apply_swaps(path2, gate_id);
                }
            }

            return swap_history;
        }

    private:
        [[nodiscard]] std::vector<int> compute_shortest_path(int start_loc_id, int target_res_id, int avoid_loc_id) const {
            // Zamiast pytać rezonator, pytamy samą lokalizację, czy "dotyka" tego rezonatora!
            if (locations.at(start_loc_id).connected_resonators.contains(target_res_id)) {
                return {start_loc_id};
            }

            std::queue<int> q;
            std::unordered_map<int, int> came_from;
            q.push(start_loc_id);
            came_from[start_loc_id] = start_loc_id;

            int found_target_loc = -1;

            while (!q.empty()) {
                int curr = q.front(); q.pop();

                if (locations.at(curr).connected_resonators.contains(target_res_id)) {
                    found_target_loc = curr;
                    break;
                }

                for (int neighbor : loc_adj.at(curr)) {
                    if (neighbor == avoid_loc_id) continue;

                    if (!came_from.contains(neighbor)) {
                        came_from[neighbor] = curr;
                        q.push(neighbor);
                    }
                }
            }

            if (found_target_loc == -1) {
                return {-1};
            }

            std::vector<int> path;
            int curr = found_target_loc;
            while (curr != start_loc_id) {
                path.push_back(curr);
                curr = came_from[curr];
            }
            path.push_back(start_loc_id);
            std::ranges::reverse(path);

            return path;
        }
    };
}


#endif //SWAP_EVALUATOR_HPP
