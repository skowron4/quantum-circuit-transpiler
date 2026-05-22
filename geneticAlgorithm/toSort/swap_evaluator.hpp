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
#include <absl/container/flat_hash_map.h>
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

        // Zastąpienie std::unordered_map przez szybkie, płaskie tablice z bezpośrednim dostępem (O(1))
        std::vector<Gate> fast_gate_map;
        std::vector<PhysicalQubit> fast_locations;
        std::vector<absl::InlinedVector<int, 8>> fast_loc_adj;

        // Ekstremalnie szybki cache: Oparty na flat_hash_map od Google + InlinedVector (ZERO alokacji na ścieżkach)
        absl::flat_hash_map<uint64_t, absl::InlinedVector<int, 16>> path_cache;

        // Pre-alokowane tablice mapowania, co omija alokacje std::unordered_map w gorącej pętli
        std::vector<int> initial_lq_to_loc;
        std::vector<int> initial_loc_to_lq;

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
            // Odszukanie najwyższych indeksów do zaalokowania płaskich wektorów
            int max_gate_id = 0;
            for (const auto& g : gates) max_gate_id = std::max(max_gate_id, g.id);

            int max_loc_id = 0;
            int max_lq_id = 0;
            for (const auto& lq : qubits_in) {
                max_loc_id = std::max(max_loc_id, lq.id);
                max_lq_id = std::max(max_lq_id, lq.id);
            }

            // Inicjalizacja super szybkich struktur bazowych
            fast_gate_map.resize(max_gate_id + 1);
            for (const auto& g : gates) fast_gate_map[g.id] = g;

            fast_locations.resize(max_loc_id + 1);
            initial_lq_to_loc.assign(max_lq_id + 1, -1);
            initial_loc_to_lq.assign(max_loc_id + 1, -1);

            for (const auto& lq : qubits_in) {
                int loc_id = lq.id;
                fast_locations[loc_id] = {loc_id, lq.connected_resonators};
                initial_lq_to_loc[lq.id] = loc_id;
                initial_loc_to_lq[loc_id] = lq.id;
            }

            // 2. Szybka budowa grafu na wektorach Abseila
            fast_loc_adj.resize(max_loc_id + 1);
            for (int loc1_id = 0; loc1_id <= max_loc_id; ++loc1_id) {
                if (fast_locations[loc1_id].connected_resonators.empty()) continue;

                for (int loc2_id = 0; loc2_id <= max_loc_id; ++loc2_id) {
                    if (loc1_id == loc2_id || fast_locations[loc2_id].connected_resonators.empty()) continue;

                    bool share_resonator = false;
                    const auto& res2_list = fast_locations[loc2_id].connected_resonators;
                    for (int res_id : fast_locations[loc1_id].connected_resonators) {
                        if (std::ranges::find(res2_list, res_id) != res2_list.end()) {
                            share_resonator = true;
                            break;
                        }
                    }
                    if (share_resonator) {
                        fast_loc_adj[loc1_id].push_back(loc2_id);
                    }
                }
            }

            // 3. PREKALKULACJA ŚCIEŻEK
            for (int start_loc_id = 0; start_loc_id <= max_loc_id; ++start_loc_id) {
                if (fast_locations[start_loc_id].connected_resonators.empty()) continue;

                for (const auto& res : resonators) {
                    for (int avoid_loc_id = 0; avoid_loc_id <= max_loc_id; ++avoid_loc_id) {
                        if (start_loc_id == avoid_loc_id || fast_locations[avoid_loc_id].connected_resonators.empty()) continue;

                        uint64_t key = make_cache_key(start_loc_id, res.id, avoid_loc_id);
                        path_cache[key] = compute_shortest_path(start_loc_id, res.id, avoid_loc_id, max_loc_id);
                    }
                }
            }
        }

        [[nodiscard]] uint64_t evaluate(const std::vector<std::pair<int, int>>& repr) const override {
            uint64_t total_swaps = 0;

            // Kopiowanie płaskiego wektora to mikrosekundy w porównaniu do kopiowania hash mapy z milisekundami
            std::vector<int> lq_to_loc = initial_lq_to_loc;
            std::vector<int> loc_to_lq = initial_loc_to_lq;

            // Używamy referencji do uniknięcia zbędnych kopii ścieżki
            auto apply_swaps = [&](const absl::InlinedVector<int, 16>& path) {
                if (path.empty() || path.size() == 1) return;

                if (path[0] == -1) {
                    total_swaps += 10000; // Kara za brak możliwości rutingu
                    printf("KURWAAAAAAAAAAAAAAA\n");
                    return;
                }

                for (size_t i = 0; i < path.size() - 1; ++i) {
                    int loc1 = path[i];
                    int loc2 = path[i+1];

                    int lq1 = loc_to_lq[loc1];
                    int lq2 = loc_to_lq[loc2];

                    lq_to_loc[lq1] = loc2;
                    lq_to_loc[lq2] = loc1;
                    loc_to_lq[loc1] = lq2;
                    loc_to_lq[loc2] = lq1;

                    total_swaps++;
                }
            };

            for (const auto& [gate_id, res_id] : repr) {
                const auto& gate = fast_gate_map[gate_id];

                // Najszybsza detekcja bramki 2-kubitowej (zamiast std::string == "2Q_LOCAL")
                if (gate.target_q2 != -1) {
                    int lq1 = gate.target_q1;
                    int lq2 = gate.target_q2;

                    int loc1 = lq_to_loc[lq1];
                    int loc2 = lq_to_loc[lq2];

                    // Omijamy .at() na rzecz iteratorów w absl (nie rzuca wyjątków pod spodem = duży profit do prędkosci)
                    auto it1 = path_cache.find(make_cache_key(loc1, res_id, loc2));
                    if (it1 != path_cache.end()) apply_swaps(it1->second);

                    // Aktualizacja lokalizacji
                    loc1 = lq_to_loc[lq1];

                    auto it2 = path_cache.find(make_cache_key(loc2, res_id, loc1));
                    if (it2 != path_cache.end()) apply_swaps(it2->second);
                }
            }

            return total_swaps;
        }

        [[nodiscard]] std::vector<std::pair<int, int>> random_repr(metah::ThreadLocalRandom& rng) const override {
            std::vector<std::pair<int, int>> repr;
            repr.reserve(gates.size());

            // Kolosalny speed-up: Odpytywanie hasha na wektorze std::vector<bool> (O(1) na bit)
            std::vector<bool> completed_gates(fast_gate_map.size(), false);
            std::vector<const Gate*> available_gates;

            for (const auto& gate : gates) {
                if (gate.prev_gate_q1 == nullptr && gate.prev_gate_q2 == nullptr) {
                    available_gates.push_back(&gate);
                }
            }

            auto is_ready = [&](const Gate* g) {
                if (!g) return false;
                if (g->prev_gate_q1 && !completed_gates[g->prev_gate_q1->id]) return false;
                if (g->prev_gate_q2 && !completed_gates[g->prev_gate_q2->id]) return false;
                return true;
            };

            while (!available_gates.empty()) {
                size_t idx = rng.random_idx(available_gates.size() - 1);
                const Gate* current = available_gates[idx];

                std::swap(available_gates[idx], available_gates.back());
                available_gates.pop_back();

                repr.emplace_back(current->id, rng.random_idx(resonators.size() - 1));
                completed_gates[current->id] = true;

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

        [[nodiscard]] std::vector<RoutingRecord> evaluate_and_track_swaps(const std::vector<std::pair<int, int>>& repr) const {
            std::vector<RoutingRecord> swap_history;

            std::vector<int> lq_to_loc = initial_lq_to_loc;
            std::vector<int> loc_to_lq = initial_loc_to_lq;

            auto apply_swaps = [&](const absl::InlinedVector<int, 16>& path, int current_gate_id) {
                if (path.empty() || path.size() == 1) return;

                if (path[0] == -1) {
                    printf("KURWAAAAAAAAAAAAAAA\n");
                    return;
                }

                for (size_t i = 0; i < path.size() - 1; ++i) {
                    int loc1 = path[i];
                    int loc2 = path[i+1];

                    int lq1 = loc_to_lq[loc1];
                    int lq2 = loc_to_lq[loc2];

                    // Zaktualizowano pod Enum (odkomentuj jeśli potrzebujesz zapisywać swapy)
                    // swap_history.push_back({OpType::SWAP, current_gate_id, lq1, lq2, loc1, loc2});

                    lq_to_loc[lq1] = loc2;
                    lq_to_loc[lq2] = loc1;
                    loc_to_lq[loc1] = lq2;
                    loc_to_lq[loc2] = lq1;
                }
            };

            for (const auto& [gate_id, res_id] : repr) {
                const auto& gate = fast_gate_map[gate_id];

                if (gate.target_q2 != -1) {
                    int lq1 = gate.target_q1;
                    int lq2 = gate.target_q2;

                    int loc1 = lq_to_loc[lq1];
                    int loc2 = lq_to_loc[lq2];

                    auto it1 = path_cache.find(make_cache_key(loc1, res_id, loc2));
                    if (it1 != path_cache.end()) apply_swaps(it1->second, gate_id);

                    loc1 = lq_to_loc[lq1];

                    auto it2 = path_cache.find(make_cache_key(loc2, res_id, loc1));
                    if (it2 != path_cache.end()) apply_swaps(it2->second, gate_id);
                }
            }

            return swap_history;
        }

    private:
        [[nodiscard]] absl::InlinedVector<int, 16> compute_shortest_path(int start_loc_id, int target_res_id, int avoid_loc_id, int max_loc_id) const {
            const auto& start_res_list = fast_locations[start_loc_id].connected_resonators;

            if (std::ranges::find(start_res_list, target_res_id) != start_res_list.end()) {
                return {start_loc_id};
            }

            std::queue<int> q;
            // Wyzbycie się narzutu na hash mape algorytmu Dijkstra na rzecz pre-alokacji
            std::vector<int> came_from(max_loc_id + 1, -1);
            std::vector<bool> visited(max_loc_id + 1, false);

            q.push(start_loc_id);
            came_from[start_loc_id] = start_loc_id;
            visited[start_loc_id] = true;

            int found_target_loc = -1;

            while (!q.empty()) {
                int curr = q.front(); q.pop();

                const auto& curr_res_list = fast_locations[curr].connected_resonators;
                if (std::ranges::find(curr_res_list, target_res_id) != curr_res_list.end()) {
                    found_target_loc = curr;
                    break;
                }

                for (int neighbor : fast_loc_adj[curr]) {
                    if (neighbor == avoid_loc_id) continue;

                    if (!visited[neighbor]) {
                        visited[neighbor] = true;
                        came_from[neighbor] = curr;
                        q.push(neighbor);
                    }
                }
            }

            if (found_target_loc == -1) {
                return {-1};
            }

            absl::InlinedVector<int, 16> path;
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
