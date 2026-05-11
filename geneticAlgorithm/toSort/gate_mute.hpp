#ifndef GATE_MUTE_HPP
#define GATE_MUTE_HPP



#include "geneticAlgorithm/rng.hpp"
#include "../templates/mutate.hpp"
#include "structures/Gate.h"

namespace metah {
    template <typename RNG = ThreadLocalRandom>
    class GateMute final : public Mutate<std::vector<std::pair<int, int>>> {
        double chance;
        RNG rng;
        const std::vector<Gate>& gates;

    public:
        template <typename... Args>
        explicit GateMute(const double chance,const std::vector<Gate>& gates ,Args&&... args) : chance(chance), gates(gates), rng(std::forward<Args>(args)...) {}

        bool mutate(std::vector<std::pair<int, int>>& vec) override {
            if (rng.random01() > chance) { return false; }

            // 1. Losujemy indeks elementu, który będziemy przesuwać
            int src_idx = rng.random_idx(vec.size() - 1);
            int target_id = vec[src_idx].first;

            // 2. Szukamy odpowiedniego obiektu Gate w referencyjnym wektorze
            const Gate* target_gate = nullptr;
            for (const auto& g : gates) {
                if (g.id == target_id) {
                    target_gate = &g;
                    break;
                }
            }

            if (!target_gate) { return false; }

            auto find_idx = [&](int id_to_find) -> int {
                for (size_t i = 0; i < vec.size(); ++i) {
                    if (vec[i].first == id_to_find) return static_cast<int>(i);
                }
                return -1;
            };

            // 3. Wyznaczamy dozwolone "okno" indeksów (granice) dla wybranej bramki
            int lower_bound = 0;
            int upper_bound = static_cast<int>(vec.size() - 1);

            // Dolna granica to największy z indeksów poprzedników + 1
            if (target_gate->prev_gate_q1) {
                int p_idx = find_idx(target_gate->prev_gate_q1->id);
                if (p_idx != -1) lower_bound = std::max(lower_bound, p_idx + 1);
            }
            if (target_gate->prev_gate_q2) {
                int p_idx = find_idx(target_gate->prev_gate_q2->id);
                if (p_idx != -1) lower_bound = std::max(lower_bound, p_idx + 1);
            }

            // Górna granica to najmniejszy z indeksów następników - 1
            if (target_gate->next_gate_q1) {
                int n_idx = find_idx(target_gate->next_gate_q1->id);
                if (n_idx != -1) upper_bound = std::min(upper_bound, n_idx - 1);
            }
            if (target_gate->next_gate_q2) {
                int n_idx = find_idx(target_gate->next_gate_q2->id);
                if (n_idx != -1) upper_bound = std::min(upper_bound, n_idx - 1);
            }

            // 4. Obliczamy liczbę dozwolonych nowych miejsc
            int choices = upper_bound - lower_bound;
            if (choices <= 0) {
                // Jeśli wybór wynosi 0, oznacza to, że bramka jest "uwięziona"
                // i jedyne legalne miejsce to jej aktualna pozycja.
                return false;
            }

            // Losujemy nowe miejsce z dostępnej puli omijając obecną pozycję (src_idx)
            // Zakładam, że random_idx(n) losuje od 0 do n włącznie.
            int random_choice = rng.random_idx(choices - 1);
            int new_idx = lower_bound + random_choice;

            if (new_idx >= src_idx) {
                new_idx++;
            }

            // 5. Przesuwamy element w wektorze używając std::rotate
            if (new_idx < src_idx) {
                // Przesunięcie w lewo
                std::rotate(vec.begin() + new_idx, vec.begin() + src_idx, vec.begin() + src_idx + 1);
            } else if (new_idx > src_idx) {
                // Przesunięcie w prawo
                std::rotate(vec.begin() + src_idx, vec.begin() + src_idx + 1, vec.begin() + new_idx + 1);
            }

            return true;
        }

        ~GateMute() override = default;
    };
}

#endif //GATE_MUTE_HPP


