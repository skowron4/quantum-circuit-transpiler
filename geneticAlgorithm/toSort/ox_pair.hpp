#ifndef OX_PAIR_HPP
#define OX_PAIR_HPP

#include <vector>
#include <unordered_set>

#include "../templates/crossover.hpp"
#include "../rng.hpp"

namespace metah {
    template <typename RNG = ThreadLocalRandom>
    class ox_pair final : public Crossover<std::vector<std::pair<int,int>>> {
        double chance;
        RNG rng;

    public:
        template <typename... Args>
        explicit ox_pair(const double chance, Args&&... args)
            : chance(chance), rng(std::forward<Args>(args)...) {}

        bool cross(const std::vector<std::pair<int,int>>& p1,
                   const std::vector<std::pair<int,int>>& p2,
                   std::vector<std::pair<int,int>>& c1,
                   std::vector<std::pair<int,int>>& c2) override {

            if (rng.random01() > chance) {
                c1 = p1;
                c2 = p2;
                return false;
            }

            const size_t size = p1.size();
            c1.resize(size);
            c2.resize(size);

            std::unordered_set<int> used1;
            std::unordered_set<int> used2;

            auto begin = rng.random_offset(0, size - 1);
            auto end   = rng.random_offset(0, size - 1);
            if (begin > end) std::swap(begin, end);
            ++end;

            // 1. kopiuj segment
            for (size_t i = begin; i < end; ++i) {
                c1[i] = p1[i];
                c2[i] = p2[i];
                used1.insert(p1[i].first);
                used2.insert(p2[i].first);
            }

            // 2. uzupełnij c1 z p2
            size_t idx = end % size;

            for (size_t i = end; i < size; ++i) {
                while (used1.contains(p2[idx].first)) {
                    if (++idx == size) idx = 0;
                }
                c1[i] = p2[idx];
                used1.insert(p2[idx].first);
            }

            for (size_t i = 0; i < begin; ++i) {
                while (used1.contains(p2[idx].first)) {
                    if (++idx == size) idx = 0;
                }
                c1[i] = p2[idx];
                used1.insert(p2[idx].first);
            }

            // 3. uzupełnij c2 z p1
            idx = end % size;

            for (size_t i = end; i < size; ++i) {
                while (used2.contains(p1[idx].first)) {
                    if (++idx == size) idx = 0;
                }
                c2[i] = p1[idx];
                used2.insert(p1[idx].first);
            }

            for (size_t i = 0; i < begin; ++i) {
                while (used2.contains(p1[idx].first)) {
                    if (++idx == size) idx = 0;
                }
                c2[i] = p1[idx];
                used2.insert(p1[idx].first);
            }

            return true;
        }
    };
}

#endif // OX_PAIR_HPP