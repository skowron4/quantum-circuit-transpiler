
#ifndef RESONATOR_MUTE_HPP
#define RESONATOR_MUTE_HPP

#include "geneticAlgorithm/rng.hpp"
#include "../templates/mutate.hpp"

namespace metah {
    template <typename RNG = ThreadLocalRandom>
    class ResonatorMute final : public Mutate<std::vector<std::pair<int, int>>> {
        double chance;
        RNG rng;
        uint64_t resonator_size;

    public:
        template <typename... Args>
        explicit ResonatorMute(const double chance,const uint64_t resonator_size ,Args&&... args) : chance(chance), resonator_size(resonator_size), rng(std::forward<Args>(args)...) {}

        bool mutate(std::vector<std::pair<int, int>>& vec) override {
            if (rng.random01() > chance) { return false; }

            auto idx = rng.random_idx(vec.size() - 1);
            auto idr = rng.random_idx(resonator_size -2);
            if (vec[idx].second <= idr) { ++idr; }

            vec[idx].second = idr;

            return true;
        }

        ~ResonatorMute() override = default;
    };
}

#endif //RESONATOR_MUTE_HPP
