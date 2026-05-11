//
// Created by igorb on 11.05.2026.
//

#ifndef RESONATOR_GATE_MUTE_HPP
#define RESONATOR_GATE_MUTE_HPP

#include <memory>

#include "geneticAlgorithm/rng.hpp"
#include "../templates/mutate.hpp"
#include "resonator_mute.hpp"
#include "gate_mute.hpp"

namespace metah {
    template <typename RNG = ThreadLocalRandom>
    class ResonatorGateMute final : public Mutate<std::vector<std::pair<int, int>>> {
        RNG rng;
        GateMute<> gateMute;
        ResonatorMute<> resonatorMute;

    public:
        template <typename... Args>
        explicit ResonatorGateMute(const double chanceGate, const double chanceResonator, const uint64_t resonator_size, const std::vector<Gate>& gates, Args&&... args)
            : rng(std::forward<Args>(args)...),
              gateMute(chanceGate, gates),
              resonatorMute(chanceResonator, resonator_size)
        {}

        bool mutate(std::vector<std::pair<int, int>>& vec) override {
            return gateMute.mutate(vec) | resonatorMute.mutate(vec);
        }

        ~ResonatorGateMute() override = default;
    };
}




#endif //RESONATOR_GATE_MUTE_HPP
