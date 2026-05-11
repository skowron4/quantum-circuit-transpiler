#ifndef METAH_TOURNAMENT_HPP
#define METAH_TOURNAMENT_HPP
#include <span>

#include "../rng.hpp"
#include "../templates/selection.hpp"

namespace metah {
    template <typename S, typename RNG = ThreadLocalRandom>
    class Tournament final : public Selection<S> {
        size_t size;
        RNG rng;

    public:
        template <typename... Args>
        explicit Tournament(const size_t size, Args&&... args) : size(size), rng(std::forward<Args>(args)...) {}

        size_t select(std::span<const S> scores) override {
            assert(!scores.empty());
            size_t best_idx = rng.random_idx(scores.size() - 1);
            for (size_t i = 1; i < size; ++i) {
                const auto idx = rng.random_idx(scores.size() - 1);
                if (scores[idx] < scores[best_idx]) {
                    best_idx = idx;
                }
            }
            return best_idx;
        }
    };
}

#endif //METAH_TOURNAMENT_HPP
