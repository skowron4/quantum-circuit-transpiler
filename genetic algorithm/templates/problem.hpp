#ifndef METAH_EVALUATOR_HPP
#define METAH_EVALUATOR_HPP
#include "../rng.hpp"
#include "solution.hpp"

namespace metah {
    template <typename R, typename S, typename RNG = ThreadLocalRandom>
    class Problem {
    public:
        [[nodiscard]] virtual S evaluate(const R& repr) const = 0;
        [[nodiscard]] virtual R random_repr(RNG& rng) const = 0;

        [[nodiscard]] virtual Solution<R, S> random_solution(RNG& rng) const {
            auto repr = random_repr(rng);
            return {repr, evaluate(repr)};
        }

        virtual ~Problem() = default;
    };
}

#endif //METAH_EVALUATOR_HPP
