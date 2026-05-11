#ifndef METAH_OPTIMIZER_HPP
#define METAH_OPTIMIZER_HPP
#include "solution.hpp"

namespace metah {
    template <typename R, typename S>
    class Optimizer {
    protected:
        std::shared_ptr<const Problem<R, S>> evaluator;

    public:
        explicit Optimizer(std::shared_ptr<const Problem<R, S>> evaluator) : evaluator(std::move(evaluator)) {}
        virtual void run() = 0;
        [[nodiscard]] virtual Solution<R, S> getBest() const = 0;
        virtual ~Optimizer() = default;
    };
}
#endif //METAH_OPTIMIZER_HPP
