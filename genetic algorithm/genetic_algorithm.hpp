#ifndef GENETIC_ALGORITHM_HPP
#define GENETIC_ALGORITHM_HPP

#include <memory>
#include <vector>
#include <omp.h>

#include <boost/align/aligned_allocator.hpp>
#include "templates/problem.hpp"
#include "templates/optimizer.hpp"
#include "templates/selection.hpp"
#include "templates/crossover.hpp"
#include "templates/mutate.hpp"

namespace metah {
    constexpr size_t POPULATION_ALIGNMENT = 64;

    template <
        typename R,
        typename S,
        typename RNG = ThreadLocalRandom>
    class GeneticAlgorithm : public Optimizer<R, S> {
        size_t population_size;
        size_t generations;
        size_t threads;
        std::unique_ptr<Selection<S>> select;
        std::unique_ptr<Crossover<R>> cross;
        std::unique_ptr<Mutate<R>> mutate;

        std::vector<R, boost::alignment::aligned_allocator<R, POPULATION_ALIGNMENT>> population;
        std::vector<S, boost::alignment::aligned_allocator<S, POPULATION_ALIGNMENT>> scores;
        Solution<R, S> best_solution;
        RNG rng;

        void update_best_solution() {
            size_t best_idx = 0;
            auto best_score = scores[0];
            for (size_t j = 1; j < population_size; ++j) {
                if (scores[j] < best_score) {
                    best_score = scores[j];
                    best_idx = j;
                }
            }
            if (best_score < best_solution.score) {
                best_solution.repr = population[best_idx];
                best_solution.score = best_score;
            }
        }

        void init_population() {
            if (population.size() < population_size) {
                population.reserve(population_size);
                scores.reserve(population_size);
                for (size_t i = population.size(); i < population_size; ++i) {
                    population.emplace_back(this->evaluator->random_repr(rng));
                    scores.emplace_back(this->evaluator->evaluate(population[i]));
                }
            }
            else if (population.size() > population_size) {
                population.resize(population_size);
                scores.resize(population_size);
            }
            update_best_solution();
        }

    public:
        template <typename... Args>
        GeneticAlgorithm(
            std::shared_ptr<const Problem<R, S>> evaluator,
            const size_t population_size,
            const size_t generations,
            const size_t threads,
            std::unique_ptr<Selection<S>> select,
            std::unique_ptr<Crossover<R>> cross,
            std::unique_ptr<Mutate<R>> mutate,
            Args&&... args
        ) : Optimizer<R, S>(std::move(evaluator)),
            population_size(population_size - population_size % 2),
            generations(generations),
            threads(threads),
            select(std::move(select)),
            cross(std::move(cross)),
            mutate(std::move(mutate)),
            best_solution(this->evaluator->random_solution(rng)),
            rng(std::forward<Args>(args)...) {}

        [[nodiscard]] Solution<R, S> getBest() const override {
            return best_solution;
        }

        void run() override {
            init_population();
            auto next_pop = population;
            auto next_scores = scores;
            #pragma omp parallel num_threads(threads) default(none) shared(next_pop, next_scores)
            {
                for (size_t i = 0; i < generations; ++i) {
                    #pragma omp for schedule(dynamic, 4)
                    for (size_t j = 0; j < population_size; j += 2) {
                        const size_t idx1 = select->select(scores);
                        const size_t idx2 = select->select(scores);
                        auto& child1 = next_pop[j];
                        auto& child2 = next_pop[j + 1];
                        auto modified1 = cross->cross(population[idx1], population[idx2],
                                                      child1, child2);
                        auto modified2{modified1};
                        modified1 = mutate->mutate(child1) || modified1;
                        modified2 = mutate->mutate(child2) || modified2;
                        next_scores[j] = modified1 ? this->evaluator->evaluate(child1) : scores[idx1];
                        next_scores[j + 1] = modified2 ? this->evaluator->evaluate(child2) : scores[idx2];
                    }
                    #pragma omp master
                    {
                        std::swap(population, next_pop);
                        std::swap(scores, next_scores);
                        update_best_solution();
                    }
                    #pragma omp barrier
                }
            }
        }

        ~GeneticAlgorithm() override = default;
    };
}

#endif //GENETIC_ALGORITHM_HPP
