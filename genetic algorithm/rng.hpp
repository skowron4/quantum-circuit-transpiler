#ifndef METAH_RNG_HPP
#define METAH_RNG_HPP
#include <random>
#include <boost/random/xoshiro.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

namespace br = boost::random;

namespace metah {
    class ThreadLocalRandom {
    public:
        ThreadLocalRandom() = default;

        static size_t random_idx(const size_t max) {
            const br::uniform_int_distribution<size_t> dist(0, max);
            return dist(randomEngine());
        }

        static std::ptrdiff_t random_offset(const std::ptrdiff_t min, const std::ptrdiff_t max) {
            const br::uniform_int_distribution dist(min, max);
            return dist(randomEngine());
        }

        static double random01() {
            const br::uniform_real_distribution dist(0.0, 1.0);
            return dist(randomEngine());
        }

        static br::xoshiro256mm& randomEngine() {
            thread_local br::xoshiro256mm rng(std::random_device{}());
            return rng;
        }
    };
}

#endif //METAH_RNG_HPP
