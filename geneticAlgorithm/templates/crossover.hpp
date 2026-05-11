#ifndef METAH_CROSSOVER_HPP
#define METAH_CROSSOVER_HPP

namespace metah {
    template <typename R>
    class Crossover {
    public:
        virtual bool cross(const R& p1, const R& p2, R& c1, R& c2) = 0;
        virtual ~Crossover() = default;
    };
}

#endif //METAH_CROSSOVER_HPP
