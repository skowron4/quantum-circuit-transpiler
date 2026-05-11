#ifndef METAH_MUTATE_HPP
#define METAH_MUTATE_HPP

namespace metah {
    template <typename R>
    class Mutate {
    public:
        virtual bool mutate(R& repr) = 0;
        virtual ~Mutate() = default;
    };
}


#endif //METAH_MUTATE_HPP
