#ifndef SELECTION_HPP
#define SELECTION_HPP

#include <span>

namespace metah {
    template <typename S>
    class Selection {
    public:
        virtual size_t select(std::span<const S> scores) = 0;

        virtual ~Selection() = default;
    };
}

#endif //SELECTION_HPP
