#ifndef QUBIT_H
#define QUBIT_H

#include <set>



class Qubit {
public:
    int id;
    std::set<int> connected_resonators;
};

#include <fmt/format.h>
#include <fmt/ranges.h>

template <>
struct fmt::formatter<Qubit> {
    constexpr auto parse(format_parse_context& ctx) {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const Qubit& q, FormatContext& ctx) const {
        return fmt::format_to(
            ctx.out(),
            "Qubit(id={}, resonators={})",
            q.id,
            q.connected_resonators
        );
    }
};

#endif //QUBIT_H
