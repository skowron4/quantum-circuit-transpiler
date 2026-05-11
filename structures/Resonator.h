#ifndef RESONATOR_H
#define RESONATOR_H

#include <set>

class Resonator {
public:
    int id;
    int y, x;
    std::set<int> connected_qubits;
};

#include <fmt/format.h>
#include <fmt/ranges.h>

template <>
struct fmt::formatter<Resonator> {
    constexpr auto parse(format_parse_context& ctx) {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const Resonator& q, FormatContext& ctx) const {
        return fmt::format_to(
            ctx.out(),
            "resonator(id={}, qubit={})",
            q.id,
            q.connected_qubits
        );
    }
};

#endif //RESONATOR_H
