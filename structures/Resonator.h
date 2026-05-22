#ifndef RESONATOR_H
#define RESONATOR_H

#include <absl/container/inlined_vector.h>

class Resonator {
public:
    int id;
    int y, x;
    absl::InlinedVector<int, 10> connected_qubits;
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
            "resonator(id={}, qubit={}, x={}, y={})",
            q.id,
            q.connected_qubits,
            q.x,
            q.y
        );
    }
};

#endif //RESONATOR_H
