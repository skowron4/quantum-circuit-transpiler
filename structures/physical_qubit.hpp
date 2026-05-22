#ifndef PHYSQUBIT_HPP
#define PHYSQUBIT_HPP
#include <absl/container/inlined_vector.h>

struct PhysicalQubit {
    int id;
    absl::InlinedVector<int, 5> connected_resonators;
};

#include <fmt/format.h>
#include <fmt/ranges.h>

template <>
struct fmt::formatter<PhysicalQubit> {
    constexpr auto parse(format_parse_context& ctx) {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const PhysicalQubit& q, FormatContext& ctx) const {
        return fmt::format_to(
            ctx.out(),
            "resonator(id={}, qubit={})",
            q.id,
            q.connected_resonators
        );
    }
};

#endif //PHYSQUBIT_HPP
