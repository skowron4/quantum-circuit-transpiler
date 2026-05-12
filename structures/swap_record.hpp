

#ifndef SWAP_RECORD_HPP
#define SWAP_RECORD_HPP

#include <fmt/format.h>
#include <fmt/ranges.h>

struct SwapRecord {
    int gate_id;
    int lq1; // Logiczny kubit 1
    int lq2; // Logiczny kubit 2
    int loc1; // Fizyczna lokalizacja 1 (opcjonalnie)
    int loc2; // Fizyczna lokalizacja 2 (opcjonalnie)
};


template <>
struct fmt::formatter<SwapRecord> {
    constexpr auto parse(fmt::format_parse_context& ctx) {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const SwapRecord& q, FormatContext& ctx) const {
        return fmt::format_to(
            ctx.out(),
            "swap(gate_id={}, lq1={}, lq2={}, loc1={}, loc2={})",
            q.gate_id,
            q.lq1,
            q.lq2,
            q.loc1,
            q.loc2
        );
    }
};


#endif //SWAP_RECORD_HPP
