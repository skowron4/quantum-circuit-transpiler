

#ifndef SWAP_RECORD_HPP
#define SWAP_RECORD_HPP

#include <string>
#include <fmt/core.h>

struct RoutingRecord {
    std::string op_type; // "move" lub "swap"
    int gate_id;
    int lq1;
    int lq2;
    int loc1;
    int loc2;
};

template <>
struct fmt::formatter<RoutingRecord> {
    constexpr auto parse(fmt::format_parse_context& ctx) { return ctx.begin(); }

    template <typename FormatContext>
    auto format(const RoutingRecord& q, FormatContext& ctx) const {
        return fmt::format_to(
            ctx.out(),
            "{}(gate_id={}, lq1={}, lq2={}, loc1={}, loc2={})",
            q.op_type, q.gate_id, q.lq1, q.lq2, q.loc1, q.loc2
        );
    }
};


#endif //SWAP_RECORD_HPP
