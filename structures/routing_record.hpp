

#ifndef SWAP_RECORD_HPP
#define SWAP_RECORD_HPP

#include <string>
#include <fmt/core.h>

enum class OpType { MOVE, SWAP, GATE_1Q, GATE_2Q_LOCAL, GATE_2Q_GLOBAL };

inline std::string to_string(OpType type) {
    switch (type) {
        case OpType::MOVE: return "MOVE";
        case OpType::SWAP: return "SWAP";
        case OpType::GATE_1Q: return "1Q";
        case OpType::GATE_2Q_LOCAL: return "2Q_LOCAL";
        case OpType::GATE_2Q_GLOBAL: return "2Q_GLOBAL";
        default: return "UNKNOWN";
    }
}

struct RoutingRecord {
    OpType op_type; // "move" lub "swap"
    int gate_id;
    int lq1;
    int lq2;
    int loc1;
    int loc2;
    int res_used = -1;       // Przekazanie wybranego rezonatora z DE do SWAP/MOVE
    double start_time = 0.0;
};

template <>
struct fmt::formatter<RoutingRecord> {
    constexpr auto parse(fmt::format_parse_context& ctx) { return ctx.begin(); }

    template <typename FormatContext>
    auto format(const RoutingRecord& q, FormatContext& ctx) const {
        return fmt::format_to(
            ctx.out(),
            "{}(gate_id={}, lq1={}, lq2={}, loc1={}, loc2={})",
            to_string(q.op_type), q.gate_id, q.lq1, q.lq2, q.loc1, q.loc2
        );
    }
};


#endif //SWAP_RECORD_HPP
