#ifndef GATE_H
#define GATE_H

#include <sstream>


class Gate {
public:
    int id;
    std::string type;     // "1Q", "2Q_LOCAL", "2Q_GLOBAL"
    int target_q1;
    int target_q2;   // -1 dla bramek 1-kubitowych

    Gate* prev_gate_q1 = nullptr;
    Gate* prev_gate_q2 = nullptr;
    Gate* next_gate_q1 = nullptr;
    Gate* next_gate_q2 = nullptr;
};

#include <fmt/format.h>

template <>
struct fmt::formatter<Gate> {
    constexpr auto parse(format_parse_context& ctx) {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const Gate& g, FormatContext& ctx) const {
        return fmt::format_to(ctx.out(), "Gate(id={}, type={}, q1={}, q2={})", g.id, g.type, g.target_q1, g.target_q2);
    }
};

#endif //GATE_H
