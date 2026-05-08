//
// Created by igorb on 8.05.2026.
//

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



#endif //GATE_H
