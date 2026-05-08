//
// Created by igorb on 8.05.2026.
//

#ifndef GATE_H
#define GATE_H

#include <sstream>


class Gate {
public:
    std::string type;     // "1Q", "2Q_LOCAL", "2Q_GLOBAL"
    int target_q1;
    int target_q2;   // -1 dla bramek 1-kubitowych
};



#endif //GATE_H
