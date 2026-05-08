//
// Created by igorb on 8.05.2026.
//

#ifndef RESONATOR_H
#define RESONATOR_H

#include <set>

class Resonator {
public:
    int id;
    int y, x;
    std::set<int> connected_qubits;
};



#endif //RESONATOR_H
