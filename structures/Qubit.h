//
// Created by igorb on 8.05.2026.
//

#ifndef QUBIT_H
#define QUBIT_H

#include <set>



class Qubit {
public:
    int id;
    std::set<int> connected_resonators;
};



#endif //QUBIT_H
