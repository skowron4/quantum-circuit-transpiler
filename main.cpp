#include <iostream>
#include <memory>

#include "generator/CircuitGenerator.h"
#include "architecture/IQMConstellation.h"
#include "geneticAlgorithm/toSort/swap_evaluator.hpp"
#include "architecture/IQMTriangleConstellation.h"
#include "fmt/ranges.h"
#include "geneticAlgorithm/toSort/gate_mute.hpp"
#include "geneticAlgorithm/toSort/resonator_mute.hpp"
#include "geneticAlgorithm/toSort/resonator_gate_mute.hpp"
#include "geneticAlgorithm/toSort/ox_pair.hpp"
#include  "geneticAlgorithm/genetic_algorithm.hpp"
#include "geneticAlgorithm/toSort/tournament.hpp"


int main() {
    // 1. Inicjalizacja Topologii (dziala tak samo dla prostokata i dla trojkata)
    // auto qpu = IQMConstellation(3,3);
    auto qpu = IQMTriangleConstellation();
    qpu.print_stats();
    // 2. Inicjalizacja Generatora (podpinamy go pod wygenerowana topologie)
    CircuitGenerator generator(qpu.get_num_qubits(), qpu.get_resonators(), 1234);

    // 3. Konfiguracja parametrów obwodu
    int TOTAL_GATES = 20;
    double PROB_2Q = 1;        // 60% szans, że bramka jest 2-kubitowa
    double PROB_GLOBAL = 0.4;    // Z tych 60%, jest 40% szans że to brmaka miedzy gwiazdami (wymaga SWAP)

    cout << "\n--- Generowanie Obwodu Kwantowego ---" << endl;
    vector<Gate> my_circuit = generator.generate(TOTAL_GATES, PROB_2Q, PROB_GLOBAL);
    generator.link_gates(my_circuit);

    // 4. Analiza wygenerowanego obwodu
    int swap_needed = 0;

    for (size_t i = 0; i < my_circuit.size(); ++i) {
        Gate g = my_circuit[i];
        if (g.type == "1Q") {
            cout << i << ": [1Q] U3 na Q_" << g.target_q1 << endl;
        } else if (g.type == "2Q_LOCAL") {
            cout << i << ": [2Q LOKALNA] CZ na Q_" << g.target_q1 << " i Q_" << g.target_q2 << " (Bezpośrednie)" << endl;
        } else if (g.type == "2Q_GLOBAL") {
            cout << i << ": [2Q GLOBALNA] CZ na Q_" << g.target_q1 << " i Q_" << g.target_q2 << " !!! WYMAGA ROUTINGU !!!" << endl;
            swap_needed++;
        }
    }

    cout << "\nPodsumowanie: Wygenerowano " << swap_needed << " bramek wymagajacych kompilacji i SWAPow." << endl;

    generator.print_circuit(my_circuit);

    int i = 0;
    cout << my_circuit[i].id << ": " << my_circuit[i].type << " na Q_" << my_circuit[i].target_q1 << ", na Q_" << my_circuit[i].target_q2 << " prev1 " <<my_circuit[i].prev_gate_q1 << " prev2 " <<my_circuit[i].prev_gate_q2<< " next1 " <<my_circuit[i].next_gate_q1<< " next2 " <<my_circuit[i].next_gate_q2<< endl;

    fmt::println("resonators {}",qpu.resonators);
    fmt::print("=====================================================\n");
    fmt::println("qubit {}",qpu.qubits);

    // Gate g1 = {0, "1Q", 1, -1};
    // Gate g2 = {1, "2Q_LOCAL", 0, 2};
    // Gate g3 = {2, "1Q", 2, -1};
    //
    // g2.next_gate_q2 = &g3;
    // g3.prev_gate_q1 = &g2;
    //
    // auto gates = std::vector<Gate>{g1,g2,g3};
    // auto x = std::vector<std::pair<int,int>>{{0,1},{1,1},{2,1},{3,1}};
    // auto xc = std::vector<std::pair<int,int>>{{0,1},{1,1},{2,1},{3,1}};
    // auto y = std::vector<std::pair<int,int>>{{3,2},{2,2},{1,2},{0,2}};
    // auto yc = std::vector<std::pair<int,int>>{{3,2},{2,2},{1,2},{0,2}};

    auto evaluator = std::make_shared<metah::SwapEvaluator>(my_circuit, qpu.resonators, qpu.qubits);
    // auto mutateRes = std::make_unique<metah::ResonatorMute<>>(0.5, qpu.resonators.size());
    // auto mutateGate = std::make_unique<metah::GateMute<>>(1, gates);
    auto mutateResGate = std::make_unique<metah::ResonatorGateMute<>>(0.05, 0.1, qpu.resonators.size(), my_circuit);
    auto ox = std::make_unique<metah::ox_pair<>>(0.7);
    auto tournament = std::make_unique<metah::Tournament<uint64_t>>(5);

    metah::GeneticAlgorithm<std::vector<std::pair<int,int>>, uint64_t> ga(
        evaluator,
        300,
        6000 ,
        12,
        std::move(tournament),
        std::move(ox),
        std::move(mutateResGate)
    );


    ga.run();
    auto [repr, score] = ga.getBest();

    fmt::println("swaps {}", evaluator->evaluate_and_track_swaps(repr));
    fmt::println("Best solution: {} with score {}", repr, score);

    // fmt::println("gates {}",evaluator->gates);
    // fmt::println("resonators {}",evaluator->resonators);
    // fmt::println("gate_map {}",evaluator->gate_map);
    // fmt::println("locations {}",evaluator->locations);
    // fmt::println("loc_adj {}",evaluator->loc_adj);
    // fmt::println("path_cache {}",evaluator->path_cache);
    // fmt::println("initial_lq_to_loc {}",evaluator->initial_lq_to_loc);

    // fmt::println("evaluate {}", evaluator->evaluate(x));
    // fmt::println("repres {}", x);

    // fmt::println("muteRes {}", mutateRes->mutate(x));
    // fmt::println("repres {}", x);
    // fmt::println("muteRes {}", mutateRes->mutate(x));
    // fmt::println("repres {}", x);

    // fmt::println("muteGate {}", mutateGate->mutate(x));
    // fmt::println("repres {}", x);
    // fmt::println("muteGate {}", mutateGate->mutate(x));
    // fmt::println("repres {}", x);

    // fmt::println("muteResGate {}", mutateResGate->mutate(x));
    // fmt::println("repres {}", x);
    // fmt::println("muteResGate {}", mutateResGate->mutate(x));
    // fmt::println("repres {}", x);

    // fmt::println("muteResGate {}", ox->cross(x,y,xc,yc));
    // fmt::println("repres x {}", x);
    // fmt::println("repres y {}", y);
    // fmt::println("repres xc {}", xc);
    // fmt::println("repres yc {}", yc);
    // fmt::println("muteResGate {}", ox->cross(xc,yc,x,y));
    // fmt::println("repres x {}", xc);
    // fmt::println("repres y {}", yc);
    // fmt::println("repres xc {}", x);
    // fmt::println("repres yc {}", y);

    return 0;
}
