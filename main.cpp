#include <iostream>

#include "CircuitGenerator.h"
#include "IQMConstellation.h"
#include "IQMTriangleConstellation.h"

// int main() {
//     std::cout << "Hello, World!" << std::endl;
//     return 0;
// }

// int main() {
//     // Inicjalizacja modelu z 3 rezonatorami
//     IQMTriangleConstellation qpu_triangle;
//
//     // Wypisanie statystyk, aby udowodnić liczbę poszczególnych typów kubitów
//     qpu_triangle.print_detailed_stats();
//
//     int start_node = 0;
//     int target_node = 2;
//
//     // O(1) zapytanie o koszt
//     int cost = qpu_triangle.get_swap_cost(start_node, target_node);
//     cout << "\nKoszt SWAP z Res 0 do Res 2: " << cost << " skok(ów)\n" << endl;
//
//     // Pobranie dostepnych kubitów dla tego skoku
//     vector<int> swap_qubits = qpu_triangle.get_shared_qubits(start_node, target_node);
//     cout << "Dostepne kubity umozliwiajace skok Res 0 -> Res 2 (" << swap_qubits.size() << " opcje):" << endl;
//     for (int q : swap_qubits) {
//         cout << " -> Kubit ID " << q << endl;
//     }
//
//     return 0;
// }

// int main() {
//     int Y = 3;
//     int X = 3;
//     IQMConstellation qpu(Y, X);
//
//     qpu.print_stats();
//
//     int start_node = 0;
//     int target_node = 8;
//
//     // O(1) zapytanie o koszt
//     int cost = qpu.get_swap_cost(start_node, target_node);
//     cout << "\nMinimalna liczba skokow (koszt SWAP): " << cost << endl;
//
//     // Błyskawiczne odtworzenie ścieżki
//     vector<int> path = qpu.get_resonator_path(start_node, target_node);
//
//     if (path.empty()) {
//         cout << "Sciezka nie istnieje!" << endl;
//     } else {
//         cout << "Trasa zrekonstruowana! Skoki:\n";
//         for (size_t i = 0; i < path.size() - 1; ++i) {
//             int r_from = path[i];
//             int r_to = path[i+1];
//
//             vector<int> swap_qubits = qpu.get_shared_qubits(r_from, r_to);
//
//             cout << "Skok " << i+1 << ": Res " << r_from << " -> Res " << r_to << " | "
//                  << "Dostepne kubity do SWAP (" << swap_qubits.size() << "): ";
//             for (int q : swap_qubits) cout << "Q_" << q << " ";
//             cout << endl;
//         }
//     }
//
//     return 0;
// }

int main() {
    // 1. Inicjalizacja Topologii (dziala tak samo dla prostokata i dla trojkata)
    IQMConstellation qpu(3, 3);
    qpu.print_stats();

    // 2. Inicjalizacja Generatora (podpinamy go pod wygenerowana topologie)
    CircuitGenerator generator(qpu.get_num_qubits(), qpu.get_resonators(), 1234);

    // 3. Konfiguracja parametrów obwodu
    int TOTAL_GATES = 20;
    double PROB_2Q = 0.6;        // 60% szans, że bramka jest 2-kubitowa
    double PROB_GLOBAL = 0.4;    // Z tych 60%, jest 40% szans że to brmaka miedzy gwiazdami (wymaga SWAP)

    cout << "\n--- Generowanie Obwodu Kwantowego ---" << endl;
    vector<Gate> my_circuit = generator.generate(TOTAL_GATES, PROB_2Q, PROB_GLOBAL);

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

    return 0;
}