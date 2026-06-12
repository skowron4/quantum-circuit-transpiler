#include <fstream>
#include <iostream>
#include <cstdint>
#include "generator/CircuitGenerator.h"
#include "geneticAlgorithm/toSort/swap_evaluator.hpp"
#include "architecture/IQMTriangleConstellation.h"
#include "dEHybridSabre/SabreDEHybridRouter.h"
#include "newDESABRE/newSabreDEHybridRouter.h"
#include "differentialEvolution/differential_evolution.h"
#include "sabre/SabreRouter.h"
#include <iomanip> // Potrzebne do formatowania (np. setw)
#include <nlohmann/json.hpp>
#include <random>
#ifdef _WIN32
#include <windows.h>
#endif
#include <filesystem>
#include "greedy/GreedyASAP.h"

class CircuitIO {
public:
    static std::vector<Gate> load_circuit(const std::string& filename, int& out_num_qubits) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Blad otwarcia pliku: " << filename << "\n";
            return {};
        }

        nlohmann::json j; file >> j;
        out_num_qubits = j["num_qubits"];
        std::vector<Gate> circuit;

        for (const auto& item : j["gates"]) {
            Gate g;
            g.id = item["id"];
            g.type = item["type"];

            auto qubits = item["qubits"];
            g.target_q1 = qubits[0];
            if (qubits.size() > 1) g.target_q2 = qubits[1];

            circuit.push_back(g);
        }
        return circuit;
    }

    static std::vector<Gate> load_circuitII(const std::string& filename, int& out_num_qubits) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Blad otwarcia pliku: " << filename << "\n";
            return {};
        }

        nlohmann::json j;
        try {
            file >> j;
        } catch (const nlohmann::json::parse_error& e) {
            std::cerr << "Blad parsowania JSON w pliku " << filename << ": " << e.what() << "\n";
            return {};
        }

        if (!j.contains("num_qubits") || !j.contains("gates")) {
            std::cerr << "Plik " << filename << " nie ma wymaganych pol (num_qubits, gates).\n";
            return {};
        }

        out_num_qubits = j["num_qubits"];
        std::vector<Gate> circuit;

        for (const auto& item : j["gates"]) {
            Gate g;
            g.id = item["id"];

            // Bezpieczne przypisanie stringa
            g.type = item["type"].get<std::string>();

            // Zawsze bezpiecznie inicjalizuj brakujące kubity na -1
            g.target_q1 = -1;
            g.target_q2 = -1;

            auto qubits = item["qubits"];

            // Zabezpieczenie: przypisuj tylko, gdy JSON faktycznie posiada te indeksy
            if (qubits.size() > 0) g.target_q1 = qubits[0];
            if (qubits.size() > 1) g.target_q2 = qubits[1];

            circuit.push_back(g);
        }

        if (circuit.empty()) {
            std::cerr << "Ostrzezenie: Plik " << filename << " wczytal sie, ale nie ma zadnych bramek!\n";
        }

        return circuit;
    }

    static void save_circuit(const std::string& filename, const std::vector<RoutingRecord>& log) {
        nlohmann::json j;
        j["gates"] = nlohmann::json::array();

        for (const auto& r : log) {
            nlohmann::json g_json;
            if (r.op_type == OpType::SWAP || r.op_type == OpType::MOVE) {
                g_json["id"] = -1;
                g_json["type"] = (r.op_type == OpType::SWAP) ? "SWAP" : "MOVE";
                // ZMIANA TUTAJ: loc1 i loc2 zamiast phys_q1 i phys_q2
                g_json["qubits"] = {r.loc1, r.loc2};
            } else {
                g_json["id"] = r.gate_id;
                g_json["type"] = (r.op_type == OpType::GATE_1Q) ? "1Q" : "2Q";
                // ZMIANA TUTAJ: loc1 i loc2
                if (r.op_type == OpType::GATE_1Q) g_json["qubits"] = {r.loc1};
                else g_json["qubits"] = {r.loc1, r.loc2};
            }
            j["gates"].push_back(g_json);
        }

        std::ofstream file(filename);
        if (file.is_open()) file << j.dump(4);
    }
};


void run_test_for_topology(const BaseTopology &qpu, const std::string &name) {
    std::cout << "\n============================================\n";
    std::cout << "URUCHAMIAM TEST DLA TOPOLOGII: " << name << "\n";
    std::cout << "============================================\n";
    qpu.print_stats();

    int logical_qubits_count = qpu.get_num_qubits() * 1.0;
    CircuitGenerator generator(qpu, logical_qubits_count, 1234);

    // Generowanie obwodu testowego
    std::vector<Gate> circuit = generator.generate(2000, 0.6, 0.4);
    generator.print_gate_statistics(circuit);
    generator.print_circuit(circuit);

    // --- Konfiguracja losowosci i liczby prob ---
    const bool use_random_seeds = false;
    std::random_device rd;
    const uint32_t sabre_seed = rd();
    const uint64_t hybrid_seed = sabre_seed;
    const int sabre_runs = 100;
    const int hybrid_sabre_seed_trials = 50; // liczba losowych ustawien startowych SABRE dla ziarna DE

    // 1. Wykonanie algorytmu bazowego (Benchmark zachłanny ASAP)
    // GreedyASAP greedy_router(qpu, circuit, qpu.get_num_qubits(), name);
    // greedy_router.run_benchmark();

    // 2. Wykonanie algorytmu optymalizacyjnego (Differential Evolution)
    // DifferentialEvolution de_router(qpu, circuit, qpu.get_num_qubits(), name);
    // de_router.optimize(100, 2000000);
    // de_router.optimize_windowed(100, 30000, 200, 1.0, qpu.get_num_qubits()); //30, 1000, 100, 0.35

    // 3. Wykonanie algorytmu SABRE
    SabreRouter sabre(qpu, circuit, qpu.get_num_qubits(), name);
    sabre.set_random_seed(use_random_seeds ? rd() : sabre_seed);
    sabre.run_benchmark_multiple(sabre_runs);

    // 4. Wykonanie algorytmu DEHybridSabre
    SabreDEHybridRouter SabreDERouter(qpu, circuit, qpu.get_num_qubits(), name, false, 3);
    SabreDERouter.set_random_seed(use_random_seeds ? rd() : hybrid_seed);
    SabreDERouter.set_sabre_seed_trials(hybrid_sabre_seed_trials);
    // SabreDERouter.optimize_windowed(30, 400, 500, 0.6, false, 0.8);
    SabreDERouter.optimize_with_feedback(30, 200, 10, false, 0.8, 0.9);
    // SabreDERouter.optimize_windowed_with_feedback(30, 300, 500, 0.7, 10, false, 0.8, 0.9);

    // 5. Wykonanie algorytmu newDEHybridSabre
    // newSabreDEHybridRouter newSabreDERouter(qpu, circuit, qpu.get_num_qubits(), name, 4);
    // newSabreDERouter.optimize_windowed(30, 100, 500, 0.7, true, false, 0.8, 0.9);

}

void run_test_for_topologyII(const BaseTopology &qpu, const std::string &name) {
    std::cout << "\n============================================\n";
    std::cout << "URUCHAMIAM TEST DLA TOPOLOGII: " << name << "\n";
    std::cout << "============================================\n";
    qpu.print_stats();

    int logical_qubits_count = qpu.get_num_qubits() * 1.0;
    CircuitGenerator generator(qpu, logical_qubits_count, 1234);

    std::vector<Gate> circuit = generator.generate(2000, 0.2, 0.1);
    generator.print_gate_statistics(circuit);
    generator.print_circuit(circuit);

    // --- Konfiguracja losowosci i liczby prob ---
    const bool use_random_seeds = true;
    std::random_device rd;
    const uint32_t sabre_seed = rd();
    const uint64_t hybrid_seed = sabre_seed;
    const int sabre_runs = 20;
    const int hybrid_runs = 5;
    const int hybrid_sabre_seed_trials = 20;

    SabreRouter sabre(qpu, circuit, qpu.get_num_qubits(), name);
    sabre.set_random_seed(use_random_seeds ? rd() : sabre_seed);
    sabre.run_benchmark_multiple(sabre_runs);

    std::cout << "\n--- Rozpoczynam " << hybrid_runs << " iteracji dla DEHybridSabre ---\n";

    std::vector<OptimizationResult> results;
    results.reserve(hybrid_runs);

    for (int i = 0; i < hybrid_runs; ++i) {
        SabreDEHybridRouter SabreDERouter(qpu, circuit, qpu.get_num_qubits(), name, false, 1);
        SabreDERouter.set_random_seed(use_random_seeds ? rd() : hybrid_seed + i);
        SabreDERouter.set_sabre_seed_trials(hybrid_sabre_seed_trials);

        // Odbieramy całą strukturę wyników
        OptimizationResult res = SabreDERouter.optimize_windowed(30, 300, 300, 0.6,  false, 0.9);
        // OptimizationResult res = SabreDERouter.optimize_windowed_with_feedback(30, 300, 500, 0.7, 6, false, 0.8, 0.9);
        results.push_back(res);

        std::cout << "Iteracja " << (i + 1) << " zakonczona. "
                  << "Wynik: " << res.score << ", "
                  << "SWAPy: " << res.swaps << ", "
                  << "Makespan: " << res.makespan << "\n";
    }

    // --- OBLICZANIE STATYSTYK ---
    OptimizationResult min_val = results[0];
    OptimizationResult max_val = results[0];
    OptimizationResult sum_val = results[0];
    OptimizationResult best_result = results[0]; // Zmienna do śledzenia pojedynczego najlepszego wyniku

    for (int i = 1; i < hybrid_runs; ++i) {
        // Zapisujemy całą instancję wyniku, jeśli ma lepszy (niższy) score
        if (results[i].score < best_result.score) {
            best_result = results[i];
        }

        min_val.score = std::min(min_val.score, results[i].score);
        max_val.score = std::max(max_val.score, results[i].score);
        sum_val.score += results[i].score;

        min_val.makespan = std::min(min_val.makespan, results[i].makespan);
        max_val.makespan = std::max(max_val.makespan, results[i].makespan);
        sum_val.makespan += results[i].makespan;

        min_val.depth = std::min(min_val.depth, results[i].depth);
        max_val.depth = std::max(max_val.depth, results[i].depth);
        sum_val.depth += results[i].depth;

        min_val.swaps = std::min(min_val.swaps, results[i].swaps);
        max_val.swaps = std::max(max_val.swaps, results[i].swaps);
        sum_val.swaps += results[i].swaps;

        min_val.moves = std::min(min_val.moves, results[i].moves);
        max_val.moves = std::max(max_val.moves, results[i].moves);
        sum_val.moves += results[i].moves;
    }

    std::cout << "\n============================================\n";
    std::cout << " PODSUMOWANIE DEHybridSabre (" << hybrid_runs << " uruchomien)\n";
    std::cout << "============================================\n";

    // --- WYPISYWANIE NAJLEPSZEGO POJEDYNCZEGO WYNIKU ---
    std::cout << " [NAJLEPSZY ZNALEZIONY WYNIK (wg min Score)]\n";
    std::cout << " Score: " << best_result.score << "\n";
    std::cout << " Makespan: " << best_result.makespan << "\n";
    std::cout << " Depth: " << best_result.depth << "\n";
    std::cout << " Swaps: " << best_result.swaps << "\n";
    std::cout << " Moves: " << best_result.moves << "\n";
    std::cout << "--------------------------------------------\n";

    std::cout << " [SCORE]\n";
    std::cout << " Min: " << min_val.score << " | Max: " << max_val.score << " | Średnia: " << (sum_val.score / hybrid_runs) << "\n";

    std::cout << " [MAKESPAN]\n";
    std::cout << " Min: " << min_val.makespan << " | Max: " << max_val.makespan << " | Średnia: " << (sum_val.makespan / hybrid_runs) << "\n";

    std::cout << " [DEPTH]\n";
    std::cout << " Min: " << min_val.depth << " | Max: " << max_val.depth << " | Średnia: " << (static_cast<double>(sum_val.depth) / hybrid_runs) << "\n";

    std::cout << " [SWAPS]\n";
    std::cout << " Min: " << min_val.swaps << " | Max: " << max_val.swaps << " | Średnia: " << (static_cast<double>(sum_val.swaps) / hybrid_runs) << "\n";

    std::cout << " [MOVES]\n";
    std::cout << " Min: " << min_val.moves << " | Max: " << max_val.moves << " | Średnia: " << (static_cast<double>(sum_val.moves) / hybrid_runs) << "\n";
    std::cout << "============================================\n";
}

 int main() {
 #ifdef _WIN32
     // Windows potrzebuje instrukcji zmiany kodowania na UTF-8
     SetConsoleOutputCP(CP_UTF8);
     // Opcjonalnie: upewnia się, że wejście (cin) też używa UTF-8
     SetConsoleCP(CP_UTF8);
 #endif
     // Topologia siatki 2x2
       IQMConstellation gridQpu(5, 5);
       // run_test_for_topology(gridQpu, "Siatka 5x5");
       run_test_for_topologyII(gridQpu, "Siatka 5x5");

     // Topologia trójkątna
//     IQMTriangleConstellation triangleQpu;
//     run_test_for_topology(triangleQpu, "Trójkąt 3-Rezonatory");
     // run_test_for_topologyII(triangleQpu, "Trójkąt 3-Rezonatory");

     return 0;
 }

// void run_test_for_topologyII(const BaseTopology &qpu, const std::string &name, int num_circuits, int runs_per_circuit) {
//     std::cout << "\n============================================\n";
//     std::cout << "URUCHAMIAM TEST DLA TOPOLOGII: " << name << "\n";
//     std::cout << "Testowanych obwodow (X): " << num_circuits << "\n";
//     std::cout << "Liczba przebiegow na obwod (Y): " << runs_per_circuit << "\n";
//     std::cout << "============================================\n";
//     qpu.print_stats();
//
//     int logical_qubits_count = qpu.get_num_qubits() * 1.0;
//     const bool use_random_seeds = true;
//     std::random_device rd;
//
//     // --- PĘTLA PO RÓŻNYCH OBWODACH (X) ---
//     for (int c = 0; c < num_circuits; ++c) {
//         std::cout << "\n********************************************\n";
//         std::cout << " GENEROWANIE I TESTOWANIE OBWODU " << (c + 1) << " / " << num_circuits << "\n";
//         std::cout << "********************************************\n";
//
//         // Zmiana seeda dla każdego obwodu, aby zawsze był inny
//         int circuit_seed = use_random_seeds ? rd() : (1234);
//         CircuitGenerator generator(qpu, logical_qubits_count, circuit_seed);
//
//         std::vector<Gate> circuit = generator.generate(2000, 0.6, 0.4);
//         generator.print_gate_statistics(circuit);
//
//         // Zostawiam zakomentowane, aby nie zaśmiecać konsoli przy wielu obwodach
//         // generator.print_circuit(circuit);
//
//         // --- Konfiguracja losowosci i liczby prob dla danego obwodu ---
//         const uint32_t sabre_seed = rd();
//         const uint64_t hybrid_seed = sabre_seed;
//         const int sabre_runs = 100;
//         const int hybrid_sabre_seed_trials = 10;
//
//         SabreRouter sabre(qpu, circuit, qpu.get_num_qubits(), name);
//         sabre.set_random_seed(use_random_seeds ? rd() : sabre_seed);
//         sabre.run_benchmark_multiple(sabre_runs);
//
//         std::cout << "\n--- Rozpoczynam " << runs_per_circuit << " iteracji dla DEHybridSabre (Obwod " << (c + 1) << ") ---\n";
//
//         std::vector<OptimizationResult> results;
//         results.reserve(runs_per_circuit);
//
//         // --- PĘTLA PO PRZEBIEGACH DLA DANEGO OBWODU (Y) ---
//         for (int i = 0; i < runs_per_circuit; ++i) {
//             SabreDEHybridRouter SabreDERouter(qpu, circuit, qpu.get_num_qubits(), name, false, 3);
//             SabreDERouter.set_random_seed(use_random_seeds ? rd() : hybrid_seed + i);
//             SabreDERouter.set_sabre_seed_trials(hybrid_sabre_seed_trials);
//
//             // Odbieramy całą strukturę wyników
//             OptimizationResult res = SabreDERouter.optimize_windowed(30, 300, 500, 0.7, false, 0.9);
//             results.push_back(res);
//
//             std::cout << " Iteracja " << (i + 1) << " zakonczona. "
//                       << "Wynik: " << res.score << ", "
//                       << "SWAPy: " << res.swaps << ", "
//                       << "Makespan: " << res.makespan << "\n";
//         }
//
//         // --- OBLICZANIE STATYSTYK DLA DANEGO OBWODU ---
//         OptimizationResult min_val = results[0];
//         OptimizationResult max_val = results[0];
//         OptimizationResult sum_val = results[0];
//         OptimizationResult best_result = results[0];
//
//         for (int i = 1; i < runs_per_circuit; ++i) {
//             if (results[i].score < best_result.score) {
//                 best_result = results[i];
//             }
//
//             min_val.score = std::min(min_val.score, results[i].score);
//             max_val.score = std::max(max_val.score, results[i].score);
//             sum_val.score += results[i].score;
//
//             min_val.makespan = std::min(min_val.makespan, results[i].makespan);
//             max_val.makespan = std::max(max_val.makespan, results[i].makespan);
//             sum_val.makespan += results[i].makespan;
//
//             min_val.depth = std::min(min_val.depth, results[i].depth);
//             max_val.depth = std::max(max_val.depth, results[i].depth);
//             sum_val.depth += results[i].depth;
//
//             min_val.swaps = std::min(min_val.swaps, results[i].swaps);
//             max_val.swaps = std::max(max_val.swaps, results[i].swaps);
//             sum_val.swaps += results[i].swaps;
//
//             min_val.moves = std::min(min_val.moves, results[i].moves);
//             max_val.moves = std::max(max_val.moves, results[i].moves);
//             sum_val.moves += results[i].moves;
//         }
//
//         std::cout << "\n--------------------------------------------\n";
//         std::cout << " PODSUMOWANIE OBWODU " << (c + 1) << " (" << runs_per_circuit << " uruchomien)\n";
//         std::cout << "--------------------------------------------\n";
//
//         std::cout << " [NAJLEPSZY ZNALEZIONY WYNIK (wg min Score)]\n";
//         std::cout << " Score: " << best_result.score << "\n";
//         std::cout << " Makespan: " << best_result.makespan << "\n";
//         std::cout << " Depth: " << best_result.depth << "\n";
//         std::cout << " Swaps: " << best_result.swaps << "\n";
//         std::cout << " Moves: " << best_result.moves << "\n";
//         std::cout << "--------------------------------------------\n";
//
//         std::cout << " [SCORE]     Min: " << min_val.score << " | Max: " << max_val.score << " | Srednia: " << (sum_val.score / runs_per_circuit) << "\n";
//         std::cout << " [MAKESPAN]  Min: " << min_val.makespan << " | Max: " << max_val.makespan << " | Srednia: " << (sum_val.makespan / runs_per_circuit) << "\n";
//         std::cout << " [DEPTH]     Min: " << min_val.depth << " | Max: " << max_val.depth << " | Srednia: " << (static_cast<double>(sum_val.depth) / runs_per_circuit) << "\n";
//         std::cout << " [SWAPS]     Min: " << min_val.swaps << " | Max: " << max_val.swaps << " | Srednia: " << (static_cast<double>(sum_val.swaps) / runs_per_circuit) << "\n";
//         std::cout << " [MOVES]     Min: " << min_val.moves << " | Max: " << max_val.moves << " | Srednia: " << (static_cast<double>(sum_val.moves) / runs_per_circuit) << "\n";
//         std::cout << "============================================\n";
//     }
// }
//
// int main() {
// #ifdef _WIN32
//     // Windows potrzebuje instrukcji zmiany kodowania na UTF-8
//     SetConsoleOutputCP(CP_UTF8);
//     // Opcjonalnie: upewnia się, że wejście (cin) też używa UTF-8
//     SetConsoleCP(CP_UTF8);
// #endif
//
//     // Definiujemy X obwodów i Y przebiegów
//     const int X_CIRCUITS = 5;
//     const int Y_RUNS = 5;
//
//     // Topologia siatki 2x2
//     IQMConstellation gridQpu(5, 5);
//     run_test_for_topologyII(gridQpu, "Siatka 5x5", X_CIRCUITS, Y_RUNS);
//
//     // Topologia trójkątna (zakomentowane jako opcja)
//     // IQMTriangleConstellation triangleQpu;
//     // run_test_for_topologyII(triangleQpu, "Trójkąt 3-Rezonatory", X_CIRCUITS, Y_RUNS);
//
//     return 0;
// }


//int main(int argc, char* argv[]) {
//    // Domyślny folder z wygenerowanymi benchmarkami
//    std::string base_dir = "../iqm_benchmarks";
//
//    if (argc >= 2) {
//        base_dir = argv[1];
//    }
//
//    if (!std::filesystem::exists(base_dir) || !std::filesystem::is_directory(base_dir)) {
//        std::cerr << "Blad: Katalog '" << base_dir << "' nie istnieje. Najpierw uruchom skrypt Pythona!\n";
//        return -1;
//    }
//
//    std::cout << "=> Rozpoczynam przetwarzanie plikow w katalogu: " << base_dir << "\n";
//
//    // --- Konfiguracja losowosci i liczby prob ---
//    bool use_random_seeds = false;
//    std::random_device rd;
//    uint32_t sabre_seed = rd();
//    uint64_t hybrid_seed = sabre_seed;
//    int sabre_runs = 300;
//    int hybrid_sabre_seed_trials = 4;
//
//    // Rekursywna iteracja po wszystkich plikach i podfolderach
//    for (const auto& entry : std::filesystem::recursive_directory_iterator(base_dir)) {
//        if (entry.is_regular_file() && entry.path().extension() == ".json") {
//
//            std::string input_file = entry.path().string();
//
//            // Zabezpieczenie: pomijamy pliki, które już są wynikiem naszego routingu
//            if (input_file.find("_routed.json") != std::string::npos) {
//                continue;
//            }
//
//            // Tworzymy nazwe pliku wyjsciowego: np. obwod_5q_50g_routed.json
//            std::string output_file = input_file.substr(0, input_file.find_last_of('.')) + "_routed.json";
//
//            std::cout << "\n--------------------------------------------------\n";
//            std::cout << "=> Wczytywanie: " << input_file << "\n";
//
//            int logical_qubits;
//            std::vector<Gate> circuit = CircuitIO::load_circuit(input_file, logical_qubits);
//
//            if(circuit.empty()) {
//                std::cerr << "Blad wczytywania obwodu. Pomijam plik.\n";
//                continue;
//            }
//
//            // =========================================================
//            // --- CAP DLA TESTÓW (DO 16 KUBITÓW) ---
//            // Zakomentuj poniższego if-a, gdy zechcesz odpalić całość
//            if (logical_qubits > 45 || logical_qubits > 44) {
//            // if (false) {
//                std::cout << "=> Pomijam (ma " << logical_qubits << " kubitow, a limit testowy to 16).\n";
//                continue; // Skaczemy do następnego pliku w folderze
//            }
//            // =========================================================
//
//            std::cout << "Bramki: " << circuit.size() << " | Kubity: " << logical_qubits << "\n";
//
//            // --- Logika Topologii ---
//            std::unique_ptr<BaseTopology> topology;
//            if (logical_qubits <= 13) {
//                topology = std::make_unique<IQMTriangleConstellation>();
//            } else {
//                int dim = 2;
//                topology = std::make_unique<IQMConstellation>(dim, dim);
//                while (topology->get_num_qubits() < logical_qubits) {
//                    dim++;
//                    topology = std::make_unique<IQMConstellation>(dim, dim);
//                }
//            }
//
//            // --- Konfiguracja Metaheurystyki ---
//            bool use_slow_mode = false;
//            int lookahead_depth = 3;
//
//            topology->print_stats();
//
//            SabreRouter sabre(*topology, circuit, logical_qubits, "sabre_run");
//            sabre.set_random_seed(use_random_seeds ? rd() : sabre_seed);
//            sabre.run_benchmark_multiple(sabre_runs);
//
//            SabreDEHybridRouter router(*topology, circuit, logical_qubits, "qiskit_run", use_slow_mode, lookahead_depth);
//            router.set_random_seed(use_random_seeds ? rd() : hybrid_seed);
//            router.set_sabre_seed_trials(hybrid_sabre_seed_trials);
//
//            int pop_size = 30;
//            int generations_per_window = 300;
//            int window_size = 16600;
//            double accept_ratio = 1.0;
//            bool lazy_offboard = true;
//
//            // --- Uruchomienie optymalizacji ---
//            OptimizationResult result = router.optimize_windowed(pop_size, generations_per_window, window_size, accept_ratio, lazy_offboard);
//
//            // --- Zapis wyjścia ---
//            // Odkomentuj i dostosuj ponizsze linijki zgodnie z Twoim API
//            // std::vector<Gate> routed_circuit = router.get_routed_circuit();
//            // CircuitIO::save_circuit(output_file, routed_circuit);
//
//            std::cout << "=> Zapisano wynik do: " << output_file << "\n";
//        }
//    }
//
//    std::cout << "\n==================================================\n";
//    std::cout << "Zakonczono przetwarzanie wszystkich obwodow z folderu: " << base_dir << "\n";
//    return 0;
//}

// Funkcja pomocnicza do zapisu wyników i kalkulacji statystyk do pliku
void save_statistics_to_file(const std::filesystem::path& file_path, const std::string& method_name, const std::vector<OptimizationResult>& results) {
    if (results.empty()) return;

    std::ofstream out(file_path);
    if (!out.is_open()) {
        std::cerr << "Blad: Nie mozna otworzyc pliku do zapisu: " << file_path << "\n";
        return;
    }

    out << "============================================\n";
    out << " RAPORT WYNIKOW: " << method_name << "\n";
    out << " Liczba uruchomien: " << results.size() << "\n";
    out << "============================================\n\n";

    // 1. Wyniki per iteracja
    out << "--- WYNIKI DLA KAZDEGO URUCHOMIENIA ---\n";
    for (size_t i = 0; i < results.size(); ++i) {
        out << "Run " << std::setw(3) << (i + 1)
            << " | Score: " << results[i].score
            << " | Makespan: " << results[i].makespan
            << " | Depth: " << results[i].depth
            << " | Swaps: " << results[i].swaps
            << " | Moves: " << results[i].moves << "\n";
    }

    // 2. Kalkulacja statystyk
    OptimizationResult min_val = results[0];
    OptimizationResult max_val = results[0];
    OptimizationResult sum_val = results[0];
    OptimizationResult best_result = results[0];

    for (size_t i = 1; i < results.size(); ++i) {
        if (results[i].score < best_result.score) {
            best_result = results[i];
        }

        min_val.score = std::min(min_val.score, results[i].score);
        max_val.score = std::max(max_val.score, results[i].score);
        sum_val.score += results[i].score;

        min_val.makespan = std::min(min_val.makespan, results[i].makespan);
        max_val.makespan = std::max(max_val.makespan, results[i].makespan);
        sum_val.makespan += results[i].makespan;

        min_val.depth = std::min(min_val.depth, results[i].depth);
        max_val.depth = std::max(max_val.depth, results[i].depth);
        sum_val.depth += results[i].depth;

        min_val.swaps = std::min(min_val.swaps, results[i].swaps);
        max_val.swaps = std::max(max_val.swaps, results[i].swaps);
        sum_val.swaps += results[i].swaps;

        min_val.moves = std::min(min_val.moves, results[i].moves);
        max_val.moves = std::max(max_val.moves, results[i].moves);
        sum_val.moves += results[i].moves;
    }

    double n = static_cast<double>(results.size());

    // 3. Podsumowanie zgodne z formatem
    out << "\n--- PODSUMOWANIE STATYSTYK ---\n";

    out << " [SCORE]\n";
    out << " Min: " << min_val.score << " | Max: " << max_val.score << " | Średnia: " << (sum_val.score / n) << "\n";

    out << " [MAKESPAN]\n";
    out << " Min: " << min_val.makespan << " | Max: " << max_val.makespan << " | Średnia: " << (sum_val.makespan / n) << "\n";

    out << " [DEPTH]\n";
    out << " Min: " << min_val.depth << " | Max: " << max_val.depth << " | Średnia: " << (sum_val.depth / n) << "\n";

    out << " [SWAPS]\n";
    out << " Min: " << min_val.swaps << " | Max: " << max_val.swaps << " | Średnia: " << (sum_val.swaps / n) << "\n";

    out << " [MOVES]\n";
    out << " Min: " << min_val.moves << " | Max: " << max_val.moves << " | Średnia: " << (sum_val.moves / n) << "\n";

    // 4. Detale najlepszego przebiegu
    out << "\n--- NAJLEPSZY ZNALEZIONY WYNIK (wg najniższego Score) ---\n";
    out << " Score: " << best_result.score << "\n";
    out << " Makespan: " << best_result.makespan << "\n";
    out << " Depth: " << best_result.depth << "\n";
    out << " Swaps: " << best_result.swaps << "\n";
    out << " Moves: " << best_result.moves << "\n";

    out.close();
}

// int main(int argc, char* argv[]) {
// #ifdef _WIN32
//     SetConsoleOutputCP(CP_UTF8);
//     SetConsoleCP(CP_UTF8);
// #endif
//
//     std::string base_dir = "../iqm_benchmarks";
//     std::string results_dir = "../results";
//
//     if (argc >= 2) base_dir = argv[1];
//     if (argc >= 3) results_dir = argv[2];
//
//     if (!std::filesystem::exists(base_dir) || !std::filesystem::is_directory(base_dir)) {
//         std::cerr << "Blad: Katalog '" << base_dir << "' nie istnieje. Najpierw uruchom skrypt Pythona!\n";
//         return -1;
//     }
//
//     std::cout << "=> Rozpoczynam przetwarzanie plikow w katalogu: " << base_dir << "\n";
//     std::cout << "=> Zapis wynikow do katalogu: " << results_dir << "\n";
//
//     // --- Konfiguracja losowosci i liczby prob ---
//     bool use_random_seeds = false;
//     std::random_device rd;
//     uint32_t sabre_seed = rd();
//     uint64_t hybrid_seed = sabre_seed;
//
//     int sabre_runs = 150;     // Liczba X (dla samego Sabre)
//     int hybrid_runs = 10;     // Liczba Y (dla hybrydy)
//     int hybrid_sabre_seed_trials = 30;
//
//     for (const auto& entry : std::filesystem::recursive_directory_iterator(base_dir)) {
//         if (entry.is_regular_file() && entry.path().extension() == ".json") {
//
//             std::string input_file = entry.path().string();
//
//             // --- Odwzorowanie struktury katalogów na ../results ---
//             std::filesystem::path in_path = entry.path();
//             std::filesystem::path rel_path = std::filesystem::relative(in_path, base_dir);
//             std::filesystem::path out_file_base = std::filesystem::path(results_dir) / rel_path;
//
//             // Stwórz drzewo folderów, jeśli nie istnieje
//             std::filesystem::create_directories(out_file_base.parent_path());
//             std::string stem = out_file_base.stem().string(); // np. "obwod_5q"
//
//             // =========================================================
//             // SCIEŻKA DLA PLIKÓW JUŻ ZROUTOWANYCH PRZEZ QISKIT
//             // =========================================================
//             if (input_file.find("_qiskit_routed.json") != std::string::npos) {
//                 std::cout << "\n--------------------------------------------------\n";
//                 std::cout << "=> Wczytywanie zroutowanego obwodu Qiskit: " << input_file << "\n";
//
//                 int num_qubits_qiskit = 0;
//                 auto qiskit_circuit = CircuitIO::load_circuitII(input_file, num_qubits_qiskit);
//
//                 if (qiskit_circuit.empty()) {
//                     std::cerr << "Blad wczytywania obwodu Qiskit. Pomijam.\n";
//                     continue;
//                 }
//
//                 std::unique_ptr<BaseTopology> topology;
//                 if (num_qubits_qiskit <= 13) {
//                     topology = std::make_unique<IQMTriangleConstellation>();
//                 } else {
//                     int dim = 2;
//                     topology = std::make_unique<IQMConstellation>(dim, dim);
//                     while (topology->get_num_qubits() < num_qubits_qiskit) {
//                         dim++;
//                         topology = std::make_unique<IQMConstellation>(dim, dim);
//                     }
//                 }
//
//                 bool use_slow_mode = false;
//                 int lookahead_depth = 4;
//                 SabreDEHybridRouter evaluator(*topology, qiskit_circuit, num_qubits_qiskit, "qiskit_run", use_slow_mode, lookahead_depth);
//                 evaluator.set_random_seed(use_random_seeds ? rd() : hybrid_seed);
//                 evaluator.set_sabre_seed_trials(hybrid_sabre_seed_trials);
//
//                 // Uruchomienie ewaluacji dla zroutowanego pliku Qiskit
//                 OptimizationResult qiskit_real_cost = evaluator.evaluate_qiskit_routed_circuit(qiskit_circuit);
//
//                 // Zapis wyników ewaluacji
//                 std::filesystem::path qiskit_out_path = out_file_base.parent_path() / (stem + "_eval.txt");
//                 save_statistics_to_file(qiskit_out_path, "Ewaluacja Qiskit", {qiskit_real_cost});
//                 std::cout << "=> Zapisano wyniki ewaluatora Qiskit do: " << qiskit_out_path << "\n";
//
//                 continue; // Przejdź do następnego pliku
//             }
//
//             // Pomijamy pliki, które router sam wyrzucał wcześniej w folderze źródłowym (Zabezpieczenie)
//             if (input_file.find("_routed.json") != std::string::npos) {
//                 continue;
//             }
//
//             // =========================================================
//             // STANDARDOWA ŚCIEŻKA (SABRE ORAZ HYBRYDA)
//             // =========================================================
//             std::cout << "\n--------------------------------------------------\n";
//             std::cout << "=> Wczytywanie bazowego obwodu: " << input_file << "\n";
//
//             int logical_qubits;
//             std::vector<Gate> circuit = CircuitIO::load_circuitII(input_file, logical_qubits);
//
//             if (circuit.empty()) {
//                 std::cerr << "Blad wczytywania obwodu. Pomijam plik.\n";
//                 continue;
//             }
//
//             std::cout << "Bramki: " << circuit.size() << " | Kubity: " << logical_qubits << "\n";
//
//             std::unique_ptr<BaseTopology> topology;
//             if (logical_qubits <= 13) {
//                 topology = std::make_unique<IQMTriangleConstellation>();
//             } else {
//                 int dim = 2;
//                 topology = std::make_unique<IQMConstellation>(dim, dim);
//                 while (topology->get_num_qubits() < logical_qubits) {
//                     dim++;
//                     topology = std::make_unique<IQMConstellation>(dim, dim);
//                 }
//             }
//
//             bool use_slow_mode = false;
//             int lookahead_depth = 4;
//
//             // ---------------------------------------------------------
//             // 1. URUCHOMIENIA BAZOWEGO SABRE (x uruchomień)
//             // ---------------------------------------------------------
//             std::cout << "=> Uruchamianie algorytmu Sabre (" << sabre_runs << " powtorzen)...\n";
//             std::vector<OptimizationResult> sabre_results;
//             sabre_results.reserve(sabre_runs);
//
//             for (int i = 0; i < sabre_runs; ++i) {
//                 SabreRouter sabre(*topology, circuit, logical_qubits, "sabre_run");
//                 sabre.set_random_seed(use_random_seeds ? rd() : sabre_seed + i);
//
//                 // Przypisujemy wynik z Sabre (SabreBenchmarkResult)
//                 SabreBenchmarkResult s_res = sabre.run_benchmark();
//
//                 // Konwertujemy na wspólną strukturę OptimizationResult dla wygody zapisu
//                 OptimizationResult res;
//                 res.score = s_res.score;
//                 res.makespan = s_res.makespan;
//                 res.depth = s_res.depth;
//                 res.swaps = s_res.swaps;
//                 res.moves = s_res.moves;
//
//                 sabre_results.push_back(res);
//             }
//             std::filesystem::path sabre_out_path = out_file_base.parent_path() / (stem + "_sabre.txt");
//             save_statistics_to_file(sabre_out_path, "Sabre Standard", sabre_results);
//             std::cout << "   Zapisano do: " << sabre_out_path << "\n";
//
//             // ---------------------------------------------------------
//             // 2. URUCHOMIENIA HYBRYDY DE+SABRE (y uruchomień)
//             // ---------------------------------------------------------
//             int pop_size = 30;
//             int generations_per_window = 300;
//             int window_size = 16000;
//             double accept_ratio = 1.0;
//             bool lazy_offboard = true;
//
//             std::cout << "=> Uruchamianie algorytmu DEHybridSabre (" << hybrid_runs << " powtorzen)...\n";
//             std::vector<OptimizationResult> hybrid_results;
//             hybrid_results.reserve(hybrid_runs);
//
//             for (int i = 0; i < hybrid_runs; ++i) {
//                 SabreDEHybridRouter router(*topology, circuit, logical_qubits, "sabre_run", use_slow_mode, lookahead_depth);
//                 router.set_random_seed(use_random_seeds ? rd() : hybrid_seed + i);
//                 router.set_sabre_seed_trials(hybrid_sabre_seed_trials);
//
//                 OptimizationResult res = router.optimize_windowed(pop_size, generations_per_window, window_size, accept_ratio, lazy_offboard);
//                 hybrid_results.push_back(res);
//             }
//             std::filesystem::path hybrid_out_path = out_file_base.parent_path() / (stem + "_hybrid.txt");
//             save_statistics_to_file(hybrid_out_path, "DEHybrid Sabre", hybrid_results);
//             std::cout << "   Zapisano do: " << hybrid_out_path << "\n";
//
//         }
//     }
//
//     std::cout << "\n==================================================\n";
//     std::cout << "Zakonczono przetwarzanie wszystkich obwodow z folderu: " << base_dir << "\n";
//     return 0;
// }