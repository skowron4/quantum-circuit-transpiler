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
#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#ifdef _WIN32
#include <windows.h>
#endif
#include <filesystem>
#include <omp.h>
#include "greedy/GreedyASAP.h"

struct CircuitFileMetadata {
    int logical_qubits = 0;
    int physical_qubits = 0;
    std::vector<int> qiskit_initial_layout_l2p;
    std::vector<int> qiskit_final_layout_l2p;
};

class CircuitIO {
public:
    static std::vector<Gate> load_circuit(const std::string& filename, int& out_num_qubits) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Blad otwarcia pliku: " << filename << "\n";
            return {};
        }

        nlohmann::json j; file >> j;
        out_num_qubits = j.value("logical_qubits", j["num_qubits"].get<int>());
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

        out_num_qubits = j.value("logical_qubits", j["num_qubits"].get<int>());
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

    static CircuitFileMetadata load_metadata(const std::string& filename) {
        CircuitFileMetadata metadata;
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Blad otwarcia pliku metadanych: " << filename << "\n";
            return metadata;
        }

        nlohmann::json j;
        try {
            file >> j;
        } catch (const nlohmann::json::parse_error& e) {
            std::cerr << "Blad parsowania metadanych JSON w pliku " << filename << ": " << e.what() << "\n";
            return metadata;
        }

        metadata.logical_qubits = j.value("logical_qubits", j.value("num_qubits", 0));
        metadata.physical_qubits = j.value("physical_qubits", j.value("num_qubits", metadata.logical_qubits));

        if (j.contains("qiskit_initial_layout_l2p") && j["qiskit_initial_layout_l2p"].is_array()) {
            metadata.qiskit_initial_layout_l2p = j["qiskit_initial_layout_l2p"].get<std::vector<int>>();
        }
        if (j.contains("qiskit_final_layout_l2p") && j["qiskit_final_layout_l2p"].is_array()) {
            metadata.qiskit_final_layout_l2p = j["qiskit_final_layout_l2p"].get<std::vector<int>>();
        }

        return metadata;
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
    const int sabre_runs = 1;
    const int hybrid_sabre_seed_trials = 1; // liczba losowych ustawien startowych SABRE dla ziarna DE
    const int sabre_lookahead_depth = 50;
    const int hybrid_lookahead_depth = 4;
    const int hybrid_pop_size = 30;
    const int hybrid_generations_per_window = 300;
    const int hybrid_window_size = 300;
    const double hybrid_accept_ratio = 0.6;
    const int hybrid_feedback_passes = 3;
    const bool hybrid_allow_lazy_offboard = false;
    const double hybrid_f = 0.8;
    const double hybrid_cr = 0.9;

    // 1. Wykonanie algorytmu bazowego (Benchmark zachłanny ASAP)
    // GreedyASAP greedy_router(qpu, circuit, qpu.get_num_qubits(), name);
    // greedy_router.run_benchmark();

    // 2. Wykonanie algorytmu optymalizacyjnego (Differential Evolution)
    // DifferentialEvolution de_router(qpu, circuit, qpu.get_num_qubits(), name);
    // de_router.optimize(100, 2000000);
    // de_router.optimize_windowed(100, 30000, 200, 1.0, qpu.get_num_qubits()); //30, 1000, 100, 0.35

    // 3. Wykonanie algorytmu SABRE
    // SabreRouter sabre(qpu, circuit, qpu.get_num_qubits(), name, sabre_lookahead_depth);
    // sabre.set_random_seed(use_random_seeds ? rd() : sabre_seed);
    // sabre.run_benchmark_multiple(sabre_runs);

    // 4. Wykonanie algorytmu DEHybridSabre
    SabreDEHybridRouter SabreDERouter(qpu, circuit, qpu.get_num_qubits(), name, false, hybrid_lookahead_depth);
    SabreDERouter.set_random_seed(use_random_seeds ? rd() : hybrid_seed);
    SabreDERouter.set_sabre_seed_trials(hybrid_sabre_seed_trials);
    SabreDERouter.optimize_windowed(30, 300, 500, 0.7, false, 0.9);
//    SabreDERouter.optimize_windowed_with_feedback(
//        hybrid_pop_size,
//        hybrid_generations_per_window,
//        hybrid_window_size,
//        hybrid_accept_ratio,
//        hybrid_feedback_passes,
//        hybrid_allow_lazy_offboard,
//        hybrid_f,
//        hybrid_cr
//    );

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

    std::vector<Gate> circuit = generator.generate(2000, 0.6, 0.4);
    generator.print_gate_statistics(circuit);
    generator.print_circuit(circuit);

    // --- Konfiguracja losowosci i liczby prob ---
    const bool use_random_seeds = true;
    std::random_device rd;
    const uint32_t sabre_seed = rd();
    const uint64_t hybrid_seed = sabre_seed;
    const int sabre_runs = 1;
    const int hybrid_runs = 5;
    const int hybrid_sabre_seed_trials = 1;
    const int sabre_lookahead_depth = 30;
    const int hybrid_lookahead_depth = 10;
    const int hybrid_pop_size = 30;
    const int hybrid_generations_per_window = 300;
    const int hybrid_window_size = 600;
    const double hybrid_accept_ratio = 0.35;
    const int hybrid_feedback_passes = 3;
    const bool hybrid_allow_lazy_offboard = false;
    const double hybrid_f = 0.8;
    const double hybrid_cr = 0.9;

    SabreRouter sabre(qpu, circuit, qpu.get_num_qubits(), name, sabre_lookahead_depth);
    sabre.set_random_seed(use_random_seeds ? rd() : sabre_seed);
    sabre.run_benchmark_multiple(sabre_runs);

    std::cout << "\n--- Rozpoczynam " << hybrid_runs << " iteracji dla DEHybridSabre ---\n";

    std::vector<OptimizationResult> results;
    results.reserve(hybrid_runs);

    for (int i = 0; i < hybrid_runs; ++i) {
        SabreDEHybridRouter SabreDERouter(qpu, circuit, qpu.get_num_qubits(), name, false, hybrid_lookahead_depth);
        SabreDERouter.set_random_seed(use_random_seeds ? rd() : hybrid_seed + i);
        SabreDERouter.set_sabre_seed_trials(hybrid_sabre_seed_trials);
        OptimizationResult res = SabreDERouter.optimize_windowed(30, 100, 500, 0.6, false, 0.2);

        // Odbieramy całą strukturę wyników
//        OptimizationResult res = SabreDERouter.optimize_windowed_with_feedback(
//            hybrid_pop_size,
//            hybrid_generations_per_window,
//            hybrid_window_size,
//            hybrid_accept_ratio,
//            hybrid_feedback_passes,
//            hybrid_allow_lazy_offboard,
//            hybrid_f,
//            hybrid_cr
//        );
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
    run_test_for_topology(gridQpu, "Siatka 5x5");
    // run_test_for_topologyII(gridQpu, "Siatka 5x5");

    // Topologia trójkątna
    // IQMTriangleConstellation triangleQpu;
    // run_test_for_topology(triangleQpu, "Trójkąt 3-Rezonatory");
    // run_test_for_topologyII(triangleQpu, "Trójkąt 3-Rezonatory");

    return 0;
}

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

struct HybridRunParams {
    int pop_size = 30;
    int generations_per_window = 180;
    int window_size = 500;
    double accept_ratio = 0.7;
    bool lazy_offboard = false;
    double f = 0.8;
    double cr = 0.9;
    int lookahead_depth = 4;
    std::string profile = "medium";
};

struct ExperimentConfig {
    std::filesystem::path input_dir = "../iqm_benchmarks";
    std::filesystem::path results_dir = "../results";
    int deterministic_runs = 10;
    int random_pilot_layouts = 8;
    int random_tuning_instances_per_bucket = 2;
    int min_random_layouts = 5;
    int max_random_layouts = 20;
    double random_target_relative_error = 0.10;
    int hyperparam_training_instances_per_bucket = 3;
    int hyperparam_training_runs = 2;
    double hyperparam_time_penalty = 0.02;
    uint64_t base_seed = 1337;
};

struct CircuitMeta {
    std::filesystem::path path;
    std::string circuit_id;
    std::string topology_bucket;
    std::string size_bucket;
    std::string training_bucket;
    int logical_qubits = 0;
    int gate_count = 0;
};

struct LayoutCandidate {
    std::string strategy;
    int layout_id = 0;
    uint64_t seed = 0;
    std::vector<int> l2p;
    double connectivity_score = 0.0;
};

struct CsvWriter {
    std::ofstream out;

    CsvWriter(const std::filesystem::path& path, const std::string& header) {
        const bool exists = std::filesystem::exists(path) && std::filesystem::file_size(path) > 0;
        std::filesystem::create_directories(path.parent_path());
        out.open(path, std::ios::app);
        if (!exists) out << header << "\n";
    }

    // Nowa funkcja: bezpieczny wrzut dużej paczki danych na raz
    void write_chunk(const std::string& chunk) {
        if (chunk.empty()) return;
        #pragma omp critical(csv_write)
        {
            out << chunk;
            // opcjonalny flush pomaga przy dużych równoległych testach
            // out.flush();
        }
    }

    static std::string escape(const std::string& value) {
        if (value.find_first_of(",\"\n\r") == std::string::npos) return value;
        std::string escaped = "\"";
        for (char c : value) {
            if (c == '"') escaped += "\"\"";
            else escaped += c;
        }
        escaped += "\"";
        return escaped;
    }
};

static std::string mapping_to_string(const std::vector<int>& mapping) {
    std::ostringstream ss;
    for (size_t i = 0; i < mapping.size(); ++i) {
        if (i > 0) ss << ' ';
        ss << i << ':' << mapping[i];
    }
    return ss.str();
}

static bool is_routed_qiskit_file(const std::filesystem::path& path) {
    return path.stem().string().find("_qiskit_routed") != std::string::npos;
}

static std::string circuit_family_from_path(const std::filesystem::path& path) {
    std::string stem = path.stem().string();
    const std::string marker = "_qiskit_routed";
    size_t marker_pos = stem.find(marker);
    if (marker_pos != std::string::npos) stem.erase(marker_pos, marker.size());
    return stem;
}

static std::unique_ptr<BaseTopology> make_topology_for_qubits(int logical_qubits, std::string& topology_name) {
    if (logical_qubits <= 13) {
        topology_name = "Trojkat 3-Rezonatory";
        return std::make_unique<IQMTriangleConstellation>();
    }

    int dim = 2;
    auto topology = std::make_unique<IQMConstellation>(dim, dim);
    while (topology->get_num_qubits() < logical_qubits) {
        ++dim;
        topology = std::make_unique<IQMConstellation>(dim, dim);
    }
    topology_name = "Siatka " + std::to_string(dim) + "x" + std::to_string(dim);
    return topology;
}

static std::string size_bucket_for_gates(int gate_count) {
    if (gate_count <= 500) return "gates_le_500";
    if (gate_count <= 2000) return "gates_le_2000";
    if (gate_count <= 8000) return "gates_le_8000";
    return "gates_gt_8000";
}

static int constellation_dimension_for_qubits(int logical_qubits) {
    int dim = 2;
    auto topology = std::make_unique<IQMConstellation>(dim, dim);
    while (topology->get_num_qubits() < logical_qubits) {
        ++dim;
        topology = std::make_unique<IQMConstellation>(dim, dim);
    }
    return dim;
}

static std::string topology_bucket_for_qubits(int logical_qubits) {
    if (logical_qubits <= 13) return "triangle";

    int dim = constellation_dimension_for_qubits(logical_qubits);
    if (dim <= 2) return "grid_2x2";
    if (dim <= 4) return "grid_3x3_4x4";
    return "grid_5x5_plus";
}

static std::string training_bucket_for(int logical_qubits, int gate_count) {
    return topology_bucket_for_qubits(logical_qubits) + "__" + size_bucket_for_gates(gate_count);
}

static HybridRunParams fallback_hybrid_params(int gate_count, int logical_qubits) {
    HybridRunParams params;
    params.profile = "fallback_" + training_bucket_for(logical_qubits, gate_count);
    params.lazy_offboard = true;
    params.window_size = std::clamp(gate_count / 2, 50, 16000);
    params.window_size = std::max(50, std::min(params.window_size, gate_count));
    return params;
}

static std::vector<int> window_candidates_for_gate_count(int gate_count) {
    std::vector<int> base = {50, 100, 250, 500, 1000, 2000, 4000, 8000, 16000};
    std::vector<int> candidates;
    for (int value : base) {
        if (value <= gate_count) candidates.push_back(value);
    }
    candidates.push_back(std::max(50, std::min(gate_count, gate_count / 2)));
    candidates.push_back(std::max(50, std::min(gate_count, gate_count)));
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    return candidates;
}

static std::vector<HybridRunParams> build_hyperparameter_candidates(int gate_count, int logical_qubits, int physical_qubits) {
    std::vector<HybridRunParams> candidates;
    const bool lazy = logical_qubits < physical_qubits;

    std::vector<int> pops = {10, 20, 30, 40, 50};
    std::vector<int> generations = {100, 250, 500, 1000};
    std::vector<double> fs = {0.1, 0.3, 0.5, 0.7, 0.9};
    std::vector<double> crs = {0.1, 0.3, 0.5, 0.7, 0.9};
    std::vector<double> accepts = {0.35, 0.5, 0.7, 0.9, 1.0};
    std::vector<int> depths = {0, 2, 4, 6, 8, 10};
    std::vector<int> windows = window_candidates_for_gate_count(gate_count);

    auto add_candidate = [&](int idx, int pop, int gen, int window, double accept, double f, double cr, int depth, const std::string& tag) {
        HybridRunParams params;
        params.pop_size = pop;
        params.generations_per_window = gen;
        params.window_size = std::max(50, std::min(window, gate_count));
        params.accept_ratio = accept;
        params.lazy_offboard = lazy;
        params.f = f;
        params.cr = cr;
        params.lookahead_depth = depth;
        params.profile = training_bucket_for(logical_qubits, gate_count) + "_" + tag + "_" + std::to_string(idx);
        candidates.push_back(params);
    };

    int idx = 0;
    for (size_t i = 0; i < pops.size(); ++i) {
        int window = windows[std::min(i, windows.size() - 1)];
        add_candidate(idx++, pops[i], generations[i % generations.size()], window,
                      accepts[(i + 2) % accepts.size()], fs[i], crs[(i + 1) % crs.size()],
                      depths[(i + 1) % depths.size()], "pop_sweep");
    }

    for (size_t i = 0; i < generations.size(); ++i) {
        int window = windows[std::min(windows.size() - 1, i + 1)];
        add_candidate(idx++, 30, generations[i], window, accepts[i % accepts.size()],
                      fs[(i + 2) % fs.size()], crs[(i + 3) % crs.size()],
                      depths[(i + 2) % depths.size()], "gen_sweep");
    }

    for (size_t i = 0; i < windows.size(); ++i) {
        add_candidate(idx++, 30, gate_count <= 500 ? 500 : 250, windows[i],
                      accepts[i % accepts.size()], fs[(i + 1) % fs.size()],
                      crs[(i + 2) % crs.size()], depths[i % depths.size()], "window_sweep");
    }

    for (size_t i = 0; i < accepts.size(); ++i) {
        int window = windows[windows.size() / 2];
        add_candidate(idx++, 40, gate_count <= 500 ? 500 : 250, window, accepts[i],
                      fs[i], crs[crs.size() - 1 - i], depths[(i + 3) % depths.size()],
                      "control_sweep");
    }

    return candidates;
}

static bool validate_mapping(const std::vector<int>& mapping, int logical_qubits, int physical_qubits) {
    if (static_cast<int>(mapping.size()) < logical_qubits) return false;
    std::set<int> used;
    for (int lq = 0; lq < logical_qubits; ++lq) {
        int p = mapping[lq];
        if (p < 0 || p >= physical_qubits || used.contains(p)) return false;
        used.insert(p);
    }
    return true;
}

static std::vector<int> make_identity_layout(int logical_qubits, int physical_qubits) {
    std::vector<int> mapping(logical_qubits, -1);
    for (int lq = 0; lq < logical_qubits && lq < physical_qubits; ++lq) mapping[lq] = lq;
    return mapping;
}

static std::vector<int> make_random_layout(int logical_qubits, int physical_qubits, std::mt19937_64& rng) {
    std::vector<int> physical(physical_qubits);
    std::iota(physical.begin(), physical.end(), 0);
    std::shuffle(physical.begin(), physical.end(), rng);

    std::vector<int> mapping(logical_qubits, -1);
    for (int lq = 0; lq < logical_qubits; ++lq) mapping[lq] = physical[lq];
    return mapping;
}

static std::vector<std::vector<double>> build_interaction_matrix(const std::vector<Gate>& circuit, int logical_qubits) {
    std::vector<std::vector<double>> weights(logical_qubits, std::vector<double>(logical_qubits, 0.0));
    for (size_t i = 0; i < circuit.size(); ++i) {
        const Gate& g = circuit[i];
        if (g.type == "1Q" || g.target_q2 < 0) continue;
        const double time_weight = 1.0 / (1.0 + static_cast<double>(i) / std::max<size_t>(1, circuit.size()));
        weights[g.target_q1][g.target_q2] += time_weight;
        weights[g.target_q2][g.target_q1] += time_weight;
    }
    return weights;
}

static std::vector<int> physical_connectivity_degrees(const BaseTopology& topology) {
    std::vector<int> degrees(topology.qubits.size(), 0);
    for (int p = 0; p < static_cast<int>(topology.qubits.size()); ++p) {
        std::set<int> neighbors;
        for (int r : topology.qubits[p].connected_resonators) {
            for (int nb : topology.resonators[r].connected_qubits) {
                if (nb != p) neighbors.insert(nb);
            }
        }
        degrees[p] = static_cast<int>(neighbors.size());
    }
    return degrees;
}

static std::vector<int> make_compact_layout(const BaseTopology& topology, const std::vector<Gate>& circuit, int logical_qubits) {
    const int physical_qubits = topology.get_num_qubits();
    std::vector<int> mapping(logical_qubits, -1);
    std::vector<bool> used_physical(physical_qubits, false);
    std::vector<bool> mapped_logical(logical_qubits, false);

    auto weights = build_interaction_matrix(circuit, logical_qubits);
    auto degrees = physical_connectivity_degrees(topology);

    int first_lq = 0;
    double best_logical_weight = -1.0;
    for (int lq = 0; lq < logical_qubits; ++lq) {
        double sum = std::accumulate(weights[lq].begin(), weights[lq].end(), 0.0);
        if (sum > best_logical_weight) {
            best_logical_weight = sum;
            first_lq = lq;
        }
    }

    int first_p = static_cast<int>(std::max_element(degrees.begin(), degrees.end()) - degrees.begin());
    mapping[first_lq] = first_p;
    mapped_logical[first_lq] = true;
    used_physical[first_p] = true;

    for (int mapped_count = 1; mapped_count < logical_qubits; ++mapped_count) {
        int next_lq = -1;
        double next_lq_score = -1.0;
        for (int lq = 0; lq < logical_qubits; ++lq) {
            if (mapped_logical[lq]) continue;
            double score = 0.0;
            for (int other = 0; other < logical_qubits; ++other) {
                if (mapped_logical[other]) score += weights[lq][other];
            }
            if (score > next_lq_score) {
                next_lq_score = score;
                next_lq = lq;
            }
        }
        if (next_lq == -1) break;

        int best_p = -1;
        double best_cost = std::numeric_limits<double>::max();
        for (int p = 0; p < physical_qubits; ++p) {
            if (used_physical[p]) continue;
            double weighted_distance = 0.0;
            for (int other = 0; other < logical_qubits; ++other) {
                if (mapped_logical[other] && weights[next_lq][other] > 0.0) {
                    weighted_distance += weights[next_lq][other] * topology.get_routing_path(p, mapping[other]).resonators.size();
                }
            }
            double compact_tie_break = static_cast<double>(degrees[p]) * 0.01;
            double cost = weighted_distance - compact_tie_break;
            if (cost < best_cost) {
                best_cost = cost;
                best_p = p;
            }
        }

        mapping[next_lq] = best_p;
        mapped_logical[next_lq] = true;
        used_physical[best_p] = true;
    }

    return mapping;
}

static double layout_connectivity_score(const BaseTopology& topology, const std::vector<Gate>& circuit, const std::vector<int>& mapping, int logical_qubits) {
    auto weights = build_interaction_matrix(circuit, logical_qubits);
    double score = 0.0;
    for (int i = 0; i < logical_qubits; ++i) {
        for (int j = i + 1; j < logical_qubits; ++j) {
            if (weights[i][j] <= 0.0) continue;
            int distance = static_cast<int>(topology.get_routing_path(mapping[i], mapping[j]).resonators.size());
            score += weights[i][j] / (1.0 + distance);
        }
    }
    return score;
}

static double estimate_random_layout_count_raw(const std::vector<double>& pilot_scores, const ExperimentConfig& config) {
    if (pilot_scores.size() < 2) return static_cast<double>(config.min_random_layouts);

    double mean = std::accumulate(pilot_scores.begin(), pilot_scores.end(), 0.0) / static_cast<double>(pilot_scores.size());
    double variance = 0.0;
    for (double value : pilot_scores) variance += (value - mean) * (value - mean);
    variance /= static_cast<double>(pilot_scores.size() - 1);
    double sd = std::sqrt(variance);

    if (std::abs(mean) < 1e-9) return static_cast<double>(config.max_random_layouts);

    return std::pow(1.96 * sd / (config.random_target_relative_error * std::abs(mean)), 2.0);
}

static int estimate_random_layout_count(const std::vector<double>& pilot_scores, const ExperimentConfig& config) {
    double required = estimate_random_layout_count_raw(pilot_scores, config);
    int count = static_cast<int>(std::ceil(required));
    return std::clamp(count, config.min_random_layouts, config.max_random_layouts);
}

static std::filesystem::path qiskit_routed_path_for_base(const std::filesystem::path& base_path) {
    return base_path.parent_path() / (base_path.stem().string() + "_qiskit_routed" + base_path.extension().string());
}

static std::vector<int> derive_qiskit_layout_from_routed_file(const std::filesystem::path& routed_path, int logical_qubits, int physical_qubits) {
    std::ifstream file(routed_path);
    if (file.is_open()) {
        try {
            nlohmann::json j;
            file >> j;
            if (j.contains("qiskit_initial_layout_l2p") && j["qiskit_initial_layout_l2p"].is_array()) {
                std::vector<int> layout(logical_qubits, -1);
                const auto& json_layout = j["qiskit_initial_layout_l2p"];
                for (int lq = 0; lq < logical_qubits && lq < static_cast<int>(json_layout.size()); ++lq) {
                    layout[lq] = json_layout[lq].get<int>();
                }
                if (validate_mapping(layout, logical_qubits, physical_qubits)) {
                    return layout;
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "Ostrzezenie: nie udalo sie odczytac initial_layout Qiskita z " << routed_path
                      << ": " << e.what() << "\n";
        }
    }

    int routed_qubits = 0;
    std::vector<Gate> routed_circuit = CircuitIO::load_circuitII(routed_path.string(), routed_qubits);
    std::vector<int> p2l(physical_qubits, -1);
    for (int p = 0; p < physical_qubits && p < logical_qubits; ++p) {
        p2l[p] = p;
    }

    for (const auto& gate : routed_circuit) {
        if (gate.type == "SWAP" &&
            gate.target_q1 >= 0 && gate.target_q1 < physical_qubits &&
            gate.target_q2 >= 0 && gate.target_q2 < physical_qubits) {
            std::swap(p2l[gate.target_q1], p2l[gate.target_q2]);
        }
    }

    std::vector<int> l2p(logical_qubits, -1);
    for (int p = 0; p < physical_qubits; ++p) {
        int lq = p2l[p];
        if (lq >= 0 && lq < logical_qubits) l2p[lq] = p;
    }

    for (int lq = 0; lq < logical_qubits; ++lq) {
        if (l2p[lq] == -1) l2p[lq] = lq;
    }

    return l2p;
}

// ZMIANA: Zamiast pisać bezpośrednio do dysku i używać OpenMP critical, piszemy do zwykłego strumienia (bufora)
static void write_layout_row(std::ostream& out, const std::string& circuit_id, const std::string& topology_name,
                             int logical_qubits, int gate_count, const LayoutCandidate& layout) {
    out << CsvWriter::escape(circuit_id) << ','
        << CsvWriter::escape(topology_name) << ','
        << logical_qubits << ','
        << gate_count << ','
        << CsvWriter::escape(layout.strategy) << ','
        << layout.layout_id << ','
        << layout.seed << ','
        << layout.connectivity_score << ','
        << CsvWriter::escape(mapping_to_string(layout.l2p)) << "\n";
}

// ZMIANA: Usunięte omp critical, operujemy na std::ostream
static void write_result_row(std::ostream& out, const std::string& circuit_id, const std::filesystem::path& circuit_path,
                             const std::string& topology_name, int logical_qubits, int gate_count,
                             const LayoutCandidate& layout, const std::string& router, int run_id,
                             uint64_t seed, const HybridRunParams& params, const OptimizationResult& result) {
    out << CsvWriter::escape(circuit_id) << ','
        << CsvWriter::escape(circuit_path.string()) << ','
        << CsvWriter::escape(topology_name) << ','
        << logical_qubits << ','
        << gate_count << ','
        << CsvWriter::escape(layout.strategy) << ','
        << layout.layout_id << ','
        << layout.seed << ','
        << layout.connectivity_score << ','
        << CsvWriter::escape(router) << ','
        << run_id << ','
        << seed << ','
        << CsvWriter::escape(params.profile) << ','
        << params.pop_size << ','
        << params.generations_per_window << ','
        << params.window_size << ','
        << params.accept_ratio << ','
        << params.lazy_offboard << ','
        << params.f << ','
        << params.cr << ','
        << params.lookahead_depth << ','
        << result.score << ','
        << result.makespan << ','
        << result.execution_time << ','
        << result.depth << ','
        << result.swaps << ','
        << result.moves << "\n";
}

static OptimizationResult to_optimization_result(const SabreBenchmarkResult& result) {
    double comparable_score = result.score + (result.makespan * 0.001) + (result.depth * 0.0005);
    return {comparable_score, result.makespan, result.execution_time, result.depth, result.swaps, result.moves};
}

static CircuitMeta make_circuit_meta(const std::filesystem::path& input_file) {
    int logical_qubits = 0;
    std::vector<Gate> circuit = CircuitIO::load_circuitII(input_file.string(), logical_qubits);
    CircuitMeta meta;
    meta.path = input_file;
    meta.circuit_id = circuit_family_from_path(input_file);
    meta.logical_qubits = logical_qubits;
    meta.gate_count = static_cast<int>(circuit.size());
    meta.topology_bucket = topology_bucket_for_qubits(logical_qubits);
    meta.size_bucket = size_bucket_for_gates(meta.gate_count);
    meta.training_bucket = meta.topology_bucket + "__" + meta.size_bucket;
    return meta;
}

static std::map<std::string, std::vector<CircuitMeta>> select_training_circuits(
    const std::vector<CircuitMeta>& base_circuits,
    const ExperimentConfig& config
) {
    std::map<std::string, std::vector<CircuitMeta>> buckets;

    for (const auto& meta : base_circuits) {
        buckets[meta.training_bucket].push_back(meta);
    }

    std::map<std::string, std::vector<CircuitMeta>> selected;
    for (auto& [bucket_name, metas] : buckets) {
        std::mt19937_64 rng(config.base_seed + static_cast<uint64_t>(std::hash<std::string>{}(bucket_name)));
        std::shuffle(metas.begin(), metas.end(), rng);

        int take = std::min(config.hyperparam_training_instances_per_bucket, static_cast<int>(metas.size()));
        selected[bucket_name].insert(selected[bucket_name].end(), metas.begin(), metas.begin() + take);
    }

    return selected;
}

static double tuning_objective(const OptimizationResult& result, const ExperimentConfig& config) {
    return result.score * (1.0 + config.hyperparam_time_penalty * std::log1p(std::max(0.0, result.execution_time)));
}

// ZMIANA: Operujemy na standardowym strumieniu wyjściowym
static void write_hyperparam_row(std::ostream& out, const std::string& bucket, const CircuitMeta& meta,
                                 const HybridRunParams& params, int candidate_id, int run_id,
                                 uint64_t seed, const OptimizationResult& result,
                                 double objective, bool selected) {
    out << CsvWriter::escape(bucket) << ','
        << CsvWriter::escape(meta.circuit_id) << ','
        << CsvWriter::escape(meta.path.string()) << ','
        << meta.logical_qubits << ','
        << meta.gate_count << ','
        << candidate_id << ','
        << CsvWriter::escape(params.profile) << ','
        << params.pop_size << ','
        << params.generations_per_window << ','
        << params.window_size << ','
        << params.accept_ratio << ','
        << params.lazy_offboard << ','
        << params.f << ','
        << params.cr << ','
        << params.lookahead_depth << ','
        << run_id << ','
        << seed << ','
        << result.score << ','
        << result.makespan << ','
        << result.execution_time << ','
        << result.depth << ','
        << result.swaps << ','
        << result.moves << ','
        << objective << ','
        << selected << "\n";
}

static HybridRunParams clamp_params_to_circuit(HybridRunParams params, int gate_count, int logical_qubits, int physical_qubits) {
    params.window_size = std::max(50, std::min(params.window_size, std::max(50, gate_count)));
    params.lazy_offboard = logical_qubits < physical_qubits;
    return params;
}

static std::map<std::string, HybridRunParams> tune_hyperparameters(
    const std::map<std::string, std::vector<CircuitMeta>>& training_circuits,
    const ExperimentConfig& config,
    CsvWriter& tuning_csv
) {
    std::map<std::string, HybridRunParams> best_params_by_bucket;

    for (const auto& [bucket_name, metas] : training_circuits) {
        if (metas.empty()) continue;

        int representative_gate_count = metas.front().gate_count;
        int representative_logical_qubits = metas.front().logical_qubits;
        std::string topology_name;
        auto representative_topology = make_topology_for_qubits(representative_logical_qubits, topology_name);

        std::vector<HybridRunParams> candidates = build_hyperparameter_candidates(
            representative_gate_count,
            representative_logical_qubits,
            representative_topology->get_num_qubits()
        );

        // Zmienne globalne dla całego "wiaderka" (bucket)
        double global_best_mean_objective = std::numeric_limits<double>::max();
        int global_best_candidate_id = -1;
        HybridRunParams global_best_params = fallback_hybrid_params(representative_gate_count, representative_logical_qubits);

        // Wspólny bufor na wszystkie wiersze wynikowe, żeby zapisać je na końcu
        using RowTuple = std::tuple<CircuitMeta, HybridRunParams, int, int, uint64_t, OptimizationResult, double>;
        std::vector<RowTuple> all_buffered_rows;

        // ZRÓWNOLEGLENIE: Tworzymy pulę wątków
        #pragma omp parallel
        {
            // Prywatne zmienne wątku, aby uniknąć konfliktów zapisu (race conditions)
            std::vector<RowTuple> local_buffered_rows;
            double local_best_mean_objective = std::numeric_limits<double>::max();
            int local_best_candidate_id = -1;
            HybridRunParams local_best_params;

            // Dynamiczny podział kandydatów między wątki
            #pragma omp for schedule(dynamic, 1)
            for (int candidate_id = 0; candidate_id < static_cast<int>(candidates.size()); ++candidate_id) {
                double objective_sum = 0.0;
                int objective_count = 0;

                for (const auto& meta : metas) {
                    int logical_qubits = 0;
                    std::vector<Gate> circuit = CircuitIO::load_circuitII(meta.path.string(), logical_qubits);
                    if (circuit.empty()) continue;

                    std::string local_topology_name;
                    auto topology = make_topology_for_qubits(logical_qubits, local_topology_name);
                    HybridRunParams params = clamp_params_to_circuit(
                        candidates[candidate_id],
                        static_cast<int>(circuit.size()),
                        logical_qubits,
                        topology->get_num_qubits()
                    );

                    std::vector<int> compact_layout = make_compact_layout(*topology, circuit, logical_qubits);
                    for (int run = 0; run < config.hyperparam_training_runs; ++run) {
                        uint64_t seed = config.base_seed + static_cast<uint64_t>(candidate_id * 10000 + run);
                        SabreDEHybridRouter hybrid(*topology, circuit, logical_qubits, local_topology_name, false, params.lookahead_depth);
                        hybrid.set_random_seed(seed);
                        hybrid.set_sabre_seed_trials(1);
                        hybrid.set_fixed_initial_mapping(compact_layout);

                        OptimizationResult result = hybrid.optimize_windowed(
                            params.pop_size,
                            params.generations_per_window,
                            params.window_size,
                            params.accept_ratio,
                            params.lazy_offboard,
                            params.f,
                            params.cr
                        );

                        double objective = tuning_objective(result, config);
                        objective_sum += objective;
                        objective_count++;

                        // Dodajemy wiersz do lokalnego bufora wątku
                        local_buffered_rows.push_back({meta, params, candidate_id, run, seed, result, objective});
                    }
                }

                if (objective_count > 0) {
                    double mean_objective = objective_sum / static_cast<double>(objective_count);
                    if (mean_objective < local_best_mean_objective) {
                        local_best_mean_objective = mean_objective;
                        local_best_candidate_id = candidate_id;
                        local_best_params = candidates[candidate_id];
                    }
                }
            } // Koniec pętli for (koniec równoległej dystrybucji)

            // REDUKCJA: Bezpieczne łączenie lokalnych wyników poszczególnych wątków do zmiennych globalnych
            #pragma omp critical(hyperparam_reduce)
            {
                if (local_best_mean_objective < global_best_mean_objective) {
                    global_best_mean_objective = local_best_mean_objective;
                    global_best_candidate_id = local_best_candidate_id;
                    global_best_params = local_best_params;
                }
                // Bezpieczne doklejenie wierszy z wątku do wspólnego wektora
                all_buffered_rows.insert(all_buffered_rows.end(), local_buffered_rows.begin(), local_buffered_rows.end());
            }
        } // Koniec sekcji #pragma omp parallel

        // Oznaczenie profilu dla zwycięskiego zestawu parametrów
        if (global_best_candidate_id != -1) {
            global_best_params.profile = bucket_name + "_best_candidate_" + std::to_string(global_best_candidate_id);
        }

        // Zrzut danych wygenerowanych z testów parametrów – jednorazowo (szybciej i bezpieczniej)
        std::ostringstream local_csv;
        for (const auto& row : all_buffered_rows) {
            const auto& meta = std::get<0>(row);
            const auto& params = std::get<1>(row);
            int candidate_id = std::get<2>(row);
            int run = std::get<3>(row);
            uint64_t seed = std::get<4>(row);
            const auto& result = std::get<5>(row);
            double objective = std::get<6>(row);

            // Flaga "selected" jest ustawiana na true tylko dla kandydata, który wygrał w globalnym rozrachunku
            write_hyperparam_row(local_csv, bucket_name, meta, params, candidate_id, run, seed, result, objective, candidate_id == global_best_candidate_id);
        }
        tuning_csv.write_chunk(local_csv.str());

        best_params_by_bucket[bucket_name] = global_best_params;
    }

    return best_params_by_bucket;
}

static HybridRunParams params_for_circuit(const std::map<std::string, HybridRunParams>& learned_params,
                                          int gate_count, int logical_qubits, int physical_qubits) {
    std::string bucket = training_bucket_for(logical_qubits, gate_count);
    auto it = learned_params.find(bucket);
    HybridRunParams params = (it == learned_params.end()) ? fallback_hybrid_params(gate_count, logical_qubits) : it->second;
    return clamp_params_to_circuit(params, gate_count, logical_qubits, physical_qubits);
}

static std::vector<LayoutCandidate> build_layout_suite(const BaseTopology& topology, const std::vector<Gate>& circuit,
                                                       int logical_qubits, const ExperimentConfig& config,
                                                       int random_layout_count, uint64_t circuit_seed,
                                                       const std::string& topology_name,
                                                       const std::filesystem::path& base_circuit_path) {
    std::vector<LayoutCandidate> layouts;

    auto add_layout = [&](std::string strategy, int layout_id, uint64_t seed, std::vector<int> mapping) {
        if (!validate_mapping(mapping, logical_qubits, topology.get_num_qubits())) {
            throw std::runtime_error("Niepoprawny layout dla strategii: " + strategy);
        }
        double score = layout_connectivity_score(topology, circuit, mapping, logical_qubits);
        layouts.push_back({std::move(strategy), layout_id, seed, std::move(mapping), score});
    };

    std::vector<int> identity = make_identity_layout(logical_qubits, topology.get_num_qubits());
    add_layout("identity_trivial", 0, config.base_seed, identity);

    std::filesystem::path qiskit_path = qiskit_routed_path_for_base(base_circuit_path);
    if (std::filesystem::exists(qiskit_path)) {
        add_layout("qiskit_sabre_initial_layout", 0, circuit_seed + 31,
                   derive_qiskit_layout_from_routed_file(qiskit_path, logical_qubits, topology.get_num_qubits()));
    }

    add_layout("compact_connectivity", 0, circuit_seed + 17, make_compact_layout(topology, circuit, logical_qubits));

    for (int i = 0; i < random_layout_count; ++i) {
        uint64_t seed = circuit_seed + 1000 + static_cast<uint64_t>(i);
        std::mt19937_64 rng(seed);
        add_layout("random", i, seed, make_random_layout(logical_qubits, topology.get_num_qubits(), rng));
    }

    return layouts;
}

static void process_qiskit_routed_file(const std::filesystem::path& input_file, const ExperimentConfig& config,
                                       const std::map<std::string, HybridRunParams>& learned_params,
                                       CsvWriter& results_csv) {
    int logical_qubits = 0;
    std::vector<Gate> qiskit_circuit = CircuitIO::load_circuitII(input_file.string(), logical_qubits);
    if (qiskit_circuit.empty()) return;
    CircuitFileMetadata metadata = CircuitIO::load_metadata(input_file.string());

    std::string topology_name;
    auto topology = make_topology_for_qubits(logical_qubits, topology_name);
    HybridRunParams params = params_for_circuit(learned_params, static_cast<int>(qiskit_circuit.size()), logical_qubits, topology->get_num_qubits());

    std::vector<Gate> empty_logical_circuit;
    SabreDEHybridRouter evaluator(*topology, empty_logical_circuit, logical_qubits, "qiskit_eval", false, params.lookahead_depth);
    evaluator.set_random_seed(config.base_seed);
    OptimizationResult qiskit_result = evaluator.evaluate_qiskit_routed_circuit(qiskit_circuit, metadata.qiskit_initial_layout_l2p);

    std::vector<int> qiskit_initial_layout = metadata.qiskit_initial_layout_l2p;
    if (!validate_mapping(qiskit_initial_layout, logical_qubits, topology->get_num_qubits())) {
        qiskit_initial_layout = make_identity_layout(logical_qubits, topology->get_num_qubits());
    }
    LayoutCandidate qiskit_layout{"qiskit_sabre_initial_layout", 0, 0,
                                  qiskit_initial_layout,
                                  0.0};

    // ZMIANA: Zapis do bufora i jednorazowy zrzut
    std::ostringstream local_results;
    write_result_row(local_results, circuit_family_from_path(input_file), input_file, topology_name,
                     logical_qubits, static_cast<int>(qiskit_circuit.size()), qiskit_layout,
                     "qiskit_sabre_export_eval", 0, config.base_seed, params, qiskit_result);
    results_csv.write_chunk(local_results.str());
}

static void process_base_circuit_file(const std::filesystem::path& input_file, const ExperimentConfig& config,
                                      const std::map<std::string, HybridRunParams>& learned_params,
                                      CsvWriter& results_csv, CsvWriter& layouts_csv, CsvWriter& tuning_csv,
                                      std::map<std::string, int>& random_tuned_instances_by_bucket,
                                      std::map<std::string, int>& random_layout_count_by_bucket) {
    int logical_qubits = 0;
    std::vector<Gate> circuit = CircuitIO::load_circuitII(input_file.string(), logical_qubits);
    if (circuit.empty()) return;

    std::string topology_name;
    auto topology = make_topology_for_qubits(logical_qubits, topology_name);
    std::string circuit_id = circuit_family_from_path(input_file);
    std::string bucket = training_bucket_for(logical_qubits, static_cast<int>(circuit.size()));
    uint64_t circuit_seed = config.base_seed + static_cast<uint64_t>(std::hash<std::string>{}(circuit_id));

    int local_random_layout_count = config.min_random_layouts;
    bool do_tuning = false;

    #pragma omp critical(tuning_maps)
    {
        if (random_layout_count_by_bucket.find(bucket) == random_layout_count_by_bucket.end() ||
            random_layout_count_by_bucket[bucket] == 0) {
            random_layout_count_by_bucket[bucket] = config.min_random_layouts;
        }

        int& tuned_instances_done = random_tuned_instances_by_bucket[bucket];

        if (tuned_instances_done < config.random_tuning_instances_per_bucket) {
            do_tuning = true;
            tuned_instances_done++;
        }

        local_random_layout_count = random_layout_count_by_bucket[bucket];
    }

    // ZMIANA: Tworzymy trzy lokalne bufory wierszy dla aktualnego wątku/obwodu
    std::ostringstream local_tuning_buffer;
    std::ostringstream local_layouts_buffer;
    std::ostringstream local_results_buffer;

    if (do_tuning) {
        std::vector<double> pilot_scores;
        for (int i = 0; i < config.random_pilot_layouts; ++i) {
            std::mt19937_64 rng(circuit_seed + 50000 + static_cast<uint64_t>(i));
            auto layout = make_random_layout(logical_qubits, topology->get_num_qubits(), rng);
            SabreRouter pilot_sabre(*topology, circuit, logical_qubits, topology_name, 20);
            pilot_sabre.set_random_seed(static_cast<uint32_t>(circuit_seed + 70000 + static_cast<uint64_t>(i)));
            pilot_scores.push_back(to_optimization_result(pilot_sabre.run_benchmark_with_initial_mapping(layout, false)).score);
        }

        double raw_estimated = estimate_random_layout_count_raw(pilot_scores, config);
        int estimated = estimate_random_layout_count(pilot_scores, config);

        #pragma omp critical(tuning_maps)
        {
            random_layout_count_by_bucket[bucket] = std::max(random_layout_count_by_bucket[bucket], estimated);
            local_random_layout_count = random_layout_count_by_bucket[bucket];
        }

        double mean = std::accumulate(pilot_scores.begin(), pilot_scores.end(), 0.0) / static_cast<double>(pilot_scores.size());
        double sq_sum = 0.0;
        for (double value : pilot_scores) sq_sum += (value - mean) * (value - mean);
        double sd = pilot_scores.size() > 1 ? std::sqrt(sq_sum / static_cast<double>(pilot_scores.size() - 1)) : 0.0;

        // ZMIANA: Brak critical podczas formatowania wiersza – piszemy do bufora w pamięci
        local_tuning_buffer << CsvWriter::escape(circuit_id) << ','
                            << CsvWriter::escape(bucket) << ','
                            << logical_qubits << ','
                            << circuit.size() << ','
                            << config.random_pilot_layouts << ','
                            << mean << ','
                            << sd << ','
                            << raw_estimated << ','
                            << local_random_layout_count << ','
                            << "sabre_fixed_layout_score" << "\n";
    }

    HybridRunParams params = params_for_circuit(learned_params, static_cast<int>(circuit.size()), logical_qubits, topology->get_num_qubits());
    std::vector<LayoutCandidate> layouts = build_layout_suite(*topology, circuit, logical_qubits, config,
                                                              local_random_layout_count, circuit_seed, topology_name, input_file);

    for (const auto& layout : layouts) {
        // ZMIANA: Zapis bezpośrednio do zadeklarowanego wcześniej bufora lokalnego
        write_layout_row(local_layouts_buffer, circuit_id, topology_name, logical_qubits, static_cast<int>(circuit.size()), layout);

        for (int run = 0; run < config.deterministic_runs; ++run) {
            uint64_t run_seed = circuit_seed + static_cast<uint64_t>(run);

            // 1. Uruchomienie Sabre
            HybridRunParams sabre_params = params;
            sabre_params.lookahead_depth = 20;
            sabre_params.profile = params.profile + "_sabre_lh20";

            SabreRouter sabre(*topology, circuit, logical_qubits, topology_name, 20);
            sabre.set_random_seed(static_cast<uint32_t>(run_seed));
            SabreBenchmarkResult sabre_result = sabre.run_benchmark_with_initial_mapping(layout.l2p, false);

            write_result_row(local_results_buffer, circuit_id, input_file, topology_name, logical_qubits,
                             static_cast<int>(circuit.size()), layout, "sabre_fixed_layout",
                             run, run_seed, sabre_params, to_optimization_result(sabre_result));

            // 2. Uruchomienie DEHybrid
            SabreDEHybridRouter hybrid(*topology, circuit, logical_qubits, topology_name, false, params.lookahead_depth);
            hybrid.set_random_seed(run_seed);
            hybrid.set_sabre_seed_trials(1);
            hybrid.set_fixed_initial_mapping(layout.l2p);
            OptimizationResult hybrid_result = hybrid.optimize_windowed(
                params.pop_size,
                params.generations_per_window,
                params.window_size,
                params.accept_ratio,
                params.lazy_offboard,
                params.f,
                params.cr
            );

            write_result_row(local_results_buffer, circuit_id, input_file, topology_name, logical_qubits,
                             static_cast<int>(circuit.size()), layout, "dehybrid_sabre_fixed_layout",
                             run, run_seed, params, hybrid_result);
        }
    }

    // ZMIANA NAJWAŻNIEJSZA: Po zakończeniu wszystkich obliczeń zrzucamy paczkę danych do głównego pliku naraz.
    tuning_csv.write_chunk(local_tuning_buffer.str());
    layouts_csv.write_chunk(local_layouts_buffer.str());
    results_csv.write_chunk(local_results_buffer.str());
}

// int main(int argc, char* argv[]) {
// #ifdef _WIN32
//     SetConsoleOutputCP(CP_UTF8);
//     SetConsoleCP(CP_UTF8);
// #endif
//
//     ExperimentConfig config;
//     if (argc >= 2) config.input_dir = argv[1];
//     if (argc >= 3) config.results_dir = argv[2];
//
//     if (!std::filesystem::exists(config.input_dir) || !std::filesystem::is_directory(config.input_dir)) {
//         std::cerr << "Blad: katalog benchmarkow nie istnieje: " << config.input_dir << "\n";
//         return -1;
//     }
//
//     CsvWriter results_csv(
//         config.results_dir / "routing_runs.csv",
//         "circuit_id,circuit_path,topology,logical_qubits,gates,layout_strategy,layout_id,layout_seed,layout_connectivity_score,router,run_id,run_seed,param_profile,pop_size,generations_per_window,window_size,accept_ratio,lazy_offboard,f,cr,lookahead_depth,score,makespan,execution_time_s,depth,swaps,moves"
//     );
//     CsvWriter layouts_csv(
//         config.results_dir / "layouts.csv",
//         "circuit_id,topology,logical_qubits,gates,layout_strategy,layout_id,layout_seed,layout_connectivity_score,mapping_l2p"
//     );
//     CsvWriter tuning_csv(
//         config.results_dir / "random_layout_tuning.csv",
//         "circuit_id,training_bucket,logical_qubits,gates,pilot_layouts,mean_score,sd_score,raw_estimated_random_layouts,used_random_layouts,metric"
//     );
//     CsvWriter hyperparam_csv(
//         config.results_dir / "hyperparameter_tuning.csv",
//         "training_bucket,circuit_id,circuit_path,logical_qubits,gates,candidate_id,param_profile,pop_size,generations_per_window,window_size,accept_ratio,lazy_offboard,f,cr,lookahead_depth,run_id,seed,score,makespan,execution_time_s,depth,swaps,moves,objective,selected"
//     );
//
//     std::vector<std::filesystem::path> files;
//     std::vector<CircuitMeta> base_metas;
//     for (const auto& entry : std::filesystem::recursive_directory_iterator(config.input_dir)) {
//         if (entry.is_regular_file() && entry.path().extension() == ".json") {
//
//             int temp_qubits = 0;
//             std::vector<Gate> check_circuit = CircuitIO::load_circuitII(entry.path().string(), temp_qubits);
//             // if (check_circuit.size() > 2000) {
//             //     continue;
//             // }
//
//             files.push_back(entry.path());
//             if (!is_routed_qiskit_file(entry.path())) {
//                 CircuitMeta meta = make_circuit_meta(entry.path());
//                 if (meta.gate_count > 0 && meta.logical_qubits > 0) {
//                     base_metas.push_back(meta);
//                 }
//             }
//         }
//     }
//
//     std::sort(files.begin(), files.end());
//     std::sort(base_metas.begin(), base_metas.end(), [](const CircuitMeta& a, const CircuitMeta& b) {
//         if (a.training_bucket != b.training_bucket) return a.training_bucket < b.training_bucket;
//         if (a.gate_count != b.gate_count) return a.gate_count < b.gate_count;
//         return a.circuit_id < b.circuit_id;
//     });
//
//     auto training_circuits = select_training_circuits(base_metas, config);
//     auto learned_params = tune_hyperparameters(training_circuits, config, hyperparam_csv);
//
//     // --- NOWE: Wypisywanie wyników tuningu do konsoli ---
//     std::cout << "\n--- ZAKONCZONO TUNING HYPERPARAMETROW ---\n";
//     for (const auto& [bucket, params] : learned_params) {
//         std::cout << "Baza [" << bucket << "] -> Zwyciezca: " << params.profile
//                   << " (Pop: " << params.pop_size
//                   << ", Gen: " << params.generations_per_window
//                   << ", Window: " << params.window_size << ")\n";
//     }
//     std::cout << "-----------------------------------------\n\n";
//
//     std::map<std::string, int> random_tuned_instances_by_bucket;
//     std::map<std::string, int> random_layout_count_by_bucket;
//
//     int max_threads = omp_get_max_threads();
//     std::cout << "\nRozpoczynam zrownoleglone przetwarzanie na " << max_threads << " watkach.\n";
//
//     std::atomic<int> processed_files{0};
//     int total_files = static_cast<int>(files.size());
//
//     #pragma omp parallel for schedule(dynamic, 1)
//     for (int i = 0; i < total_files; ++i) {
//         const auto& input_file = files[i];
//
//         if (is_routed_qiskit_file(input_file)) {
//             process_qiskit_routed_file(input_file, config, learned_params, results_csv);
//         } else {
//             process_base_circuit_file(input_file, config, learned_params, results_csv, layouts_csv, tuning_csv,
//                                       random_tuned_instances_by_bucket, random_layout_count_by_bucket);
//         }
//
//         // --- NOWE: Lekkie logowanie postępu (nadpisuje jedną linię w konsoli) ---
//         int done = ++processed_files;
//         if (done % 10 == 0 || done == total_files) {
//             // Budujemy cały napis w pamięci, dodając spacje na końcu do wyczyszczenia bufora
//             std::string msg = "\rPostep: " + std::to_string(done) + " / " + std::to_string(total_files) +
//                               " (" + std::to_string(done * 100 / total_files) + "%) plikow przetworzonych.   ";
//             #pragma omp critical(console_log)
//             {
//                 // Wysyłamy gotowy ciąg znaków w jednej atomowej operacji
//                 std::cout << msg << std::flush;
//             }
//         }
//     }
//
//     std::cout << "\nPrzetwarzanie zakonczone z sukcesem!\n";
//     return 0;
// }

// int main(int argc, char* argv[]) {
// #ifdef _WIN32
//     SetConsoleOutputCP(CP_UTF8);
//     SetConsoleCP(CP_UTF8);
// #endif
//
//     ExperimentConfig config;
//     if (argc >= 2) config.input_dir = argv[1];
//     if (argc >= 3) config.results_dir = argv[2];
//
//     if (!std::filesystem::exists(config.input_dir) || !std::filesystem::is_directory(config.input_dir)) {
//         std::cerr << "Blad: katalog benchmarkow nie istnieje: " << config.input_dir << "\n";
//         return -1;
//     }
//
//     CsvWriter results_csv(
//         config.results_dir / "routing_run.csv",
//         "circuit_id,circuit_path,topology,logical_qubits,gates,layout_strategy,layout_id,layout_seed,layout_connectivity_score,router,run_id,run_seed,param_profile,pop_size,generations_per_window,window_size,accept_ratio,lazy_offboard,f,cr,lookahead_depth,score,makespan,execution_time_s,depth,swaps,moves"
//     );
//
//     std::vector<std::filesystem::path> files;
//     for (const auto& entry : std::filesystem::recursive_directory_iterator(config.input_dir)) {
//         if (entry.is_regular_file() && entry.path().extension() == ".json") {
//
//             // Tylko pliki Qiskit
//             if (!is_routed_qiskit_file(entry.path())) continue;
//
//             int temp_qubits = 0;
//             std::vector<Gate> check_circuit = CircuitIO::load_circuitII(entry.path().string(), temp_qubits);
//
//             // POPRAWIONY WARUNEK: Od 13 kubitów w górę i do 2000 bramek włącznie
//             if (check_circuit.size() < 3500) {
//                 files.push_back(entry.path());
//             }
//         }
//     }
//
//     std::sort(files.begin(), files.end());
//     std::map<std::string, HybridRunParams> learned_params;
//
//     int max_threads = omp_get_max_threads();
//     std::cout << "\nZnaleziono " << files.size() << " plikow qiskit_routed spelniajacych kryteria (>=13q, <=2000g).\n";
//     std::cout << "Rozpoczynam zrownoleglone przetwarzanie na " << max_threads << " watkach.\n";
//
//     std::atomic<int> processed_files{0};
//     int total_files = static_cast<int>(files.size());
//
//     #pragma omp parallel for schedule(dynamic, 1)
//     for (int i = 0; i < total_files; ++i) {
//         const auto& input_file = files[i];
//
//         process_qiskit_routed_file(input_file, config, learned_params, results_csv);
//
//         int done = ++processed_files;
//         if (done % 10 == 0 || done == total_files) {
//             std::string msg = "\rPostep: " + std::to_string(done) + " / " + std::to_string(total_files) +
//                               " (" + std::to_string(done * 100 / total_files) + "%) plikow przetworzonych.   ";
//             #pragma omp critical(console_log)
//             {
//                 std::cout << msg << std::flush;
//             }
//         }
//     }
//
//     std::cout << "\nPrzetwarzanie zakonczone z sukcesem! Wyniki zostaly dopisane do routing_run.csv\n";
//     return 0;
// }