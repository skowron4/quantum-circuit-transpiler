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
    SabreRouter sabre(qpu, circuit, qpu.get_num_qubits(), name, sabre_lookahead_depth);
    sabre.set_random_seed(use_random_seeds ? rd() : sabre_seed);
    sabre.run_benchmark_multiple(sabre_runs);

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

//void run_test_for_topologyII(const BaseTopology &qpu, const std::string &name) {
//    std::cout << "\n============================================\n";
//    std::cout << "URUCHAMIAM TEST DLA TOPOLOGII: " << name << "\n";
//    std::cout << "============================================\n";
//    qpu.print_stats();
//
//    int logical_qubits_count = qpu.get_num_qubits() * 1.0;
//    CircuitGenerator generator(qpu, logical_qubits_count, 1234);
//
//    std::vector<Gate> circuit = generator.generate(2000, 0.6, 0.4);
//    generator.print_gate_statistics(circuit);
//    generator.print_circuit(circuit);
//
//    // --- Konfiguracja losowosci i liczby prob ---
//    const bool use_random_seeds = true;
//    std::random_device rd;
//    const uint32_t sabre_seed = rd();
//    const uint64_t hybrid_seed = sabre_seed;
//    const int sabre_runs = 1;
//    const int hybrid_runs = 5;
//    const int hybrid_sabre_seed_trials = 1;
//    const int sabre_lookahead_depth = 30;
//    const int hybrid_lookahead_depth = 10;
//    const int hybrid_pop_size = 30;
//    const int hybrid_generations_per_window = 300;
//    const int hybrid_window_size = 600;
//    const double hybrid_accept_ratio = 0.35;
//    const int hybrid_feedback_passes = 3;
//    const bool hybrid_allow_lazy_offboard = false;
//    const double hybrid_f = 0.8;
//    const double hybrid_cr = 0.9;
//
//    SabreRouter sabre(qpu, circuit, qpu.get_num_qubits(), name, sabre_lookahead_depth);
//    sabre.set_random_seed(use_random_seeds ? rd() : sabre_seed);
//    sabre.run_benchmark_multiple(sabre_runs);
//
//    std::cout << "\n--- Rozpoczynam " << hybrid_runs << " iteracji dla DEHybridSabre ---\n";
//
//    std::vector<OptimizationResult> results;
//    results.reserve(hybrid_runs);
//
//    for (int i = 0; i < hybrid_runs; ++i) {
//        SabreDEHybridRouter SabreDERouter(qpu, circuit, qpu.get_num_qubits(), name, false, hybrid_lookahead_depth);
//        SabreDERouter.set_random_seed(use_random_seeds ? rd() : hybrid_seed + i);
//        SabreDERouter.set_sabre_seed_trials(hybrid_sabre_seed_trials);
//        OptimizationResult res = SabreDERouter.optimize_windowed(30, 100, 500, 0.6, false, 0.2);
//
//        // Odbieramy całą strukturę wyników
////        OptimizationResult res = SabreDERouter.optimize_windowed_with_feedback(
////            hybrid_pop_size,
////            hybrid_generations_per_window,
////            hybrid_window_size,
////            hybrid_accept_ratio,
////            hybrid_feedback_passes,
////            hybrid_allow_lazy_offboard,
////            hybrid_f,
////            hybrid_cr
////        );
//        results.push_back(res);
//
//        std::cout << "Iteracja " << (i + 1) << " zakonczona. "
//                  << "Wynik: " << res.score << ", "
//                  << "SWAPy: " << res.swaps << ", "
//                  << "Makespan: " << res.makespan << "\n";
//    }
//
//    // --- OBLICZANIE STATYSTYK ---
//    OptimizationResult min_val = results[0];
//    OptimizationResult max_val = results[0];
//    OptimizationResult sum_val = results[0];
//    OptimizationResult best_result = results[0]; // Zmienna do śledzenia pojedynczego najlepszego wyniku
//
//    for (int i = 1; i < hybrid_runs; ++i) {
//        // Zapisujemy całą instancję wyniku, jeśli ma lepszy (niższy) score
//        if (results[i].score < best_result.score) {
//            best_result = results[i];
//        }
//
//        min_val.score = std::min(min_val.score, results[i].score);
//        max_val.score = std::max(max_val.score, results[i].score);
//        sum_val.score += results[i].score;
//
//        min_val.makespan = std::min(min_val.makespan, results[i].makespan);
//        max_val.makespan = std::max(max_val.makespan, results[i].makespan);
//        sum_val.makespan += results[i].makespan;
//
//        min_val.depth = std::min(min_val.depth, results[i].depth);
//        max_val.depth = std::max(max_val.depth, results[i].depth);
//        sum_val.depth += results[i].depth;
//
//        min_val.swaps = std::min(min_val.swaps, results[i].swaps);
//        max_val.swaps = std::max(max_val.swaps, results[i].swaps);
//        sum_val.swaps += results[i].swaps;
//
//        min_val.moves = std::min(min_val.moves, results[i].moves);
//        max_val.moves = std::max(max_val.moves, results[i].moves);
//        sum_val.moves += results[i].moves;
//    }
//
//    std::cout << "\n============================================\n";
//    std::cout << " PODSUMOWANIE DEHybridSabre (" << hybrid_runs << " uruchomien)\n";
//    std::cout << "============================================\n";
//
//    // --- WYPISYWANIE NAJLEPSZEGO POJEDYNCZEGO WYNIKU ---
//    std::cout << " [NAJLEPSZY ZNALEZIONY WYNIK (wg min Score)]\n";
//    std::cout << " Score: " << best_result.score << "\n";
//    std::cout << " Makespan: " << best_result.makespan << "\n";
//    std::cout << " Depth: " << best_result.depth << "\n";
//    std::cout << " Swaps: " << best_result.swaps << "\n";
//    std::cout << " Moves: " << best_result.moves << "\n";
//    std::cout << "--------------------------------------------\n";
//
//    std::cout << " [SCORE]\n";
//    std::cout << " Min: " << min_val.score << " | Max: " << max_val.score << " | Średnia: " << (sum_val.score / hybrid_runs) << "\n";
//
//    std::cout << " [MAKESPAN]\n";
//    std::cout << " Min: " << min_val.makespan << " | Max: " << max_val.makespan << " | Średnia: " << (sum_val.makespan / hybrid_runs) << "\n";
//
//    std::cout << " [DEPTH]\n";
//    std::cout << " Min: " << min_val.depth << " | Max: " << max_val.depth << " | Średnia: " << (static_cast<double>(sum_val.depth) / hybrid_runs) << "\n";
//
//    std::cout << " [SWAPS]\n";
//    std::cout << " Min: " << min_val.swaps << " | Max: " << max_val.swaps << " | Średnia: " << (static_cast<double>(sum_val.swaps) / hybrid_runs) << "\n";
//
//    std::cout << " [MOVES]\n";
//    std::cout << " Min: " << min_val.moves << " | Max: " << max_val.moves << " | Średnia: " << (static_cast<double>(sum_val.moves) / hybrid_runs) << "\n";
//    std::cout << "============================================\n";
//}
//
// int main() {
// #ifdef _WIN32
//     // Windows potrzebuje instrukcji zmiany kodowania na UTF-8
//     SetConsoleOutputCP(CP_UTF8);
//     // Opcjonalnie: upewnia się, że wejście (cin) też używa UTF-8
//     SetConsoleCP(CP_UTF8);
// #endif
//     // Topologia siatki 2x2
//       IQMConstellation gridQpu(5, 5);
////        run_test_for_topology(gridQpu, "Siatka 5x5");
//       run_test_for_topologyII(gridQpu, "Siatka 5x5");
//
//     // Topologia trójkątna
////     IQMTriangleConstellation triangleQpu;
////     run_test_for_topology(triangleQpu, "Trójkąt 3-Rezonatory");
//     // run_test_for_topologyII(triangleQpu, "Trójkąt 3-Rezonatory");
//
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
    int deterministic_runs = 3;
    int random_pilot_layouts = 8;
    int random_tuning_instances = 2;
    int min_random_layouts = 5;
    int max_random_layouts = 40;
    double random_target_relative_error = 0.10;
    uint64_t base_seed = 1337;
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

static HybridRunParams choose_hybrid_params(int gate_count, int logical_qubits) {
    HybridRunParams params;

    if (gate_count <= 250 || logical_qubits <= 8) {
        params.profile = "small_fast";
        params.pop_size = 20;
        params.generations_per_window = 100;
        params.window_size = std::max(80, std::min(250, gate_count));
        params.accept_ratio = 0.60;
    } else if (gate_count <= 1500) {
        params.profile = "medium_balanced";
        params.pop_size = 30;
        params.generations_per_window = 180;
        params.window_size = std::max(250, std::min(700, gate_count / 2));
        params.accept_ratio = 0.70;
    } else {
        params.profile = "large_windowed";
        params.pop_size = 36;
        params.generations_per_window = 140;
        params.window_size = std::max(700, std::min(1800, gate_count / 3));
        params.accept_ratio = 0.75;
    }

    return params;
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

static void write_layout_row(CsvWriter& writer, const std::string& circuit_id, const std::string& topology_name,
                             int logical_qubits, int gate_count, const LayoutCandidate& layout) {
    writer.out << CsvWriter::escape(circuit_id) << ','
               << CsvWriter::escape(topology_name) << ','
               << logical_qubits << ','
               << gate_count << ','
               << CsvWriter::escape(layout.strategy) << ','
               << layout.layout_id << ','
               << layout.seed << ','
               << layout.connectivity_score << ','
               << CsvWriter::escape(mapping_to_string(layout.l2p)) << "\n";
}

static void write_result_row(CsvWriter& writer, const std::string& circuit_id, const std::filesystem::path& circuit_path,
                             const std::string& topology_name, int logical_qubits, int gate_count,
                             const LayoutCandidate& layout, const std::string& router, int run_id,
                             uint64_t seed, const HybridRunParams& params, const OptimizationResult& result) {
    writer.out << CsvWriter::escape(circuit_id) << ','
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
               << result.depth << ','
               << result.swaps << ','
               << result.moves << "\n";
}

static OptimizationResult to_optimization_result(const SabreBenchmarkResult& result) {
    double comparable_score = result.score + (result.makespan * 0.001) + (result.depth * 0.0005);
    return {comparable_score, result.makespan, result.depth, result.swaps, result.moves};
}

static std::vector<LayoutCandidate> build_layout_suite(const BaseTopology& topology, const std::vector<Gate>& circuit,
                                                       int logical_qubits, const ExperimentConfig& config,
                                                       int random_layout_count, uint64_t circuit_seed,
                                                       const std::string& topology_name) {
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

    SabreRouter sabre_layout_proxy(topology, circuit, logical_qubits, topology_name, 20);
    sabre_layout_proxy.set_random_seed(static_cast<uint32_t>(circuit_seed));
    add_layout("sabre_layout_proxy_qiskit", 0, circuit_seed, sabre_layout_proxy.compute_bidirectional_initial_mapping(identity, false));

    add_layout("compact_connectivity", 0, circuit_seed + 17, make_compact_layout(topology, circuit, logical_qubits));

    for (int i = 0; i < random_layout_count; ++i) {
        uint64_t seed = circuit_seed + 1000 + static_cast<uint64_t>(i);
        std::mt19937_64 rng(seed);
        add_layout("random", i, seed, make_random_layout(logical_qubits, topology.get_num_qubits(), rng));
    }

    return layouts;
}

static void process_qiskit_routed_file(const std::filesystem::path& input_file, const ExperimentConfig& config,
                                       CsvWriter& results_csv) {
    int logical_qubits = 0;
    std::vector<Gate> qiskit_circuit = CircuitIO::load_circuitII(input_file.string(), logical_qubits);
    if (qiskit_circuit.empty()) return;

    std::string topology_name;
    auto topology = make_topology_for_qubits(logical_qubits, topology_name);
    HybridRunParams params = choose_hybrid_params(static_cast<int>(qiskit_circuit.size()), logical_qubits);

    SabreDEHybridRouter evaluator(*topology, qiskit_circuit, logical_qubits, "qiskit_eval", false, params.lookahead_depth);
    evaluator.set_random_seed(config.base_seed);
    OptimizationResult qiskit_result = evaluator.evaluate_qiskit_routed_circuit(qiskit_circuit);

    LayoutCandidate qiskit_layout{"qiskit_final_layout_in_file", 0, 0,
                                  make_identity_layout(logical_qubits, topology->get_num_qubits()), 0.0};
    write_result_row(results_csv, circuit_family_from_path(input_file), input_file, topology_name,
                     logical_qubits, static_cast<int>(qiskit_circuit.size()), qiskit_layout,
                     "qiskit_sabre_export_eval", 0, config.base_seed, params, qiskit_result);
}

static void process_base_circuit_file(const std::filesystem::path& input_file, const ExperimentConfig& config,
                                      CsvWriter& results_csv, CsvWriter& layouts_csv, CsvWriter& tuning_csv,
                                      int& tuned_instances_done, int& random_layout_count) {
    int logical_qubits = 0;
    std::vector<Gate> circuit = CircuitIO::load_circuitII(input_file.string(), logical_qubits);
    if (circuit.empty()) return;

    std::string topology_name;
    auto topology = make_topology_for_qubits(logical_qubits, topology_name);
    std::string circuit_id = circuit_family_from_path(input_file);
    uint64_t circuit_seed = config.base_seed + static_cast<uint64_t>(std::hash<std::string>{}(circuit_id));

    if (tuned_instances_done < config.random_tuning_instances) {
        std::vector<double> pilot_scores;
        for (int i = 0; i < config.random_pilot_layouts; ++i) {
            std::mt19937_64 rng(circuit_seed + 50000 + static_cast<uint64_t>(i));
            auto layout = make_random_layout(logical_qubits, topology->get_num_qubits(), rng);
            SabreRouter pilot_sabre(*topology, circuit, logical_qubits, topology_name, 4);
            pilot_sabre.set_random_seed(static_cast<uint32_t>(circuit_seed + 70000 + static_cast<uint64_t>(i)));
            pilot_scores.push_back(to_optimization_result(pilot_sabre.run_benchmark_with_initial_mapping(layout, false)).score);
        }

        double raw_estimated = estimate_random_layout_count_raw(pilot_scores, config);
        int estimated = estimate_random_layout_count(pilot_scores, config);
        random_layout_count = std::max(random_layout_count, estimated);

        double mean = std::accumulate(pilot_scores.begin(), pilot_scores.end(), 0.0) / static_cast<double>(pilot_scores.size());
        double sq_sum = 0.0;
        for (double value : pilot_scores) sq_sum += (value - mean) * (value - mean);
        double sd = pilot_scores.size() > 1 ? std::sqrt(sq_sum / static_cast<double>(pilot_scores.size() - 1)) : 0.0;

        tuning_csv.out << CsvWriter::escape(circuit_id) << ','
                       << logical_qubits << ','
                       << circuit.size() << ','
                       << config.random_pilot_layouts << ','
                       << mean << ','
                       << sd << ','
                       << raw_estimated << ','
                       << random_layout_count << ','
                       << "sabre_fixed_layout_score" << "\n";
        tuned_instances_done++;
    }

    HybridRunParams params = choose_hybrid_params(static_cast<int>(circuit.size()), logical_qubits);
    std::vector<LayoutCandidate> layouts = build_layout_suite(*topology, circuit, logical_qubits, config,
                                                              random_layout_count, circuit_seed, topology_name);

    for (const auto& layout : layouts) {
        write_layout_row(layouts_csv, circuit_id, topology_name, logical_qubits, static_cast<int>(circuit.size()), layout);

        for (int run = 0; run < config.deterministic_runs; ++run) {
            uint64_t run_seed = circuit_seed + static_cast<uint64_t>(run);

            SabreRouter sabre(*topology, circuit, logical_qubits, topology_name, params.lookahead_depth);
            sabre.set_random_seed(static_cast<uint32_t>(run_seed));
            SabreBenchmarkResult sabre_result = sabre.run_benchmark_with_initial_mapping(layout.l2p, false);
            write_result_row(results_csv, circuit_id, input_file, topology_name, logical_qubits,
                             static_cast<int>(circuit.size()), layout, "sabre_fixed_layout",
                             run, run_seed, params, to_optimization_result(sabre_result));

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

            write_result_row(results_csv, circuit_id, input_file, topology_name, logical_qubits,
                             static_cast<int>(circuit.size()), layout, "dehybrid_sabre_fixed_layout",
                             run, run_seed, params, hybrid_result);
        }
    }
}

int main(int argc, char* argv[]) {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif

    ExperimentConfig config;
    if (argc >= 2) config.input_dir = argv[1];
    if (argc >= 3) config.results_dir = argv[2];

    if (!std::filesystem::exists(config.input_dir) || !std::filesystem::is_directory(config.input_dir)) {
        std::cerr << "Blad: katalog benchmarkow nie istnieje: " << config.input_dir << "\n";
        return -1;
    }

    CsvWriter results_csv(
        config.results_dir / "routing_runs.csv",
        "circuit_id,circuit_path,topology,logical_qubits,gates,layout_strategy,layout_id,layout_seed,layout_connectivity_score,router,run_id,run_seed,param_profile,pop_size,generations_per_window,window_size,accept_ratio,lazy_offboard,f,cr,lookahead_depth,score,makespan,depth,swaps,moves"
    );
    CsvWriter layouts_csv(
        config.results_dir / "layouts.csv",
        "circuit_id,topology,logical_qubits,gates,layout_strategy,layout_id,layout_seed,layout_connectivity_score,mapping_l2p"
    );
    CsvWriter tuning_csv(
        config.results_dir / "random_layout_tuning.csv",
        "circuit_id,logical_qubits,gates,pilot_layouts,mean_score,sd_score,raw_estimated_random_layouts,used_random_layouts,metric"
    );

    int tuned_instances_done = 0;
    int random_layout_count = config.min_random_layouts;

    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(config.input_dir)) {
        if (entry.is_regular_file() && entry.path().extension() == ".json") {
            files.push_back(entry.path());
        }
    }
    std::sort(files.begin(), files.end());

    for (const auto& input_file : files) {
        std::cout << "\n=> Przetwarzam: " << input_file << "\n";
        if (is_routed_qiskit_file(input_file)) {
            process_qiskit_routed_file(input_file, config, results_csv);
        } else {
            process_base_circuit_file(input_file, config, results_csv, layouts_csv, tuning_csv,
                                      tuned_instances_done, random_layout_count);
        }
    }

    std::cout << "\nZakonczono. CSV zapisane w: " << config.results_dir << "\n";
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