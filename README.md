# Hybrid SABRE-DE Quantum Circuit Router for Resonator-Based Topologies

![C++20](https://img.shields.io/badge/C%2B%2B20-Combinatorial_Optimization-00599C) ![Domain](https://img.shields.io/badge/Domain-Quantum_Compilation_%26_Routing-purple) ![Python](https://img.shields.io/badge/Python-Qiskit_Benchmarks-3776AB) ![Platform](https://img.shields.io/badge/Platform-x86__64_AVX2-informational) ![Type](https://img.shields.io/badge/Type-Master's_Thesis_R%26D-orange)

## Overview

This repository contains a C++20 combinatorial optimization engine for routing and scheduling logical quantum circuits onto resonator-based ("star/crystal") physical topologies modeled on IQM Constellation hardware. Unlike standard nearest-neighbor grid architectures, physical qubits in this model do not couple directly; executing a two-qubit gate requires moving a state onto a shared microwave resonator (`MOVE`), executing the interaction, and lazily off-boarding the state back to a physical node when required by subsequent operations. The core solver (`SabreDEHybridRouter`) combines a bidirectional SABRE heuristic with a sliding-window Differential Evolution (DE/rand/1/bin) metaheuristic to simultaneously optimize gate execution order, resonator assignment, routing paths, and lazy qubit relocation.

## System Architecture & Pipeline

The compilation and optimization pipeline executes in six stages:

1. **Topology & All-Pairs Path Precomputation (`BaseTopology`, `IQMConstellation`):**
   Constructs a hexagonal resonator grid (up to 5x5 resonators / 70 physical qubits) where resonators are adjacent if they share physical qubits. Computes all-pairs shortest paths over the resonator graph in O(R^3) via Floyd-Warshall over flat row-major matrices (`distance_matrix`, `next_node_matrix`), and precomputes O(Q^2) routing paths between all physical qubit pairs for O(1) reference queries (`get_routing_path`) during routing.

2. **DAG Construction & Gate Folding (`HybridSabreDAG`):**
   Parses the input circuit and folds consecutive single-qubit (1Q) operations on the same wire into the subsequent two-qubit (2Q) gate node (`q1_1q_count`, `q2_1q_count`, `absl::InlinedVector<int, 4>`), reducing graph vertices and scheduling overhead. Builds both a forward DAG and a reversed DAG (`reversed_dag`) for bidirectional passes.

3. **Multi-Trial SABRE Seeding (`prepare_sabre_seed`):**
   Generates initial logical-to-physical mappings using a cascade of heuristics: identity mapping, degree-based mapping, and a time-decay interaction graph (`ALPHA = 0.8`, weight 0.8^layer). Executes forward-backward-forward SABRE passes to extract a baseline execution trace and ranks seed profiles using a composite metric (`0.5 * DE_score + 0.5 * SABRE_cost`), storing the top 4 profiles in `sabre_seed_profile_bank`.

4. **Genotype Encoding & Population Initialization (`initialize_de_population`):**
   Encodes the scheduling and routing problem into a continuous genotype `std::vector<double>` in [0, 1] divided into 4 or 5 contiguous blocks of size `N_gates`:
   - `PRIORITY_BLOCK`: Topological ready-set execution priority.
   - `RESONATOR_BLOCK`: Target resonator selection when multiple shared resonators exist.
   - `QUBIT_IN_BLOCK`: Selection of which qubit moves onto the resonator.
   - `ROUTE_BIAS_BLOCK`: Physical edge weight modifier in [0.95, 1.10] (`0.95 + 0.15 * gene`) applied during heuristic Swap/Move evaluation.
   - `LAZY_TARGET_BLOCK`: Landing physical node selection during lazy off-boarding (active when logical qubits < physical qubits).
   The population is initialized with exact SABRE seeds, uniformly perturbed seed variants (+/- 0.05 noise), and uniform random individuals.

5. **Sliding-Window Differential Evolution with State Snapshots (`optimize_windowed`):**
   Rather than evolving all `5 * N_gates` genes simultaneously—which degrades convergence on deep circuits—the optimizer slides a window (`window_size`, `accept_ratio`) across SABRE's chronological execution order (`gate_time_map`). Upon completing a window, `propagate_window_success` locks the accepted prefix genes across the entire population and captures a `DEStateSnapshot` (complete hardware mappings, resonator occupancy, ready times, DAG dependencies, and decay counters). Subsequent window evaluations resume directly from `load_snapshot`, eliminating redundant prefix simulation. Optional iterative feedback passes (`optimize_windowed_with_feedback`) run the best final mapping backward through the reversed DAG to refine the initial layout with inter-pass elitism.

6. **Hardware Timing Simulation & Verification (`validate_logical_circuit`):**
   Simulates hardware execution with `MOVE_COST = 1.0`, `SWAP_COST = 3.14`, relocation via empty physical qubit at `2.0 * MOVE_COST`, and 1Q/2Q gate durations of `1.0`. Evaluates individuals via `routing_cost + 0.001 * makespan + 0.0005 * logical_depth`. Final schedules are verified by a deterministic replay validator checking DAG causality, resonator topology constraints, resource collisions, and time-overlap invariants.

## Engineering & Performance Highlights

- **Zero-Allocation Hot Loop & Workspace Reuse:**
  Fitness evaluation (`evaluate_fitness_impl`) executes thousands of times per window. To avoid heap contention and allocator overhead, `HybridWorkspace` pre-allocates all simulation buffers once during router construction based on physical qubit count (P), resonator count (R), and gate count (G). `DERunContext` binds references directly to the workspace.
- **Data-Oriented Design (SoA & CSR Layouts):**
  - Physical swap adjacency is stored in a flat Compressed Sparse Row (CSR) layout (`flat_precomputed_swaps` and `precomputed_swaps_offsets`) of 8-byte `SwapEdge` structs.
  - Qubit interaction lists use a Structure-of-Arrays (SoA) layout (`inter_targets` and `inter_weights`) with a fixed stride of `MAX_INTER_PER_Q = 64` per qubit.
  - Lookahead layers and population matrices (`pop_size * n_genes`) are flattened into contiguous 1D `std::vector` buffers to maximize L1/L2 cache line utilization.
- **Delta Heuristic Scoring & O(1) Epoch Resets:**
  In `de_route_step` (`slow = false`), candidate SWAP/MOVE edges (u, v) are scored by computing distance deltas strictly over the active interaction lists of u and v rather than re-evaluating the entire front and lookahead layers. BFS lookahead traversal and edge deduplication avoid `std::memset` / `std::fill` overhead by comparing integer arrays (`visited_epoch`, `seen_edges`) against monotonically increasing step counters (`current_bfs_epoch`, `route_step_counter`) for O(1) state resets.
- **Small-Buffer Optimization & Avoidance of `std::vector<bool>`:**
  Topology and DAG nodes use `absl::InlinedVector` (`<int, 5>` for qubit-to-resonator links, `<int, 10>` for resonator-to-qubit links, `<int, 2>` for `HybridGateNode::next_gates`) to keep adjacency lists inline within the parent allocation. Bit-packed `std::vector<bool>` is replaced with byte-addressed `std::vector<uint8_t>` (`in_active`, `in_extended`, `seen_pairs`) to eliminate bit-masking instructions in tight loops.
- **Compile-Time Branch Elimination & Sub-Linear Crossover:**
  `evaluate_fitness_impl` and `de_force_down` are templated on `<bool LOGGING, bool CAN_LAZY>` using `if constexpr` to strip `RoutingRecord` (40-byte struct) construction and vector pushes from the evolutionary hot path, enabling logging only during final schedule reconstruction. Binomial crossover in `generate_trial_vector_windowed` uses `__restrict` pointers and geometric distribution skipping (`std::geometric_distribution<int>(cr)`) to jump directly between mutated gene indices in O(CR * W) time rather than iterating over the full window W.
- **Parallel Batch Evaluation & Qiskit Interoperability:**
  The batch evaluation harness parallelizes hyperparameter grid sweeps (`tune_hyperparameters`) and multi-circuit benchmarks via OpenMP (`#pragma omp parallel for schedule(dynamic, 1)`), using thread-local `std::ostringstream` buffers flushed under `#pragma omp critical`. A Python pipeline (`qiskitScript.py`) generates VQE, QAOA, Quantum Volume, and Shor proxy circuits, runs Qiskit's `sabre` transpiler (`optimization_level=0`, 20 trials), and exports JSON traces that C++ replays on the exact resonator timing model (`evaluate_qiskit_routed_circuit`) for direct baseline comparison.

## Benchmarks & Results

Benchmarks compare three routing pipelines evaluated under the identical resonator hardware timing model:
- `qiskit_sabre_export_eval`: Qiskit SABRE transpiler output replayed with resonator MOVE/SWAP timing.
- `sabre_fixed_layout`: Standalone bidirectional C++ `SabreRouter` (lookahead depth 20).
- `dehybrid_sabre_fixed_layout`: `SabreDEHybridRouter` with windowed DE (`pop_size = 30`, `window_size = 500`, `accept_ratio = 0.7`, `F = 0.9`, `CR = 0.9`, lookahead depth 4).

*Test environment: AMD RYZEN 7 5800x`, `16GB 3200MHz`, Windows x86_64, MSVC 14.44 (`/O2 /Ob2 /Oi /Ot /arch:AVX2 /fp:fast`, LTO enabled).*

| Topology | Qubits (Q) | Gates (G) | Router | Routing Cost | SWAPs / MOVEs | Makespan | Wall Time (s) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Triangle (3 Res) | 13 | 500 | Qiskit SABRE | `[INSERT]` | `[INSERT]` | `[INSERT]` | `[INSERT]` |
| Triangle (3 Res) | 13 | 500 | C++ SABRE | `[INSERT]` | `[INSERT]` | `[INSERT]` | `[INSERT]` |
| Triangle (3 Res) | 13 | 500 | **DE+SABRE Hybrid** | `[INSERT]` | `[INSERT]` | `[INSERT]` | `[INSERT]` |
| Grid 3x3 (9 Res) | 30 | 2000 | Qiskit SABRE | `[INSERT]` | `[INSERT]` | `[INSERT]` | `[INSERT]` |
| Grid 3x3 (9 Res) | 30 | 2000 | C++ SABRE | `[INSERT]` | `[INSERT]` | `[INSERT]` | `[INSERT]` |
| Grid 3x3 (9 Res) | 30 | 2000 | **DE+SABRE Hybrid** | `[INSERT]` | `[INSERT]` | `[INSERT]` | `[INSERT]` |
| Grid 5x5 (25 Res) | 70 | 2000 | Qiskit SABRE | `[INSERT]` | `[INSERT]` | `[INSERT]` | `[INSERT]` |
| Grid 5x5 (25 Res) | 70 | 2000 | C++ SABRE | `[INSERT]` | `[INSERT]` | `[INSERT]` | `[INSERT]` |
| Grid 5x5 (25 Res) | 70 | 2000 | **DE+SABRE Hybrid** | `[INSERT]` | `[INSERT]` | `[INSERT]` | `[INSERT]` |

### Performance Metrics Summary
- **Routing Cost Reduction:** `[INSERT: e.g., X.X% mean reduction]` in total routing cost (`SWAP_COST = 3.14`, `MOVE_COST = 1.0`) achieved by `SabreDEHybridRouter` over standalone SABRE across 216 benchmark circuits.
- **Makespan & Depth Improvement:** `[INSERT: e.g., X.X% reduction]` in physical execution makespan via concurrent resonator scheduling and lazy off-boarding.
- **Snapshot Windowing Speedup:** `[INSERT: e.g., X.Xx speedup]` in windowed DE evaluation time using `DEStateSnapshot` prefix resumption compared to full-circuit re-simulation.
- **Evaluation Throughput:** `[INSERT: e.g., X,XXX fitness evaluations/sec]` at 2,000 gates on a 70-qubit (5x5 resonator) topology.

## Build & Usage

### Prerequisites
- CMake >= 3.25
- C++20 compatible compiler (MSVC 19.30+, GCC 12+, or Clang 15+) with OpenMP support
- `vcpkg` (manifest mode dependencies: `abseil`, `boost-random`, `boost-align`, `fmt`, `nlohmann-json`)
- Python 3.10+ with `qiskit`, `pandas`, `numpy`, `scipy`, `seaborn`, `matplotlib` (for benchmark generation and plotting)

*(Note for Linux / case-sensitive filesystems: update `Generator/CircuitGenerator.*` to lowercase `generator/CircuitGenerator.*` in `CMakeLists.txt` before configuring).*

### Compilation

Configure and build in `Release` mode with Interprocedural Optimization (LTO) and AVX2 / native architecture flags enabled:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE="$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake"
cmake --build build --config Release
```

### Execution

1. **Generate Qiskit Benchmark Suite (Optional, for batch mode):**
   Generates 216 benchmark circuits (5 to 70 qubits, 50 to 8,000 gates across VQE, QAOA, Quantum Volume, and Shor styles) along with Qiskit SABRE baseline routes in `iqm_benchmarks/`:

```bash
python qiskitScript.py
```

2. **Run the C++ Optimizer:**
   By default, the active entry point generates a 2,000-gate circuit on a 5x5 resonator grid (70 physical qubits), executes `SabreDEHybridRouter::optimize_windowed`, validates the resulting schedule, and exports the physical operation log to `../hybrid_schedule_moves_Siatka 5x5.txt`:

```bash
./build/magisterka
```

3. **Generate Post-Processing Plots & LaTeX Tables:**
   Reads batch metrics from `results/routing_runs.csv` and outputs 300 DPI figures and LaTeX tables to `wykresy/`:

```bash
python latexDiagrams.py
```
