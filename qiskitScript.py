# import json
# import os
# import random
# from qiskit import QuantumCircuit, transpile
# from qiskit.circuit.library import quantum_volume
#
# class QiskitCppBridge:
#     def __init__(self):
#         self.num_qubits = 0
#
#     def export_to_cpp(self, qc: QuantumCircuit, filename: str):
#         self.num_qubits = qc.num_qubits
#         export_data = []
#
#         # Przechodzimy przez instrukcje obwodu
#         for gate_id, instr in enumerate(qc.data):
#             # Pobieramy indeksy kubitów, na których działa operacja
#             qubits = [qc.find_bit(q).index for q in instr.qubits]
#
#             # Wymuszamy format: TYLKO "1Q" lub "2Q" zgodnie z Twoim routerem
#             gate_type = "1Q" if len(qubits) == 1 else "2Q"
#
#             export_data.append({
#                 "id": gate_id,
#                 "type": gate_type,
#                 "qubits": qubits
#             })
#
#         with open(filename, "w") as f:
#             json.dump({"num_qubits": self.num_qubits, "gates": export_data}, f, indent=4)
#
# # --- TOPOLOGIA IQM CONSTELLATION ---
#
# def get_iqm_constellation_edges(num_qubits):
#     """
#     Generuje graf połączeń dla architektury IQM Constellation (Star Lattice).
#     Symuluje klastry 6-kubitowe (gwiazdy) wokół centralnego rezonatora.
#     """
#     edges = set()
#     clusters = []
#     current_q = 0
#
#     while current_q < num_qubits:
#         if not clusters:
#             cluster = list(range(min(6, num_qubits)))
#             clusters.append(cluster)
#             current_q = 6
#         else:
#             overlap = 2
#             shared = clusters[-1][-overlap:]
#             remaining_slots = 6 - overlap
#
#             if current_q >= num_qubits:
#                 break
#
#             new_qubits = list(range(current_q, min(current_q + remaining_slots, num_qubits)))
#             cluster = shared + new_qubits
#             clusters.append(cluster)
#             current_q += len(new_qubits)
#
#     for cluster in clusters:
#         for i in range(len(cluster)):
#             for j in range(i + 1, len(cluster)):
#                 q1, q2 = sorted([cluster[i], cluster[j]])
#                 edges.add((q1, q2))
#
#     return list(edges)
#
#
# # --- GENERATOR OBWODÓW ---
#
# def generate_circuit_for_style(num_qubits, target_gates, circuit_style):
#     """
#     Generuje czysto strukturalne obwody pod testy routingu.
#     Kąty i typy bramek nie mają znaczenia - służą tylko jako zaślepki 1Q i 2Q.
#     """
#     if circuit_style == "Constellation_VQE":
#         qc = QuantumCircuit(num_qubits)
#         hw_edges = get_iqm_constellation_edges(num_qubits)
#         if not hw_edges and num_qubits > 1:
#             hw_edges = [(0, 1)]
#
#         gates_per_layer = num_qubits + len(hw_edges)
#         num_layers = max(1, target_gates // gates_per_layer) if gates_per_layer > 0 else 1
#
#         for layer in range(num_layers):
#             for q in range(num_qubits):
#                 qc.rx(0, q)
#
#             selected_edges = hw_edges if layer % 2 == 0 else reversed(hw_edges)
#             for edge in selected_edges:
#                 qc.cx(edge[0], edge[1])
#
#         return qc
#
#     elif circuit_style == "Constellation_QAOA":
#         qc = QuantumCircuit(num_qubits)
#
#         problem_edges = []
#         for i in range(num_qubits):
#             for j in range(i + 1, num_qubits):
#                 if random.random() < 0.40:
#                     problem_edges.append((i, j))
#         if not problem_edges and num_qubits > 1:
#             problem_edges = [(0, 1)]
#
#         gates_per_layer = num_qubits + len(problem_edges)
#         num_layers = max(1, target_gates // gates_per_layer) if gates_per_layer > 0 else 1
#
#         for layer in range(num_layers):
#             for q in range(num_qubits):
#                 qc.rx(0, q)
#
#             for edge in problem_edges:
#                 qc.cx(edge[0], edge[1])
#
#         return qc
#
#     elif circuit_style == "QuantumVolume":
#         # Używamy nowej funkcji (quantum_volume) zamiast zdeprecjonowanej klasy
#         qc = quantum_volume(num_qubits, num_qubits)
#         # Baza 'u' i 'cx' gwarantuje, że Qiskit rozbije wszystko bez rzucania błędów
#         qc = transpile(qc, basis_gates=['u', 'cx'], optimization_level=0)
#         return qc
#
#     elif circuit_style == "Shor_Proxy":
#         qc = QuantumCircuit(num_qubits)
#         if num_qubits > 2:
#             counting_qubits = num_qubits // 2
#             target_qubits = num_qubits - counting_qubits
#
#             # Zaślepki 1Q
#             for i in range(counting_qubits):
#                 qc.rx(0, i)
#
#             # Splątanie między dwoma rejestrami
#             for i in range(counting_qubits):
#                 for j in range(target_qubits):
#                     qc.cx(i, counting_qubits + j)
#
#             # Ręczna struktura imitująca QFT (dużo krawędzi każdego z każdym)
#             # Używamy czystych bramek z Qiskita, żeby nie polegać na wycofanych klasach
#             for j in reversed(range(counting_qubits)):
#                 qc.h(j)
#                 for m in range(j):
#                     qc.cp(3.14 / (2 ** (j - m)), m, j)
#             for j in range(counting_qubits // 2):
#                 qc.swap(j, counting_qubits - j - 1)
#
#         # Qiskit teraz poprawnie rozbije bramki 'swap' i 'cp' na standardowe 'u' i 'cx'
#         qc = transpile(qc, basis_gates=['u', 'cx'], optimization_level=0)
#         return qc
#
#
# # --- GŁÓWNY SKRYPT ---
#
# def generate_iqm_benchmark_suite(bridge: QiskitCppBridge, base_dir: str = "iqm_benchmarks"):
#     QUBITS = [5, 10, 13, 15, 25, 30, 35, 38, 45, 55, 65, 70]
#     GATES = [50, 200, 500, 1000, 2000, 4000, 5000, 8000]
#
#     LAYERED_STYLES = ["Constellation_VQE", "Constellation_QAOA"]
#     ALGO_STYLES = ["QuantumVolume", "Shor_Proxy"]
#
#     print("🚀 Odpalam generowanie układów pod IQM Constellation...")
#
#     os.makedirs(base_dir, exist_ok=True)
#     total_exported = 0
#
#     for q in QUBITS:
#         q_dir = os.path.join(base_dir, f"{q}_qubits")
#         os.makedirs(q_dir, exist_ok=True)
#
#         # 1. Generowanie algorytmów z określoną liczbą bramek (VQE, QAOA)
#         for g in GATES:
#             g_dir = os.path.join(q_dir, f"{g}_gates")
#             os.makedirs(g_dir, exist_ok=True)
#
#             for style in LAYERED_STYLES:
#                 try:
#                     qc = generate_circuit_for_style(q, g, style)
#                     filename = f"{style}_{q}q_{g}g.json"
#                     filepath = os.path.join(g_dir, filename)
#
#                     bridge.export_to_cpp(qc, filepath)
#                     total_exported += 1
#                 except Exception as e:
#                     print(f"❌ Błąd przy generowaniu {style} dla {q}q i {g}g: {e}")
#
#         # 2. Generowanie algorytmów natywnych, gdzie rozrost zależy od liczby kubitów (QV, Shor)
#         algo_dir = os.path.join(q_dir, "algo_native_gates")
#         os.makedirs(algo_dir, exist_ok=True)
#
#         for style in ALGO_STYLES:
#             try:
#                 qc = generate_circuit_for_style(q, target_gates=None, circuit_style=style)
#                 actual_gates = len(qc.data) # Pobranie faktycznej wygenerowanej liczby bramek
#
#                 filename = f"{style}_{q}q_{actual_gates}g.json"
#                 filepath = os.path.join(algo_dir, filename)
#
#                 bridge.export_to_cpp(qc, filepath)
#                 total_exported += 1
#             except Exception as e:
#                 print(f"❌ Błąd przy generowaniu natywnego {style} dla {q}q: {e}")
#
#         print(f"✅ Przetworzono układ dla {q} kubitów.")
#
#     print("\n" + "="*60)
#     print(f"🔥 GOTOWE! Wyeksportowano {total_exported} plików do: {base_dir}")
#     print("="*60)
#
#
# if __name__ == "__main__":
#     bridge = QiskitCppBridge()
#     generate_iqm_benchmark_suite(bridge, base_dir="iqm_benchmarks")

import json
import os
import random
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import quantum_volume
from qiskit.transpiler import CouplingMap  # <-- DODANE: Wymagane dla Sabre

class QiskitCppBridge:
    def __init__(self):
        self.num_qubits = 0

    def export_to_cpp(self, qc: QuantumCircuit, filename: str):
        self.num_qubits = qc.num_qubits
        export_data = []

        for gate_id, instr in enumerate(qc.data):
            qubits = [qc.qubits.index(q) for q in instr.qubits]
            gate_type = "1Q" if len(qubits) == 1 else "2Q"

            export_data.append({
                "id": gate_id,
                "type": gate_type,
                "qubits": qubits
            })

        with open(filename, "w") as f:
            json.dump({"num_qubits": self.num_qubits, "gates": export_data}, f, indent=4)


# --- TOPOLOGIA IQM CONSTELLATION ---

def get_iqm_constellation_edges(num_qubits):
    edges = set()
    clusters = []
    current_q = 0

    while current_q < num_qubits:
        if not clusters:
            cluster = list(range(min(6, num_qubits)))
            clusters.append(cluster)
            current_q = 6
        else:
            overlap = 2
            shared = clusters[-1][-overlap:]
            remaining_slots = 6 - overlap

            if current_q >= num_qubits:
                break

            new_qubits = list(range(current_q, min(current_q + remaining_slots, num_qubits)))
            cluster = shared + new_qubits
            clusters.append(cluster)
            current_q += len(new_qubits)

    for cluster in clusters:
        for i in range(len(cluster)):
            for j in range(i + 1, len(cluster)):
                q1, q2 = sorted([cluster[i], cluster[j]])
                edges.add((q1, q2))

    return list(edges)


# --- GENERATOR OBWODÓW ---

def generate_circuit_for_style(num_qubits, target_gates, circuit_style, hw_edges=None):
    if circuit_style == "Constellation_VQE":
        qc = QuantumCircuit(num_qubits)
        valid_edges = [e for e in hw_edges if e[0] < num_qubits and e[1] < num_qubits] if hw_edges else []
        if not valid_edges and num_qubits > 1:
            valid_edges = [(0, 1)]

        gates_per_layer = num_qubits + len(valid_edges)
        num_layers = max(1, target_gates // gates_per_layer) if gates_per_layer > 0 else 1

        for layer in range(num_layers):
            for q in range(num_qubits):
                qc.rx(1.57, q)  # ZMIANA: Dodany niezerowy kąt (np. pi/2)
            selected_edges = valid_edges if layer % 2 == 0 else reversed(valid_edges)
            for edge in selected_edges:
                qc.cx(edge[0], edge[1])

        return qc

    elif circuit_style == "Constellation_QAOA":
        qc = QuantumCircuit(num_qubits)
        problem_edges = []
        for i in range(num_qubits):
            for j in range(i + 1, num_qubits):
                if random.random() < 0.40:
                    problem_edges.append((i, j))
        if not problem_edges and num_qubits > 1:
            problem_edges = [(0, 1)]

        gates_per_layer = num_qubits + len(problem_edges)
        num_layers = max(1, target_gates // gates_per_layer) if gates_per_layer > 0 else 1

        for layer in range(num_layers):
            for q in range(num_qubits):
                qc.rx(1.57, q)  # ZMIANA: Dodany niezerowy kąt (np. pi/2)
            for edge in problem_edges:
                qc.cx(edge[0], edge[1])

        return qc

    elif circuit_style == "QuantumVolume":
        qc = quantum_volume(num_qubits, num_qubits)
        qc = transpile(qc, basis_gates=['u', 'cx'], optimization_level=0)
        return qc

    elif circuit_style == "Shor_Proxy":
        qc = QuantumCircuit(num_qubits)
        if num_qubits > 2:
            counting_qubits = num_qubits // 2
            target_qubits = num_qubits - counting_qubits

            for i in range(counting_qubits):
                qc.rx(1.57, i) # <--- TUTAJ ZMIANA Z 0 NA 1.57
            for i in range(counting_qubits):
                for j in range(target_qubits):
                    qc.cx(i, counting_qubits + j)

            for j in reversed(range(counting_qubits)):
                qc.h(j)
                for m in range(j):
                    qc.cp(3.14 / (2 ** (j - m)), m, j)
            for j in range(counting_qubits // 2):
                qc.swap(j, counting_qubits - j - 1)

        qc = transpile(qc, basis_gates=['u', 'cx'], optimization_level=0)
        return qc


# --- GŁÓWNY SKRYPT ---

import json
import os
import random
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import quantum_volume
from qiskit.transpiler import CouplingMap

class QiskitCppBridge:
    def __init__(self):
        self.num_qubits = 0

    def export_to_cpp(self, qc: QuantumCircuit, filename: str):
        self.num_qubits = qc.num_qubits
        export_data = []

        for gate_id, instr in enumerate(qc.data):
            qubits = [qc.qubits.index(q) for q in instr.qubits]
            gate_type = "1Q" if len(qubits) == 1 else "2Q"

            export_data.append({
                "id": gate_id,
                "type": gate_type,
                "qubits": qubits
            })

        with open(filename, "w") as f:
            json.dump({"num_qubits": self.num_qubits, "gates": export_data}, f, indent=4)

    def export_transpiled_to_cpp(self, qc: QuantumCircuit, filename: str):
        self.num_qubits = qc.num_qubits
        export_data = []

        for gate_id, instr in enumerate(qc.data):
            qubits = [qc.qubits.index(q) for q in instr.qubits]
            op_name = instr.operation.name.lower()

            if op_name == "swap":
                gate_type = "SWAP"
            elif len(qubits) == 1:
                gate_type = "1Q"
            else:
                gate_type = "2Q"

            export_data.append({
                "id": gate_id,
                "type": gate_type,
                "qubits": qubits
            })

        with open(filename, "w") as f:
            json.dump({"num_qubits": self.num_qubits, "gates": export_data}, f, indent=4)


# =====================================================================
# --- DOKŁADNE ODWZOROWANIE TOPOLOGII C++ (ŹRÓDŁO PRAWDY DLA QISKITA) ---
# =====================================================================

def _build_edges_from_resonators(resonator_coords):
    dy = [0, -1, -1, 0, 1, 1]
    dx = [1, 1, 0, -1, -1, 0]

    qubit_map = {}
    qubits_list = []
    resonators_connected_qubits = {r_id: [] for r_id in resonator_coords.values()}

    # 2. Budowa połączeń (dokładnie jak w C++ build_topology)
    for (y, x), r_id in resonator_coords.items():
        for i in range(6):
            n1_coord = (y + dy[i], x + dx[i])
            next_i = (i + 1) % 6
            n2_coord = (y + dy[next_i], x + dx[next_i])

            n1_id = resonator_coords.get(n1_coord, -1)
            n2_id = resonator_coords.get(n2_coord, -1)

            res_group = {r_id}
            if n1_id != -1: res_group.add(n1_id)
            if n2_id != -1: res_group.add(n2_id)

            res_group_frozenset = frozenset(res_group)

            if res_group_frozenset not in qubit_map:
                new_q_id = len(qubits_list)
                qubits_list.append(res_group)
                qubit_map[res_group_frozenset] = new_q_id

            q_id = qubit_map[res_group_frozenset]
            if q_id not in resonators_connected_qubits[r_id]:
                resonators_connected_qubits[r_id].append(q_id)

    # 3. Dopełnianie (padding) rezonatorów do 6 kubitów (brzegowe)
    for r_id in resonator_coords.values():
        while len(resonators_connected_qubits[r_id]) < 6:
            new_q_id = len(qubits_list)
            qubits_list.append({r_id})
            resonators_connected_qubits[r_id].append(new_q_id)

    # 4. Generowanie krawędzi dla Qiskita: dwa kubity mogą zrobić CX, jeśli dzielą rezonator
    edges = set()
    for connected_qubits in resonators_connected_qubits.values():
        for i in range(len(connected_qubits)):
            for j in range(i + 1, len(connected_qubits)):
                q1, q2 = sorted([connected_qubits[i], connected_qubits[j]])
                edges.add((q1, q2))

    return list(edges), len(qubits_list)


def build_iqm_constellation_edges(Y_rows: int, X_cols: int):
    resonator_coords = {}
    for y in range(Y_rows):
        for x in range(X_cols):
            r_id = y * X_cols + x
            resonator_coords[(y, x)] = r_id
    return _build_edges_from_resonators(resonator_coords)


def build_iqm_triangle_edges():
    resonator_coords = {
        (0, 0): 0,
        (0, 1): 1,
        (-1, 1): 2
    }
    return _build_edges_from_resonators(resonator_coords)


def get_iqm_topology_edges_and_size(target_logical_qubits):
    """
    Replikuje C++ main.cpp:
    - Dla <= 13 kubitów wybiera trójkąt.
    - Dla > 13 skaluje macierz dim x dim aż znajdzie wystarczającą pojemność.
    """
    if target_logical_qubits <= 13:
        return build_iqm_triangle_edges()
    else:
        dim = 2
        edges, phys_q_count = build_iqm_constellation_edges(dim, dim)
        while phys_q_count < target_logical_qubits:
            dim += 1
            edges, phys_q_count = build_iqm_constellation_edges(dim, dim)
        return edges, phys_q_count


# =====================================================================
# --- GENERATOR OBWODÓW ---
# =====================================================================

def generate_circuit_for_style(num_qubits, target_gates, circuit_style, hw_edges=None):
    if circuit_style == "Constellation_VQE":
        qc = QuantumCircuit(num_qubits)
        # Bierzemy hw_edges z faktycznej topologii, przycinamy do rozmiaru logicznego
        valid_edges = [e for e in hw_edges if e[0] < num_qubits and e[1] < num_qubits] if hw_edges else []
        if not valid_edges and num_qubits > 1:
            valid_edges = [(0, 1)]

        gates_per_layer = num_qubits + len(valid_edges)
        num_layers = max(1, target_gates // gates_per_layer) if gates_per_layer > 0 else 1

        for layer in range(num_layers):
            for q in range(num_qubits):
                qc.rx(0, q)
            selected_edges = valid_edges if layer % 2 == 0 else reversed(valid_edges)
            for edge in selected_edges:
                qc.cx(edge[0], edge[1])

        return qc

    elif circuit_style == "Constellation_QAOA":
        qc = QuantumCircuit(num_qubits)
        problem_edges = []
        for i in range(num_qubits):
            for j in range(i + 1, num_qubits):
                if random.random() < 0.40:
                    problem_edges.append((i, j))
        if not problem_edges and num_qubits > 1:
            problem_edges = [(0, 1)]

        gates_per_layer = num_qubits + len(problem_edges)
        num_layers = max(1, target_gates // gates_per_layer) if gates_per_layer > 0 else 1

        for layer in range(num_layers):
            for q in range(num_qubits):
                qc.rx(0, q)
            for edge in problem_edges:
                qc.cx(edge[0], edge[1])

        return qc

    elif circuit_style == "QuantumVolume":
        qc = quantum_volume(num_qubits, num_qubits)
        qc = transpile(qc, basis_gates=['u', 'cx'], optimization_level=0)
        return qc

    elif circuit_style == "Shor_Proxy":
        qc = QuantumCircuit(num_qubits)
        if num_qubits > 2:
            counting_qubits = num_qubits // 2
            target_qubits = num_qubits - counting_qubits

            for i in range(counting_qubits):
                qc.rx(0, i)
            for i in range(counting_qubits):
                for j in range(target_qubits):
                    qc.cx(i, counting_qubits + j)

            for j in reversed(range(counting_qubits)):
                qc.h(j)
                for m in range(j):
                    qc.cp(3.14 / (2 ** (j - m)), m, j)
            for j in range(counting_qubits // 2):
                qc.swap(j, counting_qubits - j - 1)

        qc = transpile(qc, basis_gates=['u', 'cx'], optimization_level=0)
        return qc


# =====================================================================
# --- GŁÓWNY SKRYPT ---
# =====================================================================

def generate_iqm_benchmark_suite(bridge: QiskitCppBridge, base_dir: str = "iqm_benchmarks"):
    QUBITS = [5, 10, 13, 15, 25, 30, 35, 38, 45, 55, 65, 70]
    GATES = [50, 200, 500, 1000, 2000, 4000, 5000, 8000]

    LAYERED_STYLES = ["Constellation_VQE", "Constellation_QAOA"]
    ALGO_STYLES = ["QuantumVolume", "Shor_Proxy"]

    print("🚀 Odpalam generowanie układów pod IQM Constellation i estymację Sabre...\n")

    os.makedirs(base_dir, exist_ok=True)
    total_exported = 0

    for q in QUBITS:
        # Generujemy lustrzane krawędzie i pobieramy FAKTYCZNY fizyczny limit
        hw_edges, phys_qubits_count = get_iqm_topology_edges_and_size(q)

        print(f"--- 🖥️  KLASTER: {q} KUBITÓW LOGICZNYCH (Hardware ma {phys_qubits_count} kubitów fizycznych) ---")
        q_dir = os.path.join(base_dir, f"{q}_qubits")
        os.makedirs(q_dir, exist_ok=True)

        cmap = CouplingMap(hw_edges)
        cmap.make_symmetric()

        for g in GATES:
            g_dir = os.path.join(q_dir, f"{g}_gates")
            os.makedirs(g_dir, exist_ok=True)

            for style in LAYERED_STYLES:
                try:
                    qc = generate_circuit_for_style(q, g, style, hw_edges)

                    filename = f"{style}_{q}q_{g}g.json"
                    filepath = os.path.join(g_dir, filename)
                    bridge.export_to_cpp(qc, filepath)
                    total_exported += 1

                    qc_sabre = transpile(qc, coupling_map=cmap, routing_method='sabre', optimization_level=0, basis_gates=['u', 'cx', 'swap'])
                    swaps_added = qc_sabre.count_ops().get('swap', 0)

                    filepath_qiskit = os.path.join(g_dir, f"{style}_{q}q_{g}g_qiskit_routed.json")
                    bridge.export_transpiled_to_cpp(qc_sabre, filepath_qiskit)

                    print(f"  ├── [ {q}Q | {g}G | {style} ] Sabre dodał: {swaps_added} SWAP")

                except Exception as e:
                    print(f"  ├── ❌ Błąd przy generowaniu {style} dla {q}q i {g}g: {e}")

        algo_dir = os.path.join(q_dir, "algo_native_gates")
        os.makedirs(algo_dir, exist_ok=True)

        for style in ALGO_STYLES:
            try:
                qc = generate_circuit_for_style(q, target_gates=None, circuit_style=style)
                actual_gates = len(qc.data)

                filename = f"{style}_{q}q_{actual_gates}g.json"
                filepath = os.path.join(algo_dir, filename)
                bridge.export_to_cpp(qc, filepath)
                total_exported += 1

                qc_sabre = transpile(qc, coupling_map=cmap, routing_method='sabre', optimization_level=0, basis_gates=['u', 'cx', 'swap'])
                swaps_added = qc_sabre.count_ops().get('swap', 0)

                filepath_qiskit = os.path.join(algo_dir, f"{style}_{q}q_{actual_gates}g_qiskit_routed.json")
                bridge.export_transpiled_to_cpp(qc_sabre, filepath_qiskit)

                print(f"  └── [ {q}Q | {actual_gates}G (Natywne) | {style} ] Sabre dodał: {swaps_added} SWAP")

            except Exception as e:
                print(f"  └── ❌ Błąd przy generowaniu natywnego {style} dla {q}q: {e}")

        print()

    print("="*60)
    print(f"🔥 GOTOWE! Wyeksportowano {total_exported} plików bazowych (i ich wersje qiskit_routed) do: {base_dir}")
    print("="*60)


if __name__ == "__main__":
    bridge = QiskitCppBridge()
    generate_iqm_benchmark_suite(bridge, base_dir="iqm_benchmarks")