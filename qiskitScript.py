import json
import os
import random
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import quantum_volume

class QiskitCppBridge:
    def __init__(self):
        self.num_qubits = 0

    def export_to_cpp(self, qc: QuantumCircuit, filename: str):
        self.num_qubits = qc.num_qubits
        export_data = []

        # Przechodzimy przez instrukcje obwodu
        for gate_id, instr in enumerate(qc.data):
            # Pobieramy indeksy kubitów, na których działa operacja
            qubits = [qc.find_bit(q).index for q in instr.qubits]

            # Wymuszamy format: TYLKO "1Q" lub "2Q" zgodnie z Twoim routerem
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
    """
    Generuje graf połączeń dla architektury IQM Constellation (Star Lattice).
    Symuluje klastry 6-kubitowe (gwiazdy) wokół centralnego rezonatora.
    """
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

def generate_circuit_for_style(num_qubits, target_gates, circuit_style):
    """
    Generuje czysto strukturalne obwody pod testy routingu.
    Kąty i typy bramek nie mają znaczenia - służą tylko jako zaślepki 1Q i 2Q.
    """
    if circuit_style == "Constellation_VQE":
        qc = QuantumCircuit(num_qubits)
        hw_edges = get_iqm_constellation_edges(num_qubits)
        if not hw_edges and num_qubits > 1:
            hw_edges = [(0, 1)]

        gates_per_layer = num_qubits + len(hw_edges)
        num_layers = max(1, target_gates // gates_per_layer) if gates_per_layer > 0 else 1

        for layer in range(num_layers):
            for q in range(num_qubits):
                qc.rx(0, q)

            selected_edges = hw_edges if layer % 2 == 0 else reversed(hw_edges)
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
        # Używamy nowej funkcji (quantum_volume) zamiast zdeprecjonowanej klasy
        qc = quantum_volume(num_qubits, num_qubits)
        # Baza 'u' i 'cx' gwarantuje, że Qiskit rozbije wszystko bez rzucania błędów
        qc = transpile(qc, basis_gates=['u', 'cx'], optimization_level=0)
        return qc

    elif circuit_style == "Shor_Proxy":
        qc = QuantumCircuit(num_qubits)
        if num_qubits > 2:
            counting_qubits = num_qubits // 2
            target_qubits = num_qubits - counting_qubits

            # Zaślepki 1Q
            for i in range(counting_qubits):
                qc.rx(0, i)

            # Splątanie między dwoma rejestrami
            for i in range(counting_qubits):
                for j in range(target_qubits):
                    qc.cx(i, counting_qubits + j)

            # Ręczna struktura imitująca QFT (dużo krawędzi każdego z każdym)
            # Używamy czystych bramek z Qiskita, żeby nie polegać na wycofanych klasach
            for j in reversed(range(counting_qubits)):
                qc.h(j)
                for m in range(j):
                    qc.cp(3.14 / (2 ** (j - m)), m, j)
            for j in range(counting_qubits // 2):
                qc.swap(j, counting_qubits - j - 1)

        # Qiskit teraz poprawnie rozbije bramki 'swap' i 'cp' na standardowe 'u' i 'cx'
        qc = transpile(qc, basis_gates=['u', 'cx'], optimization_level=0)
        return qc


# --- GŁÓWNY SKRYPT ---

def generate_iqm_benchmark_suite(bridge: QiskitCppBridge, base_dir: str = "iqm_benchmarks"):
    QUBITS = [5, 10, 13, 15, 25, 30, 35, 38, 45, 55, 65, 70]
    GATES = [50, 200, 500, 1000, 2000, 4000, 5000, 8000]

    LAYERED_STYLES = ["Constellation_VQE", "Constellation_QAOA"]
    ALGO_STYLES = ["QuantumVolume", "Shor_Proxy"]

    print("🚀 Odpalam generowanie układów pod IQM Constellation...")

    os.makedirs(base_dir, exist_ok=True)
    total_exported = 0

    for q in QUBITS:
        q_dir = os.path.join(base_dir, f"{q}_qubits")
        os.makedirs(q_dir, exist_ok=True)

        # 1. Generowanie algorytmów z określoną liczbą bramek (VQE, QAOA)
        for g in GATES:
            g_dir = os.path.join(q_dir, f"{g}_gates")
            os.makedirs(g_dir, exist_ok=True)

            for style in LAYERED_STYLES:
                try:
                    qc = generate_circuit_for_style(q, g, style)
                    filename = f"{style}_{q}q_{g}g.json"
                    filepath = os.path.join(g_dir, filename)

                    bridge.export_to_cpp(qc, filepath)
                    total_exported += 1
                except Exception as e:
                    print(f"❌ Błąd przy generowaniu {style} dla {q}q i {g}g: {e}")

        # 2. Generowanie algorytmów natywnych, gdzie rozrost zależy od liczby kubitów (QV, Shor)
        algo_dir = os.path.join(q_dir, "algo_native_gates")
        os.makedirs(algo_dir, exist_ok=True)

        for style in ALGO_STYLES:
            try:
                qc = generate_circuit_for_style(q, target_gates=None, circuit_style=style)
                actual_gates = len(qc.data) # Pobranie faktycznej wygenerowanej liczby bramek

                filename = f"{style}_{q}q_{actual_gates}g.json"
                filepath = os.path.join(algo_dir, filename)

                bridge.export_to_cpp(qc, filepath)
                total_exported += 1
            except Exception as e:
                print(f"❌ Błąd przy generowaniu natywnego {style} dla {q}q: {e}")

        print(f"✅ Przetworzono układ dla {q} kubitów.")

    print("\n" + "="*60)
    print(f"🔥 GOTOWE! Wyeksportowano {total_exported} plików do: {base_dir}")
    print("="*60)


if __name__ == "__main__":
    bridge = QiskitCppBridge()
    generate_iqm_benchmark_suite(bridge, base_dir="iqm_benchmarks")

# import json
# import os
# import random
# from qiskit import QuantumCircuit, transpile
# from qiskit.circuit.library import QuantumVolume, QFT
#
# class QiskitCppBridge:
#     def __init__(self):
#         self.num_qubits = 0
#
#     def export_to_cpp(self, qc: QuantumCircuit, filename: str):
#         self.num_qubits = qc.num_qubits
#         export_data = []
#
#         for gate_id, instr in enumerate(qc.data):
#             qubits = [qc.find_bit(q).index for q in instr.qubits]
#             # Router wymaga tylko 1Q i 2Q - wymuszamy to bez względu na nazwę operacji
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
# # --- FUNKCJE POMOCNICZE (bez zmian) ---
# def get_iqm_constellation_edges(num_qubits):
#     edges = set()
#     clusters = []
#     current_q = 0
#     while current_q < num_qubits:
#         if not clusters:
#             cluster = list(range(min(6, num_qubits)))
#             clusters.append(cluster)
#             current_q = 6
#         else:
#             overlap = 2
#             shared = clusters[-1][-overlap:]
#             remaining_slots = 6 - overlap
#             if current_q >= num_qubits: break
#             new_qubits = list(range(current_q, min(current_q + remaining_slots, num_qubits)))
#             cluster = shared + new_qubits
#             clusters.append(cluster)
#             current_q += len(new_qubits)
#     for cluster in clusters:
#         for i in range(len(cluster)):
#             for j in range(i + 1, len(cluster)):
#                 edges.add(tuple(sorted([cluster[i], cluster[j]])))
#     return list(edges)
#
# def generate_circuit_for_style(num_qubits, target_gates, circuit_style):
#     # WERSJA NIEZOPTYMALIZOWANA (tylko rozbicie na basis_gates)
#     if "opt" not in circuit_style:
#         level = 0
#     else:
#         level = 3 # Qiskit robi co może, by zredukować liczbę bramek
#
#     clean_style = circuit_style.replace("_opt", "")
#
#     if clean_style in ["Constellation_VQE", "Constellation_QAOA"]:
#         qc = QuantumCircuit(num_qubits)
#         gates_added = 0
#         edges = get_iqm_constellation_edges(num_qubits)
#         layer_idx = 0
#         while gates_added < target_gates:
#             for q in range(num_qubits):
#                 if gates_added >= target_gates: break
#                 qc.ry(0.5, q) if layer_idx % 2 == 0 else qc.rx(0.5, q)
#                 gates_added += 1
#             for edge in edges:
#                 if gates_added >= target_gates: break
#                 qc.cx(edge[0], edge[1])
#                 gates_added += 1
#             layer_idx += 1
#         return transpile(qc, basis_gates=['rx', 'ry', 'rz', 'cx'], optimization_level=level)
#
#     elif clean_style == "QuantumVolume":
#         qc = QuantumVolume(num_qubits, num_qubits)
#         return transpile(qc, basis_gates=['rx', 'ry', 'rz', 'cx'], optimization_level=level)
#
#     elif clean_style == "Shor_Proxy":
#         qc = QuantumCircuit(num_qubits)
#         if num_qubits > 2:
#             mid = num_qubits // 2
#             for i in range(mid): qc.h(i)
#             for i in range(mid):
#                 for j in range(num_qubits - mid):
#                     qc.cx(i, mid + j)
#             qc.append(QFT(mid, inverse=True), range(mid))
#         return transpile(qc, basis_gates=['rx', 'ry', 'rz', 'cx'], optimization_level=level)
#
# # --- GŁÓWNA FUNKCJA ---
# def generate_iqm_benchmark_suite(bridge: QiskitCppBridge, base_dir: str = "iqm_benchmarks"):
#     QUBITS = [5, 10, 13] # Zmniejsz dla testu, potem zwiększ
#     GATES = [50, 200]
#
#     # Dodajemy wersje "_opt"
#     STYLES = ["Constellation_VQE", "Constellation_QAOA", "QuantumVolume", "Shor_Proxy"]
#
#     os.makedirs(base_dir, exist_ok=True)
#
#     for q in QUBITS:
#         q_dir = os.path.join(base_dir, f"{q}_qubits")
#         os.makedirs(q_dir, exist_ok=True)
#
#         for style in STYLES:
#             # Generujemy wersję "surową" i "zoptymalizowaną"
#             for s_name in [style, f"{style}_opt"]:
#                 try:
#                     # Dla algorytmów algorytmicznych (QV/Shor) target_gates ignorujemy
#                     target = 500 if "Constellation" in s_name else 0
#                     qc = generate_circuit_for_style(q, target, s_name)
#
#                     filename = f"{s_name}_{q}q_{len(qc.data)}g.json"
#                     bridge.export_to_cpp(qc, os.path.join(q_dir, filename))
#                     print(f"✅ Wygenerowano: {filename}")
#                 except Exception as e:
#                     print(f"❌ Błąd {s_name} ({q}q): {e}")
#
# if __name__ == "__main__":
#     generate_iqm_benchmark_suite(QiskitCppBridge())