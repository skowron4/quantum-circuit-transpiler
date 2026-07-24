import pandas as pd
import seaborn as sns
import math
import numpy as np
import os
import re

def sanitize_filename(name):
    return re.sub(r'[^\w\-]', '_', str(name))

# =====================================================================
# NOWE FUNKCJE DO ZAGREGOWANYCH TABEL LATEX
# =====================================================================

def write_legend_table(circuit_map, filepath):
    """Tworzy globalną legendę w ułożeniu dwukolumnowym (ID | Nazwa | ID | Nazwa)."""
    # Sortujemy po numerze C
    sorted_circs = sorted(circuit_map.keys(), key=lambda c: int(circuit_map[c][1:]))
    n = len(sorted_circs)

    # Ustawiamy liczbę "dużych" kolumn (2 oznacza: ID|Nazwa || ID|Nazwa)
    num_cols = 2
    rows = math.ceil(n / num_cols)

    lines = [
        "\\begingroup",
        "\\footnotesize",
        "\\renewcommand{\\arraystretch}{0.95} % Lekkie zwężenie wierszy w pionie",
        "\\begin{longtable}{|c|l|c|l|}",
        "\\caption{Legenda badanych obwodów kwantowych}\\label{tab:legenda_obwodow} \\\\",
        "\\hline",
        "\\textbf{ID} & \\textbf{Oryginalna nazwa obwodu} & \\textbf{ID} & \\textbf{Oryginalna nazwa obwodu} \\\\",
        "\\hline",
        "\\endfirsthead",
        "\\multicolumn{4}{c}{\\tablename\\ \\thetable{} -- ciąg dalszy z poprzedniej strony} \\\\",
        "\\hline",
        "\\textbf{ID} & \\textbf{Oryginalna nazwa obwodu} & \\textbf{ID} & \\textbf{Oryginalna nazwa obwodu} \\\\",
        "\\hline",
        "\\endhead",
        "\\hline",
        "\\endfoot",
        "\\hline",
        "\\endlastfoot"
    ]

    # Generowanie wierszy dla wielokolumnowego układu
    for i in range(rows):
        row_items = []
        for j in range(num_cols):
            idx = i + j * rows
            if idx < n:
                circ = sorted_circs[idx]
                c_id = circuit_map[circ]
                c_name_safe = circ.replace('_', '\\_')
                row_items.append(f"\\textbf{{{c_id}}} & {c_name_safe}")
            else:
                row_items.append(" & ") # Puste komórki, jeśli brakuje elementów do pary

        lines.append(" & ".join(row_items) + " \\\\")

    lines.extend([
        "\\end{longtable}",
        "\\endgroup",
        "\\vspace{1em}",
        "% ---------------------------------------------------------\n\n"
    ])

    with open(filepath, 'a', encoding='utf-8') as f:
        f.write("\n".join(lines))

def append_aggregated_tables(df_topo, filepath, layout, topo, circuit_map):
    routers = [r for r in ['Sabre', 'DE+SABRE', 'Qiskit'] if r in df_topo['router_short'].unique()]
    if not routers: return

    df_mean = df_topo.groupby(['circuit_clean', 'router_short'], observed=True)[['swaps', 'moves', 'weighted_cost', 'execution_time_s']].mean().unstack()
    circuits_in_topo = sorted(df_topo['circuit_clean'].dropna().unique(), key=lambda c: int(circuit_map[c][1:]))

    comparisons = []
    if 'DE+SABRE' in routers and 'Sabre' in routers: comparisons.append(('DE+SABRE', 'Sabre'))
    if 'DE+SABRE' in routers and 'Qiskit' in routers: comparisons.append(('DE+SABRE', 'Qiskit'))

    lines = []
    is_wide = 'Qiskit' in routers

    lines.append("\\begingroup")
    if is_wide: lines.append("\\begin{landscape}")

    lines.append("\\footnotesize")
    lines.append("\\renewcommand{\\arraystretch}{1.1}")
    lines.append("\\setlength{\\tabcolsep}{2pt}")

    # ZMIANA: Usunięto 1 kolumnę z każdego bloku porównania (3 kolumny zamiast 4)
    blocks = ["cccc" for _ in routers] + ["ccc" for _ in comparisons]
    col_format = "|c|" + "|".join(blocks) + "|"

    lines.append(f"\\begin{{longtable}}{{{col_format}}}")
    lines.append(f"\\caption{{Wyniki i porównanie: {topo.replace('_', '\\_')}}} \\\\")
    lines.append("\\hline")

    h1 = ["\\textbf{ID}"]
    for r in routers: h1.append(f"\\multicolumn{{4}}{{c|}}{{\\textbf{{{r}}}}}")
    for i, c in enumerate(comparisons):
        pipe = "|" if i == len(comparisons)-1 else "|"
        h1.append(f"\\multicolumn{{3}}{{c{pipe}}}{{\\textbf{{{c[0]} vs {c[1]}}}}}")
    lines.append(" & ".join(h1) + " \\\\")

    # ZMIANA: Nagłówek dla porównania bez '% t'
    h2 = [""] + ["S & M & K & t"] * len(routers) + ["\\% S & \\% M & \\% K"] * len(comparisons)
    lines.append(" & ".join(h2) + " \\\\")
    lines.append("\\hline\\endfirsthead")
    lines.append(f"\\multicolumn{{{1 + 4*len(routers) + 3*len(comparisons)}}}{{c}}{{-- cd. --}} \\\\ \\hline")
    lines.append(" & ".join(h1) + " \\\\")
    lines.append(" & ".join(h2) + " \\\\ \\hline\\endhead\\hline\\endfoot\\endlastfoot")

    # Funkcja do obliczeń: 0 vs 0 to "-", procenty całkowite
    def calc_pct(a, b):
        if pd.isna(a) or pd.isna(b): return "-"
        if abs(a - b) < 0.001: return "-"  # ZMIANA: 0 vs 0 to kreska
        if b < 0.01: return "$\\infty$"
        val = ((a - b) / b) * 100
        return ">9999" if abs(val) > 9999 else f"{int(round(val)):+d}"

    for circ in circuits_in_topo:
        c_id = circuit_map[circ]
        row = [f"\\textbf{{{c_id}}}"]
        # Wyniki surowe
        for r in routers:
            if r in df_mean['swaps'].columns and pd.notna(df_mean['swaps'][r].get(circ)):
                s, m, k, t = df_mean['swaps'][r][circ], df_mean['moves'][r][circ], df_mean['weighted_cost'][r][circ], df_mean['execution_time_s'][r][circ]
                row.extend([str(int(s)), str(int(m)), f"{k:.1f}", "$\\infty$" if t == np.inf else f"{t:.2f}"])
            else: row.extend(["-", "-", "-", "-"])

        # Porównanie bez czasu
        for c_A, c_B in comparisons:
            s_A, s_B = df_mean['swaps'][c_A].get(circ, np.nan), df_mean['swaps'][c_B].get(circ, np.nan)
            m_A, m_B = df_mean['moves'][c_A].get(circ, np.nan), df_mean['moves'][c_B].get(circ, np.nan)
            k_A, k_B = df_mean['weighted_cost'][c_A].get(circ, np.nan), df_mean['weighted_cost'][c_B].get(circ, np.nan)
            row.extend([calc_pct(s_A, s_B), calc_pct(m_A, m_B), calc_pct(k_A, k_B)])

        lines.append(" & ".join(row) + " \\\\")

    lines.append("\\hline\\end{longtable}")
    if is_wide: lines.append("\\end{landscape}")
    lines.append("\\endgroup\n\\vspace{1em}\n")

    with open(filepath, 'a', encoding='utf-8') as f:
        f.write("\n".join(lines))

# =====================================================================
# ZMODYFIKOWANA GŁÓWNA FUNKCJA KONTROLNA
# =====================================================================

def getDiagrams():
    csv_path = 'results/routing_runs.csv'
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"❌ Nie znaleziono pliku: {csv_path}")
        return

    cols_to_numeric = ['logical_qubits', 'gates', 'swaps', 'moves', 'execution_time_s', 'makespan']
    for col in cols_to_numeric:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    df = df[df['gates'] <= 3500].copy()

    router_mapping = {
        'sabre_fixed_layout': 'Sabre',
        'dehybrid_sabre_fixed_layout': 'DE+SABRE',
        'qiskit_sabre_export_eval': 'Qiskit'
    }
    df['router_short'] = df['router'].replace(router_mapping)

    qiskit_mask = df['router'] == 'qiskit_sabre_export_eval'
    df.loc[qiskit_mask, 'layout_strategy'] = 'qiskit_sabre_initial_layout'
    df['circuit_clean'] = df['circuit_id'].str.replace('Constellation_', '').str.replace('_Proxy', '')

    # =================================================================
    # FILTROWANIE OBWODÓW
    # =================================================================
    c1_qubits = [25, 30, 35, 38, 45, 55, 65, 70]
    m1 = df['logical_qubits'].isin(c1_qubits) & ((df['gates'] == 50) | df['circuit_clean'].str.contains('_50g|50g', regex=True))

    c2_qubits = [30, 35, 38, 45, 55, 65, 70]
    m2 = df['logical_qubits'].isin(c2_qubits) & df['circuit_clean'].str.contains('QAOA', case=False) & ((df['gates'] == 200) | df['circuit_clean'].str.contains('_200g|200g', regex=True))

    c3_qubits = [55, 65, 70]
    m3_vqe = df['logical_qubits'].isin(c3_qubits) & df['circuit_clean'].str.contains('VQE', case=False) & ((df['gates'] == 200) | df['circuit_clean'].str.contains('_200g|200g', regex=True))
    m3_qaoa = df['logical_qubits'].isin(c3_qubits) & df['circuit_clean'].str.contains('QAOA', case=False) & ((df['gates'] == 500) | df['circuit_clean'].str.contains('_500g|500g', regex=True))
    m3 = m3_vqe | m3_qaoa

    c4_qubits = [70]
    m4 = df['logical_qubits'].isin(c4_qubits) & df['circuit_clean'].str.contains('QAOA', case=False) & ((df['gates'] == 2000) | df['circuit_clean'].str.contains('_2000g|2000g', regex=True))

    df = df[~(m1 | m2 | m3 | m4)].copy()
    # =================================================================

    df.sort_values(by=['logical_qubits', 'gates'], inplace=True)
    df['swaps'] = df['swaps'].fillna(0)
    df['moves'] = df['moves'].fillna(0)
    df['weighted_cost'] = (df['swaps'] * 3.14) + (df['moves'] * 1.0)

    # 1. GENEROWANIE GLOBALNEGO SŁOWNIKA IDENTYFIKATORÓW (C1, C2...)
    circuit_order_df = df[['circuit_clean', 'logical_qubits', 'gates']].dropna().drop_duplicates(subset=['circuit_clean'])
    circuit_order_df = circuit_order_df.sort_values(by=['logical_qubits', 'gates'])
    ordered_circuits = circuit_order_df['circuit_clean'].tolist()

    circuit_map = {name: f"C{i+1}" for i, name in enumerate(ordered_circuits)}

    # 2. INICJALIZACJA PLIKU LATEX I LEGENY
    os.makedirs('wykresy', exist_ok=True)
    master_tex_path = os.path.join('wykresy', 'wszystkie_tabele_zagregowane.tex')
    with open(master_tex_path, 'w', encoding='utf-8') as f:
        f.write("% Wysoce zagregowany plik z tabelami wg standardów naukowych\n\n")

    write_legend_table(circuit_map, master_tex_path)

    layouts = df['layout_strategy'].unique()
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
    palette_algos = "Set1"

    print(f"🔍 Rozpoczynam generowanie wykresów i zagregowanych tabel...")

    for layout in layouts:
        layout_safe_name = sanitize_filename(layout)
        layout_dir = os.path.join('wykresy', layout_safe_name)
        os.makedirs(layout_dir, exist_ok=True)

        df_layout = df[df['layout_strategy'] == layout]
        topologies = df_layout['topology'].unique()

        for topo in topologies:
            df_topo = df_layout[df_layout['topology'] == topo].copy()

            # Wymuszamy kategoryzację obwodów by zachować globalną kolejność
            df_topo['circuit_clean'] = pd.Categorical(
                df_topo['circuit_clean'],
                categories=ordered_circuits,
                ordered=True
            )

            # ... TUTAJ PONIŻEJ MOŻESZ WSTAWIĆ SWOJE WYKRESY (ETAPY 1-4 Z POPRZEDNIEGO KODU) ...
            # Skrypt wykresów pozostawiam nienaruszony z wcześniejszego kodu - rysowanie wykresów się nie zmienia.

            # ---------------------------------------------------------
            # GENEROWANIE NOWYCH ZAGREGOWANYCH TABEL (ZAMIAST ETAPU 5)
            # ---------------------------------------------------------
            append_aggregated_tables(df_topo, master_tex_path, layout, topo, circuit_map)

            # Opcjonalnie: Zrzut szerokich "brudnych" danych do pliku CSV (jako załącznik do internetu)
            tables_dir = os.path.join(layout_dir, 'tabele_csv_surowe')
            os.makedirs(tables_dir, exist_ok=True)
            df_mean_csv = df_topo.groupby(['circuit_clean', 'router_short'], observed=True)[['swaps', 'moves', 'weighted_cost', 'execution_time_s']].mean().unstack()
            df_mean_csv.to_csv(os.path.join(tables_dir, f'{sanitize_filename(topo)}_pelne_dane.csv'))

        print(f"✅ Przetworzono mapowanie: {layout}")

    # generate_convergence_plot()
    # generate_amdahl_plot(df)
    # generate_layout_comparison_plot(df)
    # generate_scalability_plot(df)

    print(f"\n🎉 Gotowe! Agregacja zakończona. Wyniki znajdziesz w: {master_tex_path}")

if __name__ == "__main__":
    getDiagrams()