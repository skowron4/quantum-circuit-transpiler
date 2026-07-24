import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import re

def sanitize_filename(name):
    """Usuwa niedozwolone znaki z nazw plików/folderów"""
    return re.sub(r'[^\w\-]', '_', str(name))

def append_latex_table(df, filepath, caption, label):
    """Generuje ostatecznie i poprawnie wyrównaną tabelę longtable bez przesunięć tekstowych."""
    df_tex = df.copy()

    df_tex.columns = [str(col).replace('%', '\\%').replace('_', '\\_') for col in df_tex.columns]

    if df_tex.index.name:
        index_name = str(df_tex.index.name).replace('_', '\\_')
    else:
        index_name = "Obwod"

    for col in df_tex.select_dtypes(include=['object']):
        df_tex[col] = df_tex[col].astype(str).str.replace('_', '\\_')

    num_cols = len(df.columns) + 1

    # --- STABILNE SKALOWANIE W EMBRACH (em) ---
    if num_cols <= 4:
        data_col_width = "7.5em"
    elif num_cols <= 5:
        data_col_width = "6.5em"
    else:
        data_col_width = "4.8em"  # Dla 6 kolumn, idealnie zapobiega rozjeżdżaniu tekstu

    # POPRAWKA: Dodany brakujący cudzysłów po {data_col_width}}}
    col_format = "l " + " ".join([f">{{\\centering\\arraybackslash}}m{{{data_col_width}}}"] * (num_cols - 1))

    lines = []
    lines.append("\\begingroup")
    lines.append("\\footnotesize")
    lines.append("\\setlength{\\tabcolsep}{4pt}") # Bezpieczny, stały margines wewnętrzny

    lines.append(f"\\begin{{longtable}}{{{col_format}}}")
    lines.append(f"\\caption{{{caption}}}\\label{{{label}}} \\\\")
    lines.append("\\hline")

    # --- NAGŁÓWEK NA PIERWSZEJ STRONIE ---
    headers = [index_name] + list(df_tex.columns)
    lines.append(" & ".join(headers) + " \\\\")
    lines.append("\\hline")
    lines.append("\\endfirsthead")

    # --- NAGŁÓWEK NA KOLEJNYCH STRONACH ---
    lines.append(f"\\multicolumn{{{num_cols}}}{{c}}{{\\tablename\\ \\thetable{{}} -- ciąg dalszy z poprzedniej strony}} \\\\")
    lines.append("\\hline")
    lines.append(" & ".join(headers) + " \\\\")
    lines.append("\\hline")
    lines.append("\\endhead")

    lines.append("\\hline")
    lines.append(f"\\multicolumn{{{num_cols}}}{{r}}{{Ciąg dalszy na następnej stronie...}} \\\\")
    lines.append("\\endfoot")

    lines.append("\\hline")
    lines.append("\\endlastfoot")

    # --- DANE ---
    # POPRAWKA: Naprawione wcięcia w pętlach i instrukcjach warunkowych
    for index, row in df_tex.iterrows():
        idx_str = str(index).replace('_', '\\_')

        row_vals = []
        for val in row:
            if pd.isna(val):
                row_vals.append("-")
            elif isinstance(val, (int, float)):
                if isinstance(val, float) and np.isinf(val):
                    row_vals.append("$\\infty$" if val > 0 else "$-\\infty$")
                elif isinstance(val, float) and val.is_integer():
                    row_vals.append(f"{int(val)}")
                elif isinstance(val, float):
                    row_vals.append(f"{val:.2f}")
                else:
                    row_vals.append(str(val))
            else:
                row_vals.append(str(val))

        line_str = f"{idx_str} & " + " & ".join(row_vals) + " \\\\"
        lines.append(line_str)

    lines.append("\\end{longtable}")
    lines.append("\\endgroup")
    lines.append("\n% ---------------------------------------------------------\n\n")

    tabular_body = "\n".join(lines)

    with open(filepath, 'a', encoding='utf-8') as f:
        f.write(tabular_body)

def generate_convergence_plot():
    csv_path = 'results/convergence.csv'
    if not os.path.exists(csv_path):
        return

    df_conv = pd.read_csv(csv_path)
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df_conv, x='generation', y='best_score', hue='circuit_size', palette='tab10')
    plt.title('Krzywa Zbieżności (Convergence Plot) dla DE+SABRE')
    plt.xlabel('Numer Generacji')
    plt.ylabel('Wartość funkcji celu (Best Score)')
    plt.yscale('log')
    plt.legend(title="Rozmiar obwodu")
    plt.tight_layout()
    os.makedirs('wykresy/dodatkowe', exist_ok=True)
    plt.savefig('wykresy/dodatkowe/01_convergence.png', dpi=300)
    plt.close()

def generate_amdahl_plot(df_time):
    cores = [1, 2, 4, 8, 16, 32]
    P = 0.95
    mean_times = df_time[df_time['router_short'] == 'DE+SABRE']['execution_time_s'].mean()
    if pd.isna(mean_times): return
    times = [mean_times * ((1 - P) + P / c) for c in cores]

    plt.figure(figsize=(10, 6))
    sns.lineplot(x=cores, y=times, marker='o', linewidth=2, color='darkred')
    plt.title("Symulowany Wykres Przyspieszenia (Prawo Amdahla) dla DE+SABRE", pad=15)
    plt.xlabel("Liczba rdzeni procesora")
    plt.ylabel("Estymowany czas kompilacji [s]")
    plt.xticks(cores)
    plt.grid(True, which="both", ls="--")
    plt.tight_layout()
    os.makedirs('wykresy/dodatkowe', exist_ok=True)
    plt.savefig('wykresy/dodatkowe/02_amdahl_speedup.png', dpi=300)
    plt.close()

def generate_layout_comparison_plot(df):
    layout_rename = {
        'identity_trivial': 'Identity (Trywialne)',
        'random': 'Random (Losowe)',
        'compact_connectivity': 'Heurystyka (Degree + Time-decay)',
        'qiskit_sabre_initial_layout': 'Qiskit Sabre Layout'
    }
    df_filtered = df[
        (df['circuit_clean'].str.contains('QAOA', case=False, na=False)) &
        (df['gates'] >= 800) & (df['gates'] <= 1200) &
        (df['topology'].str.contains('3x3', na=False)) &
        (df['router_short'] == 'DE+SABRE')
        ].copy()

    if df_filtered.empty: return

    df_filtered['layout_name'] = df_filtered['layout_strategy'].map(layout_rename).fillna(df_filtered['layout_strategy'])
    plt.figure(figsize=(10, 6))
    sns.barplot(data=df_filtered, x='layout_name', y='weighted_cost', errorbar='sd', palette='viridis', capsize=0.1)
    plt.title("Wpływ początkowego mapowania na końcowy Ważony Koszt Trasowania\n(Obwód QAOA ~1000g, Topologia: Siatka 3x3, Router: DE+SABRE)", pad=15)
    plt.xlabel("Metoda Mapowania Początkowego (Initial Layout)")
    plt.ylabel("Ważony Koszt Trasowania")
    plt.xticks(rotation=15)
    plt.tight_layout()
    os.makedirs('wykresy/dodatkowe', exist_ok=True)
    plt.savefig('wykresy/dodatkowe/03_layout_comparison.png', dpi=300)
    plt.close()

def generate_scalability_plot(df):
    topo_order = ['Trojkat 3-Rezonatory', 'Siatka 2x2', 'Siatka 3x3', 'Siatka 4x4', 'Siatka 5x5']
    df_scale = df[df['topology'].isin(topo_order)].copy()
    if df_scale.empty: return
    df_scale['topology'] = pd.Categorical(df_scale['topology'], categories=topo_order, ordered=True)
    df_grouped = df_scale.groupby(['topology', 'router_short'], observed=False)[['weighted_cost', 'execution_time_s']].mean().reset_index()

    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df_grouped, x='topology', y='weighted_cost', hue='router_short', style='router_short', markers=True, dashes=False, linewidth=2, markersize=10, palette='Set1')
    plt.title("Skalowalność algorytmów: Ważony Koszt Trasowania vs Rozmiar QPU", pad=15)
    plt.xlabel("Rozmiar Topologii")
    plt.ylabel("Średni Ważony Koszt Trasowania")
    plt.grid(True, ls="--")
    plt.legend(title="Algorytm")
    plt.tight_layout()
    os.makedirs('wykresy/dodatkowe', exist_ok=True)
    plt.savefig('wykresy/dodatkowe/04_scalability_cost.png', dpi=300)
    plt.close()

    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df_grouped, x='topology', y='execution_time_s', hue='router_short', style='router_short', markers=True, dashes=False, linewidth=2, markersize=10, palette='Set1')
    plt.title("Skalowalność czasowa: Czas Kompilacji vs Rozmiar QPU", pad=15)
    plt.xlabel("Rozmiar Topologii")
    plt.ylabel("Średni Czas Kompilacji [s] (skala log)")
    plt.yscale('log')
    plt.grid(True, which="both", ls="--")
    plt.legend(title="Algorytm")
    plt.tight_layout()
    plt.savefig('wykresy/dodatkowe/05_scalability_time.png', dpi=300)
    plt.close()

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
    # FILTROWANIE OBWODÓW (zgodnie z nowymi wytycznymi zliczania bramek)
    # =================================================================
    # 1. Dla 25, 30, 35, 38, 45, 55, 65, 70 kubitów nie umieszczamy 50 bramek
    c1_qubits = [25, 30, 35, 38, 45, 55, 65, 70]
    m1 = df['logical_qubits'].isin(c1_qubits) & ((df['gates'] == 50) | df['circuit_clean'].str.contains('_50g|50g', regex=True))

    # 2. Dla 30, 35, 38, 45, 55, 65, 70 kubitów nie umieszczamy QAOA 200 bramek
    c2_qubits = [30, 35, 38, 45, 55, 65, 70]
    m2 = df['logical_qubits'].isin(c2_qubits) & df['circuit_clean'].str.contains('QAOA', case=False) & ((df['gates'] == 200) | df['circuit_clean'].str.contains('_200g|200g', regex=True))

    # 3. Dla 55, 65, 70 kubitów nie umieszczamy VQE 200 bramek i QAOA 500 bramek
    c3_qubits = [55, 65, 70]
    m3_vqe = df['logical_qubits'].isin(c3_qubits) & df['circuit_clean'].str.contains('VQE', case=False) & ((df['gates'] == 200) | df['circuit_clean'].str.contains('_200g|200g', regex=True))
    m3_qaoa = df['logical_qubits'].isin(c3_qubits) & df['circuit_clean'].str.contains('QAOA', case=False) & ((df['gates'] == 500) | df['circuit_clean'].str.contains('_500g|500g', regex=True))
    m3 = m3_vqe | m3_qaoa

    # 4. Dla 70 kubitów nie umieszczamy QAOA 2000 bramek
    c4_qubits = [70]
    m4 = df['logical_qubits'].isin(c4_qubits) & df['circuit_clean'].str.contains('QAOA', case=False) & ((df['gates'] == 2000) | df['circuit_clean'].str.contains('_2000g|2000g', regex=True))

    # Zastosowanie masek (odrzucamy te, które spełniają powyższe warunki)
    df = df[~(m1 | m2 | m3 | m4)].copy()
    # =================================================================

    df.sort_values(by=['logical_qubits', 'gates'], inplace=True)

    df['swaps'] = df['swaps'].fillna(0)
    df['moves'] = df['moves'].fillna(0)
    df['weighted_cost'] = (df['swaps'] * 3.14) + (df['moves'] * 1.0)

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
    palette_algos = "Set1"

    layouts = df['layout_strategy'].unique()
    print(f"🔍 Znaleziono {len(layouts)} strategii layoutu. Rozpoczynam generowanie...")

    os.makedirs('wykresy', exist_ok=True)
    master_tex_path = os.path.join('wykresy', 'wszystkie_tabele.tex')
    with open(master_tex_path, 'w', encoding='utf-8') as f:
        f.write("% Plik zbiorczy ze wszystkimi tabelami LaTeX\n")
        f.write("% Wygenerowano automatycznie\n\n")

    for layout in layouts:
        layout_safe_name = sanitize_filename(layout)
        layout_dir = os.path.join('wykresy', layout_safe_name)
        os.makedirs(layout_dir, exist_ok=True)

        df_layout = df[df['layout_strategy'] == layout]
        topologies = df_layout['topology'].unique()

        for topo in topologies:
            df_topo = df_layout[df_layout['topology'] == topo].copy()

            # =========================================================
            # DODANE: SORTOWANIE OD NAJMNIEJSZEJ LICZBY KUBITÓW (POPRAWKA V2)
            # =========================================================
            # Krok 1: Wyciągamy tylko potrzebne kolumny i usuwamy puste nazwy (NaN)
            circuit_order_df = df_topo[['circuit_clean', 'logical_qubits', 'gates']].copy()
            circuit_order_df = circuit_order_df.dropna(subset=['circuit_clean'])

            # Krok 2: Sortujemy rosnąco najpierw po kubitach, potem po bramkach
            circuit_order_df = circuit_order_df.sort_values(by=['logical_qubits', 'gates'])

            # Krok 3: Usuwamy wszystkie powtórzenia nazw obwodów, zostawiając tylko pierwsze
            # (czyli to posortowane najwyżej, z najmniejszą liczbą kubitów/bramek)
            circuit_order_df = circuit_order_df.drop_duplicates(subset=['circuit_clean'], keep='first')

            # Tworzymy ostateczną, absolutnie unikalną listę kategorii
            ordered_circuits = circuit_order_df['circuit_clean'].tolist()

            # Nadajemy kolumnie 'circuit_clean' typ kategoryczny
            df_topo['circuit_clean'] = pd.Categorical(
                df_topo['circuit_clean'],
                categories=ordered_circuits,
                ordered=True
            )
            # =========================================================
            # =========================================================

            num_circuits = len(df_topo['circuit_clean'].unique())
            plot_height = max(5.0, num_circuits * 0.4)
            topo_safe_name = sanitize_filename(topo)
            title_suffix = f"Topologia: {topo}\nLayout: {layout}"

            # ---------------------------------------------------------
            # ETAP 1: WYKRESY ABSOLUTNE (STANDARDOWE)
            # ---------------------------------------------------------
            metrics = [
                ("swaps", "Liczba bramek SWAP", "01_swaps", False),
                ("moves", "Liczba bramek MOVE", "02_moves", False),
                ("weighted_cost", "Ważony Koszt Trasowania", "03_weighted_cost", False),
                ("execution_time_s", "Czas wykonania [s]", "04_time", True)
            ]

            for col_name, xlabel, file_suffix, is_log in metrics:
                plt.figure(figsize=(10, plot_height))
                ax = sns.barplot(
                    data=df_topo, y="circuit_clean", x=col_name,
                    hue="router_short", palette=palette_algos, errorbar=None
                )
                plt.title(f"{xlabel}\n{title_suffix}", pad=15)
                plt.xlabel(xlabel + (" (skala log)" if is_log else ""))
                plt.ylabel("Obwód")

                if is_log:
                    ax.set_xscale('log')

                plt.legend(title="Algorytm", bbox_to_anchor=(1.05, 1), loc='upper left')
                plt.tight_layout()
                plt.savefig(os.path.join(layout_dir, f'{topo_safe_name}_{file_suffix}.png'), dpi=300)
                plt.close()

            # ---------------------------------------------------------
            # ETAP 2: WYKRESY RÓŻNICY PROCENTOWEJ (%)
            # ---------------------------------------------------------
            df_mean = df_topo.groupby(['circuit_clean', 'router_short'], observed=True)[['swaps', 'moves', 'weighted_cost', 'execution_time_s']].mean().unstack()

            pairs_to_compare = [
                ('DE+SABRE', 'Sabre'),
                ('DE+SABRE', 'Qiskit'),
                ('Sabre', 'Qiskit')
            ]

            metric_labels = {'swaps': 'SWAP', 'moves': 'MOVE', 'weighted_cost': 'Koszt Ważony', 'execution_time_s': 'Czas'}

            for metric in ['swaps', 'moves', 'weighted_cost', 'execution_time_s']:
                if metric not in df_mean.columns.levels[0]:
                    continue

                for r_A, r_B in pairs_to_compare:
                    if r_A not in df_mean[metric] or r_B not in df_mean[metric]:
                        continue

                    val_A = df_mean[metric][r_A]
                    val_B = df_mean[metric][r_B]

                    # Tu zamieniamy 0 w mianowniku na NaN, żeby wykres ominął nieskończoność i się nie rozjechał
                    diff_series = ((val_A - val_B) / val_B.replace(0, np.nan)) * 100
                    diff_series = diff_series.where(~((val_A == 0) & (val_B == 0)), 0)
                    diff_series = diff_series.dropna()

                    if diff_series.empty:
                        continue

                    df_diff = diff_series.reset_index(name='diff')
                    label_better = f"{r_A} Lepszy (<= 0%)"
                    label_worse = f"{r_A} Gorszy (> 0%)"
                    df_diff['Kategoria'] = np.where(df_diff['diff'] <= 0, label_better, label_worse)

                    plt.figure(figsize=(10, plot_height))
                    ax = sns.barplot(
                        data=df_diff, y='circuit_clean', x='diff',
                        hue='Kategoria', palette={label_better: '#2ca02c', label_worse: '#d62728'},
                        dodge=False
                    )

                    plt.axvline(0, color='black', linewidth=1.5, linestyle='--')

                    if metric == 'execution_time_s':
                        xmin, xmax = plt.xlim()
                        plt.xlim(min(-100, xmin), min(xmax, 2000))
                    else:
                        xmin, xmax = plt.xlim()
                        plt.xlim(min(-10, xmin), max(10, xmax))

                    plt.legend(title="Status", bbox_to_anchor=(1.05, 1), loc='upper left')
                    metric_upper = metric_labels[metric]
                    plt.title(f"% Różnica: {r_A} vs {r_B} ({metric_upper})\n{title_suffix}", pad=15)
                    plt.xlabel(f"Różnica procentowa [%]\n(wartość < 0 oznacza zielony pasek w lewo)")
                    plt.ylabel("Obwód")
                    plt.tight_layout()

                    diff_filename = f'{topo_safe_name}_05_diff_{metric}_{r_A}_vs_{r_B}.png'
                    plt.savefig(os.path.join(layout_dir, diff_filename), dpi=300)
                    plt.close()

            # ---------------------------------------------------------
            # ETAP 3 i 4 ... (bez zmian, pominięte dla czytelności kodu w odpowiedzi, zakładam że wkleisz swoją resztę)
            # ---------------------------------------------------------
            plt.figure(figsize=(12, plot_height))

            ax = sns.boxplot(
                data=df_topo, y="circuit_clean", x="weighted_cost",
                hue="router_short", palette=palette_algos, showfliers=False
            )

            sns.stripplot(
                data=df_topo, y="circuit_clean", x="weighted_cost",
                hue="router_short", palette=palette_algos,
                dodge=True, alpha=0.8, size=6, linewidth=1, ax=ax
            )

            handles, labels = ax.get_legend_handles_labels()
            n_algos = len(df_topo['router_short'].unique())
            ax.legend(handles[:n_algos], labels[:n_algos], title="Algorytm", bbox_to_anchor=(1.05, 1), loc='upper left')

            plt.title(f"Rozkład Wyników z naniesionymi próbami - Ważony Koszt (SWAP*3.14 + MOVE*1.0)\n{title_suffix}", pad=15)
            plt.xlabel("Koszt operacji trasowania (Niżej = Lepiej)")
            plt.ylabel("Obwód")
            plt.tight_layout()
            plt.savefig(os.path.join(layout_dir, f'{topo_safe_name}_06_boxplot_weighted_cost.png'), dpi=300)
            plt.close()

            # ---------------------------------------------------------
            # ETAP 4: FRONTY PARETO (W pętli dla SWAP, MOVE i Kosztu Ważonego)
            # ---------------------------------------------------------
            pareto_metrics = [
                ("swaps", "Liczba bramek SWAP", "07_pareto_swaps"),
                ("moves", "Liczba bramek MOVE", "08_pareto_moves"),
                ("weighted_cost", "Ważony Koszt Trasowania", "09_pareto_weighted_cost")
            ]

            for y_metric, y_label, file_suffix in pareto_metrics:
                # Omijamy rysowanie jeśli dane są zupełnie puste dla danej topologii
                if df_topo[y_metric].isna().all():
                    continue

                g = sns.relplot(
                    data=df_topo, x="execution_time_s", y=y_metric,
                    hue="router_short", style="router_short",
                    col="circuit_clean", col_wrap=3,
                    palette=palette_algos, s=120, alpha=0.8, kind="scatter",
                    facet_kws={'sharey': False, 'sharex': False},
                    height=4, aspect=1.2
                )

                g.set_axis_labels("Czas kompilacji [s] (log)", y_label)
                g.set(xscale="log")
                g.fig.subplots_adjust(top=0.85, hspace=0.4, wspace=0.3)
                g.fig.suptitle(f"Front Pareto: Czas vs {y_label}\n{title_suffix}")

                plt.savefig(os.path.join(layout_dir, f'{topo_safe_name}_{file_suffix}.png'), dpi=300)
                plt.close()

            # ---------------------------------------------------------
            # ETAP 5: GENEROWANIE TABEL
            # ---------------------------------------------------------
            tables_dir = os.path.join(layout_dir, 'tabele_csv')
            os.makedirs(tables_dir, exist_ok=True)

            metrics_info = {
                'swaps': 'Dodane bramki SWAP',
                'moves': 'Dodane bramki MOVE',
                'weighted_cost': 'Koszt trasowania',
                'execution_time_s': 'Czas kompilacji [s]'
            }

            # Bezpieczna funkcja do liczenia % na potrzeby tabel (obsługuje Nieskończoność i ułamki bliskie zeru)
            def calc_pct_diff(A, B):
                res = []
                for a, b in zip(A, B):
                    if pd.isna(a) or pd.isna(b):
                        res.append(np.nan)
                    # Ucinamy absurdy: jeśli baza jest mniejsza niż 0.01, uznajemy ją za zero
                    elif b < 0.01:
                        res.append(0.0 if a < 0.01 else np.inf)
                    else:
                        res.append(((a - b) / b) * 100)
                return pd.Series(res, index=A.index)

            for metric_col, metric_label in metrics_info.items():
                if metric_col not in df_topo.columns or df_topo[metric_col].isna().all():
                    continue

                df_metric = df_topo.groupby(['circuit_clean', 'router_short'], observed=True)[metric_col].mean().unstack()

                for r in ['Sabre', 'DE+SABRE', 'Qiskit']:
                    if r not in df_metric.columns:
                        df_metric[r] = np.nan

                tb = pd.DataFrame(index=df_metric.index)
                tb.index.name = "Obwod"

                # Używamy ZNACZNIE krótszych nazw - bez powtarzania "Dodane bramki..."
                tb['Sabre'] = df_metric['Sabre']
                tb['DE+SABRE'] = df_metric['DE+SABRE']

                has_qiskit = df_metric['Qiskit'].notna().any()

                if has_qiskit:
                    tb['Qiskit'] = df_metric['Qiskit']
                    # Zawijanie wyrazów (spacje) zamiast nawiasów dla lepszego łamania
                    tb['Różnica % DE+SABRE vs Qiskit'] = calc_pct_diff(df_metric['DE+SABRE'], df_metric['Qiskit'])
                    tb['Różnica % Sabre vs Qiskit'] = calc_pct_diff(df_metric['Sabre'], df_metric['Qiskit'])
                else:
                    tb['Różnica % DE+SABRE vs Sabre'] = calc_pct_diff(df_metric['DE+SABRE'], df_metric['Sabre'])

                tb_rounded = tb.round(2)
                tb_rounded.to_csv(os.path.join(tables_dir, f'{topo_safe_name}_{metric_col}.csv'))

                caption = f"Porównanie: {metric_label} dla topologii {topo} (Layout: {layout})."
                if metric_col == 'weighted_cost':
                    caption += " Koszt ważony (SWAP $\\times$ 3.14 + MOVE $\\times$ 1.0)."

                append_latex_table(
                    tb_rounded,
                    master_tex_path,
                    caption=caption.replace('_', '\\_'),
                    label=f"tab:{metric_col}_{topo_safe_name}_{layout_safe_name}"
                )

        print(f"✅ Przetworzono wykresy i nowe tabele dla układu: {layout}")

    generate_convergence_plot()
    generate_amdahl_plot(df)
    generate_layout_comparison_plot(df)
    generate_scalability_plot(df)

    print(f"\n🎉 Gotowe! Wszystkie tabele w jednym pliku LaTeX: {master_tex_path}")

if __name__ == "__main__":
    getDiagrams()