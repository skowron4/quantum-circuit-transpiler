import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import scipy.stats as stats

# Globalny słownik kolorów dla stałego przypisania barw do algorytmów
PALETTE_ALGOS = {
    'Sabre': '#e41a1c',     # Czerwony
    'DE+SABRE': '#377eb8',  # Niebieski
    'Qiskit': '#4daf4a'     # Zielony
}

def sanitize_filename(name):
    """Usuwa niedozwolone znaki z nazw plików/folderów"""
    return re.sub(r'[^\w\-]', '_', str(name))

def generate_convergence_plot():
    """Krzywe Zbieżności (Convergence Plots)."""
    csv_path = 'results/convergence.csv'
    if not os.path.exists(csv_path):
        print(f"⚠️ Pominięto Krzywe Zbieżności: Brak pliku {csv_path}.")
        return

    df_conv = pd.read_csv(csv_path)
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df_conv, x='generation', y='best_score', hue='circuit_size', palette='tab10')
    # Tytuły usunięte dla spójności z LaTeX
    plt.xlabel('Numer Generacji')
    plt.ylabel('Wartość funkcji celu (Best Score)')
    plt.yscale('log')
    plt.legend(title="Rozmiar obwodu")
    plt.tight_layout()
    os.makedirs('wykresy/dodatkowe', exist_ok=True)
    plt.savefig('wykresy/dodatkowe/01_convergence.png', dpi=300)
    plt.close()

def generate_amdahl_plot(df_time):
    """Symulowanego Wykres Przyspieszenia (Amdahl's Law)."""
    cores = [1, 2, 4, 8, 16, 32]
    P = 0.95

    mean_times = df_time[df_time['router_short'] == 'DE+SABRE']['execution_time_s'].mean()
    if pd.isna(mean_times): return

    times = [mean_times * ((1 - P) + P / c) for c in cores]

    plt.figure(figsize=(10, 6))
    sns.lineplot(x=cores, y=times, marker='o', linewidth=2, color='darkred')
    # Tytuły usunięte dla spójności z LaTeX
    plt.xlabel("Liczba rdzeni procesora")
    plt.ylabel("Estymowany czas kompilacji [s]")
    plt.xticks(cores)
    plt.grid(True, which="both", ls="--")
    plt.tight_layout()

    os.makedirs('wykresy/dodatkowe', exist_ok=True)
    plt.savefig('wykresy/dodatkowe/02_amdahl_speedup.png', dpi=300)
    plt.close()

def generate_layout_comparison_plot(df):
    """Wykres wpływu początkowego mapowania (Layout) na końcowy koszt."""
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

    if df_filtered.empty:
        print("⚠️ Pominięto wykres wpływu layoutu: Brak danych dla QAOA ~1000g na Siatce 3x3.")
        return

    df_filtered['layout_name'] = df_filtered['layout_strategy'].map(layout_rename).fillna(df_filtered['layout_strategy'])

    plt.figure(figsize=(10, 6))
    ax = sns.barplot(
        data=df_filtered, x='layout_name', y='weighted_cost',
        errorbar='sd', palette='viridis', capsize=0.1
    )
    # Tytuły usunięte dla spójności z LaTeX
    plt.xlabel("Metoda Mapowania Początkowego (Initial Layout)")
    plt.ylabel("Ważony Koszt Trasowania (SWAP*3.14 + MOVE*1.0)")
    plt.xticks(rotation=15)
    plt.tight_layout()

    os.makedirs('wykresy/dodatkowe', exist_ok=True)
    plt.savefig('wykresy/dodatkowe/03_layout_comparison.png', dpi=300)
    plt.close()

def generate_scalability_plot(df):
    """Wykres skalowalności ze wzrostem rozmiaru topologii sprzętowej."""
    topo_order = ['Trojkat 3-Rezonatory', 'Siatka 2x2', 'Siatka 3x3', 'Siatka 4x4', 'Siatka 5x5']

    df_scale = df[df['topology'].isin(topo_order)].copy()
    if df_scale.empty:
        return

    df_scale['topology'] = pd.Categorical(df_scale['topology'], categories=topo_order, ordered=True)
    df_grouped = df_scale.groupby(['topology', 'router_short'], observed=False)[['weighted_cost', 'execution_time_s']].mean().reset_index()

    plt.figure(figsize=(10, 6))
    sns.lineplot(
        data=df_grouped, x='topology', y='weighted_cost',
        hue='router_short', style='router_short',
        markers=True, dashes=False, linewidth=2, markersize=10, palette=PALETTE_ALGOS
    )
    # Tytuły usunięte dla spójności z LaTeX
    plt.xlabel("Rozmiar Topologii")
    plt.ylabel("Średni Ważony Koszt Trasowania")
    plt.grid(True, ls="--")
    plt.legend(title="Algorytm")
    plt.tight_layout()
    os.makedirs('wykresy/dodatkowe', exist_ok=True)
    plt.savefig('wykresy/dodatkowe/04_scalability_cost.png', dpi=300)
    plt.close()

    plt.figure(figsize=(10, 6))
    sns.lineplot(
        data=df_grouped, x='topology', y='execution_time_s',
        hue='router_short', style='router_short',
        markers=True, dashes=False, linewidth=2, markersize=10, palette=PALETTE_ALGOS
    )
    # Tytuły usunięte dla spójności z LaTeX
    plt.xlabel("Rozmiar Topologii")
    plt.ylabel("Średni Czas Kompilacji [s] (skala log)")
    plt.yscale('log')
    plt.grid(True, which="both", ls="--")
    plt.legend(title="Algorytm")
    plt.tight_layout()
    plt.savefig('wykresy/dodatkowe/05_scalability_time.png', dpi=300)
    plt.close()

def generate_latex_plots(df):
    """
    Generuje wykresy przeznaczone bezpośrednio do artykułu w LaTeX.
    Zapisuje je w folderze 'latex/nazwa_layoutu/'.
    """
    base_latex_dir = 'latex'
    os.makedirs(base_latex_dir, exist_ok=True)

    target_filters = {
        'Trojkat 3-Rezonatory': [10, 13],
        'Siatka 5x5': [55, 70]
    }

    print("📊 Rozpoczynam generowanie dedykowanych wykresów do folderu 'latex'...")

    for topo, qubits_list in target_filters.items():
        df_topo = df[df['topology'] == topo].copy()
        if df_topo.empty:
            continue

        for q_val in qubits_list:
            df_filtered = df_topo[df_topo['logical_qubits'] == q_val].copy()

            # --- FILTROWANIE PO NAZWACH OBWODÓW (REGEX) ---
            if topo == 'Siatka 5x5' and q_val == 55:
                mask_50 = df_filtered['circuit_clean'].str.contains(r'(?<!\d)50(?!\d)', regex=True)
                mask_200 = df_filtered['circuit_clean'].str.contains(r'(?<!\d)200(?!\d)', regex=True)
                mask_qaoa_1000 = df_filtered['circuit_clean'].str.contains(r'QAOA.*(?<!\d)1000(?!\d)', case=False, regex=True)

                df_filtered = df_filtered[~(mask_50 | mask_200 | mask_qaoa_1000)]

            if topo == 'Siatka 5x5' and q_val == 70:
                mask_200 = df_filtered['circuit_clean'].str.contains(r'(?<!\d)200(?!\d)', regex=True)
                mask_qaoa_500 = df_filtered['circuit_clean'].str.contains(r'QAOA.*(?<!\d)500(?!\d)', case=False, regex=True)
                mask_qaoa_2000 = df_filtered['circuit_clean'].str.contains(r'QAOA.*(?<!\d)2000(?!\d)', case=False, regex=True)

                df_filtered = df_filtered[~(mask_200 | mask_qaoa_500 | mask_qaoa_2000)]
            # --------------------------------------------------------

            if df_filtered.empty:
                print(f"⚠️ Brak danych (po filtracji) dla topologii '{topo}' z liczbą kubitów: {q_val}")
                continue

            layouts = df_filtered['layout_strategy'].unique()
            for layout in layouts:
                df_layout = df_filtered[df_filtered['layout_strategy'] == layout].copy()

                if df_layout.empty:
                    continue

                layout_safe = sanitize_filename(layout)
                topo_safe = sanitize_filename(topo)

                # Tworzymy podfolder dla danego mapowania (Initial Layout)
                layout_dir = os.path.join(base_latex_dir, layout_safe)
                os.makedirs(layout_dir, exist_ok=True)

                num_circuits = len(df_layout['circuit_clean'].unique())
                plot_height = max(4.0, num_circuits * 0.45)

                metrics = [
                    ("swaps", "Liczba bramek SWAP", "01_swaps", False),
                    ("moves", "Liczba bramek MOVE", "02_moves", False),
                    ("weighted_cost", "Ważony Koszt Trasowania", "03_weighted_cost", False),
                    ("execution_time_s", "Czas wykonania [s]", "04_time", True)
                ]

                for col_name, xlabel, file_suffix, is_log in metrics:
                    plt.figure(figsize=(10, plot_height))
                    ax = sns.barplot(
                        data=df_layout, y="circuit_clean", x=col_name,
                        hue="router_short", palette=PALETTE_ALGOS, errorbar=None
                    )

                    # BRAK TYTUŁÓW!
                    plt.xlabel(xlabel + (" (skala log)" if is_log else ""))
                    plt.ylabel("Obwód")

                    if is_log:
                        ax.set_xscale('log')

                    plt.legend(title="Algorytm", bbox_to_anchor=(1.05, 1), loc='upper left')
                    plt.tight_layout()

                    # Zapisujemy do folderu layoutu, nazwa jest prostsza
                    filename = f"{topo_safe}_Q{q_val}_{file_suffix}.png"
                    plt.savefig(os.path.join(layout_dir, filename), dpi=300)
                    plt.close()

                # WYKRES PUDEŁKOWY
                plt.figure(figsize=(11, plot_height))
                ax = sns.boxplot(
                    data=df_layout, y="circuit_clean", x="weighted_cost",
                    hue="router_short", palette=PALETTE_ALGOS, showfliers=False
                )
                sns.stripplot(
                    data=df_layout, y="circuit_clean", x="weighted_cost",
                    hue="router_short", palette=PALETTE_ALGOS,
                    dodge=True, alpha=0.8, size=5, linewidth=1, ax=ax
                )
                handles, labels = ax.get_legend_handles_labels()
                unique_algos = df_layout['router_short'].dropna().unique()
                n_algos = len(unique_algos)
                ax.legend(handles[:n_algos], labels[:n_algos], title="Algorytm", bbox_to_anchor=(1.05, 1), loc='upper left')

                # BRAK TYTUŁÓW!
                plt.xlabel("Koszt operacji trasowania (Niżej = Lepiej)")
                plt.ylabel("Obwód")
                plt.tight_layout()

                box_filename = f"{topo_safe}_Q{q_val}_06_boxplot.png"
                plt.savefig(os.path.join(layout_dir, box_filename), dpi=300)
                plt.close()

    print("✅ Zakończono generowanie wykresów do folderu 'latex/'.")

def getDiagrams():
    csv_path = 'results/routing_runs.csv'
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"❌ Nie znaleziono pliku: {csv_path}")
        return

    cols_to_numeric = ['gates', 'swaps', 'moves', 'execution_time_s', 'makespan']
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
    df.sort_values(by=['logical_qubits', 'gates'], inplace=True)

    df['swaps'] = df['swaps'].fillna(0)
    df['moves'] = df['moves'].fillna(0)
    df['weighted_cost'] = (df['swaps'] * 3.14) + (df['moves'] * 1.0)

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)

    layouts = df['layout_strategy'].unique()
    print(f"🔍 Znaleziono {len(layouts)} strategii layoutu. Rozpoczynam generowanie standardowe...")

    for layout in layouts:
        layout_safe_name = sanitize_filename(layout)
        layout_dir = os.path.join('wykresy', layout_safe_name)
        os.makedirs(layout_dir, exist_ok=True)

        df_layout = df[df['layout_strategy'] == layout]
        topologies = df_layout['topology'].unique()

        for topo in topologies:
            df_topo = df_layout[df_layout['topology'] == topo].copy()
            num_circuits = len(df_topo['circuit_clean'].unique())
            plot_height = max(5.0, num_circuits * 0.4)
            topo_safe_name = sanitize_filename(topo)

            # ---------------------------------------------------------
            # ETAP 1: WYKRESY ABSOLUTNE
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
                    hue="router_short", palette=PALETTE_ALGOS, errorbar=None
                )
                plt.xlabel(xlabel + (" (skala log)" if is_log else ""))
                plt.ylabel("Obwód")

                if is_log:
                    ax.set_xscale('log')

                plt.legend(title="Algorytm", bbox_to_anchor=(1.05, 1), loc='upper left')
                plt.tight_layout()
                plt.savefig(os.path.join(layout_dir, f'{topo_safe_name}_{file_suffix}.png'), dpi=300)
                plt.close()

            # ---------------------------------------------------------
            # ETAP 2: RÓŻNICE PROCENTOWE
            # ---------------------------------------------------------
            df_mean = df_topo.groupby(['circuit_clean', 'router_short'])[['swaps', 'moves', 'weighted_cost', 'execution_time_s']].mean().unstack()
            pairs_to_compare = [('DE+SABRE', 'Sabre'), ('DE+SABRE', 'Qiskit'), ('Sabre', 'Qiskit')]
            metric_labels = {'swaps': 'SWAP', 'moves': 'MOVE', 'weighted_cost': 'Koszt Ważony', 'execution_time_s': 'Czas'}

            for metric in ['swaps', 'moves', 'weighted_cost', 'execution_time_s']:
                if metric not in df_mean.columns.levels[0]: continue

                for r_A, r_B in pairs_to_compare:
                    if r_A not in df_mean[metric] or r_B not in df_mean[metric]: continue

                    val_A = df_mean[metric][r_A]
                    val_B = df_mean[metric][r_B]
                    diff_series = ((val_A - val_B) / val_B.replace(0, 1e-9)) * 100
                    diff_series = diff_series.where(~((val_A == 0) & (val_B == 0)), 0).dropna()

                    if diff_series.empty: continue

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
                    plt.xlabel("Różnica procentowa [%]")
                    plt.ylabel("Obwód")
                    plt.tight_layout()
                    plt.savefig(os.path.join(layout_dir, f'{topo_safe_name}_05_diff_{metric}_{r_A}_vs_{r_B}.png'), dpi=300)
                    plt.close()

            # ---------------------------------------------------------
            # ETAP 3: WYKRESY PUDEŁKOWE
            # ---------------------------------------------------------
            plt.figure(figsize=(12, plot_height))
            ax = sns.boxplot(
                data=df_topo, y="circuit_clean", x="weighted_cost",
                hue="router_short", palette=PALETTE_ALGOS, showfliers=False
            )
            sns.stripplot(
                data=df_topo, y="circuit_clean", x="weighted_cost",
                hue="router_short", palette=PALETTE_ALGOS,
                dodge=True, alpha=0.8, size=6, linewidth=1, ax=ax
            )
            handles, labels = ax.get_legend_handles_labels()
            unique_algos = df_topo['router_short'].dropna().unique()
            n_algos = len(unique_algos)
            ax.legend(handles[:n_algos], labels[:n_algos], title="Algorytm", bbox_to_anchor=(1.05, 1), loc='upper left')

            plt.xlabel("Koszt operacji trasowania (Niżej = Lepiej)")
            plt.ylabel("Obwód")
            plt.tight_layout()
            plt.savefig(os.path.join(layout_dir, f'{topo_safe_name}_06_boxplot_weighted_cost.png'), dpi=300)
            plt.close()

            # ---------------------------------------------------------
            # ETAP 4: FRONTY PARETO
            # ---------------------------------------------------------
            pareto_metrics = [
                ("swaps", "Liczba bramek SWAP", "07_pareto_swaps"),
                ("moves", "Liczba bramek MOVE", "08_pareto_moves"),
                ("weighted_cost", "Ważony Koszt Trasowania", "09_pareto_weighted_cost")
            ]

            for y_metric, y_label, file_suffix in pareto_metrics:
                if df_topo[y_metric].isna().all(): continue

                g = sns.relplot(
                    data=df_topo, x="execution_time_s", y=y_metric,
                    hue="router_short", style="router_short",
                    col="circuit_clean", col_wrap=3,
                    palette=PALETTE_ALGOS, s=120, alpha=0.8, kind="scatter",
                    facet_kws={'sharey': False, 'sharex': False},
                    height=4, aspect=1.2
                )
                g.set_axis_labels("Czas kompilacji [s] (log)", y_label)
                g.set(xscale="log")
                g.fig.subplots_adjust(top=0.85, hspace=0.4, wspace=0.3)
                # Tytuły Pareto też usunięte z poziomu obrazka
                g.fig.suptitle("")
                plt.savefig(os.path.join(layout_dir, f'{topo_safe_name}_{file_suffix}.png'), dpi=300)
                plt.close()

            # --- TABELE ---
            tables_dir = os.path.join(layout_dir, 'tabele')
            os.makedirs(tables_dir, exist_ok=True)
            df_means = df_topo.groupby(['circuit_clean', 'router_short'])['swaps'].mean().unstack().fillna(0)
            df_base_gates = df_topo.groupby('circuit_clean')['gates'].first()
            tb1 = pd.DataFrame(index=df_means.index)
            tb1['Oryginalna liczba bramek'] = df_base_gates
            for r in ['Qiskit', 'Sabre', 'DE+SABRE']:
                if r in df_means.columns: tb1[f'Bramki po {r}'] = df_base_gates + (df_means[r] * 3)
            if 'DE+SABRE' in df_means.columns:
                tb1['Narzut DE+SABRE (%)'] = ((tb1['Bramki po DE+SABRE'] - tb1['Oryginalna liczba bramek']) / tb1['Oryginalna liczba bramek']) * 100
            tb1.round(2).to_csv(os.path.join(tables_dir, f'{topo_safe_name}_Tabela_1_Narzut.csv'))

    generate_convergence_plot()
    generate_amdahl_plot(df)
    generate_layout_comparison_plot(df)
    generate_scalability_plot(df)

    # Wywołanie nowej funkcji LaTeXowej
    generate_latex_plots(df)

    print("\n🎉 Gotowe! Wszystkie zaktualizowane wykresy znajdziesz w podfolderach 'wykresy/' oraz 'latex/'.")

if __name__ == "__main__":
    getDiagrams()