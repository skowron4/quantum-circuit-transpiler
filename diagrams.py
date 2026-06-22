import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import re

def sanitize_filename(name):
    """Usuwa niedozwolone znaki z nazw plików/folderów"""
    return re.sub(r'[^\w\-]', '_', str(name))

def getDiagrams():
    # 1. Wczytanie danych
    csv_path = 'results/routing_runs.csv'
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"❌ Nie znaleziono pliku: {csv_path}")
        return

    # --- POPRAWKA: Wymuszenie konwersji wszystkich kolumn liczbowych na float/int ---
    cols_to_numeric = ['gates', 'swaps', 'moves', 'execution_time_s']
    for col in cols_to_numeric:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # Odfiltrowanie gigantycznych obwodów
    df = df[df['gates'] <= 3500].copy()

    # 2. Przygotowanie czytelnych etykiet algorytmów
    router_mapping = {
        'sabre_fixed_layout': 'Sabre',
        'dehybrid_sabre_fixed_layout': 'DEHybrid',
        'qiskit_sabre_export_eval': 'Qiskit'
    }
    df['router_short'] = df['router'].replace(router_mapping)

    # Przerzucenie Qiskita do 'qiskit_sabre_initial_layout', żeby był na jednym wykresie z innymi
    qiskit_mask = df['router'] == 'qiskit_sabre_export_eval'
    df.loc[qiskit_mask, 'layout_strategy'] = 'qiskit_sabre_initial_layout'

    # Oczyszczenie nazw obwodów
    df['circuit_clean'] = df['circuit_id'].str.replace('Constellation_', '').str.replace('_Proxy', '')

    # Sortowanie obwodów od najmniejszych do największych
    df.sort_values(by=['logical_qubits', 'gates'], inplace=True)

    # Ustawienia stylów Seaborn
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
    palette_algos = "Set1"

    # 3. Pętla główna: lecimy po każdym LAYOUCIE
    layouts = df['layout_strategy'].unique()
    print(f"🔍 Znaleziono {len(layouts)} strategii layoutu. Rozpoczynam generowanie...")

    for layout in layouts:
        layout_safe_name = sanitize_filename(layout)

        layout_dir = os.path.join('wykresy', layout_safe_name)
        os.makedirs(layout_dir, exist_ok=True)

        df_layout = df[df['layout_strategy'] == layout]
        topologies = df_layout['topology'].unique()

        for topo in topologies:
            df_topo = df_layout[df_layout['topology'] == topo]

            # Dynamiczna wysokość
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
                ("execution_time_s", "Czas wykonania [s]", "03_time", True)
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
            # ETAP 2: WYKRESY RÓŻNICY PROCENTOWEJ (%) - POPRAWIONE
            # ---------------------------------------------------------
            df_mean = df_topo.groupby(['circuit_clean', 'router_short'])[['swaps', 'moves', 'execution_time_s']].mean().unstack()

            pairs_to_compare = [
                ('DEHybrid', 'Sabre'),
                ('DEHybrid', 'Qiskit'),
                ('Sabre', 'Qiskit')
            ]

            metric_labels = {'swaps': 'SWAP', 'moves': 'MOVE', 'execution_time_s': 'Czas'}

            for metric in ['swaps', 'moves', 'execution_time_s']:
                if metric not in df_mean.columns.levels[0]:
                    continue

                for r_A, r_B in pairs_to_compare:
                    if r_A not in df_mean[metric] or r_B not in df_mean[metric]:
                        continue

                    val_A = df_mean[metric][r_A]
                    val_B = df_mean[metric][r_B]

                    # Obliczenie % różnicy
                    diff_series = ((val_A - val_B) / val_B.replace(0, 1e-9)) * 100
                    diff_series = diff_series.where(~((val_A == 0) & (val_B == 0)), 0)
                    diff_series = diff_series.dropna()

                    if diff_series.empty:
                        continue

                    df_diff = diff_series.reset_index(name='diff')

                    # POPRAWKA: Jawne kategorie tekstowe zamiast True/False
                    label_better = f"{r_A} Lepszy (<= 0%)"
                    label_worse = f"{r_A} Gorszy (> 0%)"

                    df_diff['Kategoria'] = np.where(df_diff['diff'] <= 0, label_better, label_worse)

                    plt.figure(figsize=(10, plot_height))

                    # POPRAWKA: dodge=False - wyłącza przesuwanie/zwężanie słupków
                    # POPRAWKA: twarde przypisanie stringów do palety kolorów
                    ax = sns.barplot(
                        data=df_diff, y='circuit_clean', x='diff',
                        hue='Kategoria',
                        palette={label_better: '#2ca02c', label_worse: '#d62728'},
                        dodge=False
                    )

                    plt.axvline(0, color='black', linewidth=1.5, linestyle='--')

                    # Zabezpieczenie przed gigantycznymi różnicami czasu
                    if metric == 'execution_time_s':
                        xmin, xmax = plt.xlim()
                        # Pozwala wykresowi pokazać do -100% z lewej i uciąć gigantyczne +% z prawej
                        plt.xlim(min(-100, xmin), min(xmax, 2000))
                    else:
                        # Dodajemy lekki margines dla czytelności bramek
                        xmin, xmax = plt.xlim()
                        plt.xlim(min(-10, xmin), max(10, xmax))

                    # Legenda ląduje ładnie z prawej strony
                    plt.legend(title="Status", bbox_to_anchor=(1.05, 1), loc='upper left')

                    metric_upper = metric_labels[metric]
                    plt.title(f"% Różnica: {r_A} vs {r_B} ({metric_upper})\n{title_suffix}", pad=15)
                    plt.xlabel(f"Różnica procentowa [%]\n(wartość < 0 oznacza zielony pasek w lewo)")
                    plt.ylabel("Obwód")
                    plt.tight_layout()

                    diff_filename = f'{topo_safe_name}_04_diff_{metric}_{r_A}_vs_{r_B}.png'
                    plt.savefig(os.path.join(layout_dir, diff_filename), dpi=300)
                    plt.close()

        print(f"✅ Przetworzono layout: {layout}")

    print("\n🎉 Gotowe! Wykresy (również te ze słupkami na plus!) znajdziesz w folderze 'wykresy/'.")

if __name__ == "__main__":
    getDiagrams()