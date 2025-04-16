import os
import sys
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

# === HANDLE INPUT ARGUMENT ===
if len(sys.argv) != 2:
    print("Usage: python generate_combo_plots.py <SOURCE>")
    print("Example: python generate_combo_plots.py TCGA")
    sys.exit(1)

SOURCE = sys.argv[1]
CSV_PATH = f"GNN_labels_{SOURCE}.csv"
NETWORKS_DIR = f"./networks_{SOURCE}"
OUTPUT_FILE = f"mutation_combo_plots_{SOURCE}.png"
DPI = 600

# === LOAD CSV WITH LABELS ===
if not os.path.exists(CSV_PATH):
    print(f"CSV file not found: {CSV_PATH}")
    sys.exit(1)

if not os.path.exists(NETWORKS_DIR):
    print(f"Network folder not found: {NETWORKS_DIR}")
    sys.exit(1)

labels_df = pd.read_csv(CSV_PATH)
filename_to_label = dict(zip(labels_df["filename"], labels_df["labels"]))


def format_p(p):
    try:
        if p <= 0:
            return "< 1e-300"
        return f"{p:.2e}"
    except:
        return "N/A"


# === EXTRACT FEATURES FROM NETWORKS ===
def extract_features(graph):
    degrees = dict(graph.degree())
    clustering = nx.clustering(graph)

    try:
        avg_path_len = nx.average_shortest_path_length(graph)
    except Exception:
        avg_path_len = float('nan')  # handle disconnected graphs

    return {
        "num_nodes": graph.number_of_nodes(),
        "num_edges": graph.number_of_edges(),
        "density": nx.density(graph),
        "avg_clustering": nx.average_clustering(graph),
        "max_clustering": max(clustering.values()) if clustering else 0,
        "degree_variance": pd.Series(degrees.values()).var(),
        "num_triangles": sum(nx.triangles(graph).values()) // 3,
        "avg_path_length": avg_path_len
    }


data = []
for file in os.listdir(NETWORKS_DIR):
    if file.endswith(".graphml") and file in filename_to_label:
        path = os.path.join(NETWORKS_DIR, file)
        try:
            G = nx.read_graphml(path)
            G = nx.Graph(G)
            features = extract_features(G)
            features["filename"] = file
            features["label"] = filename_to_label[file]
            data.append(features)
        except Exception as e:
            print(f"Error processing {file}: {e}")

df = pd.DataFrame(data)

# === FILTER TO TOP 7 MUTATION TYPES ===
top_labels = df["label"].value_counts().head(7).index.tolist()
df_top = df[df["label"].isin(top_labels)].copy()

# === COMBO PLOTS: VIOLIN + BOX + SWARM (2 per row) ===
features_to_plot = [
    "avg_clustering", "degree_variance", "max_clustering", 
    "num_edges", "num_nodes", "density", "num_triangles",
    "avg_path_length" 
]

sns.set(style="whitegrid", palette="Set2", font_scale=1.1)
n_feats = len(features_to_plot)
ncols = 2
nrows = (n_feats + 1) // ncols

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(16, 5 * nrows))
axes = axes.flatten()

import pingouin as pg
from scipy.stats import kruskal

for i, feature in enumerate(features_to_plot):
    ax = axes[i]
    
    # Group values for this feature
    groups = [df_top[df_top["label"] == lbl][feature].dropna() for lbl in top_labels]



    # Welch's ANOVA using pingouin
    try:
        welch_df = pg.welch_anova(dv=feature, between='label', data=df_top)
        anova_p = welch_df['p-unc'].iloc[0]
    except Exception:
        anova_p = float("nan")

    # Kruskal-Wallis test
    try:
        groups = [df_top[df_top["label"] == lbl][feature].dropna() for lbl in top_labels]
        kruskal_stat, kruskal_p = kruskal(*groups)
    except Exception:
        kruskal_p = float("nan")



    # Combo plot: violin + box + strip
    sns.violinplot(data=df_top, x="label", y=feature, ax=ax, inner=None, linewidth=1.2, cut=0)
    sns.boxplot(data=df_top, x="label", y=feature, ax=ax, whis=1.5, width=0.2, showcaps=True, boxprops={'zorder': 2}, showfliers=False)
    sns.stripplot(data=df_top, x="label", y=feature, ax=ax, color='k', alpha=0.3, size=1.5, jitter=0.2, dodge=True)

    # Titles and formatting
    ax.set_title(
    f"{feature.replace('_', ' ').title()} by Mutation Type\n"
    f"Welch's ANOVA p={format_p(anova_p)} | Kruskal p={format_p(kruskal_p)}"
    )
    ax.set_xlabel("")
    ax.set_ylabel(feature.replace('_', ' ').title())
    ax.tick_params(axis='x', rotation=45)

# Remove any unused axes
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.savefig(OUTPUT_FILE, dpi=DPI)
print(f"Combo plots saved to {OUTPUT_FILE}")
