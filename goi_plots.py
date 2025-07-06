import pandas as pd
import numpy as np
from plotnine import *
from scipy.stats import f_oneway, kruskal

data = "CCLE"

# Load and rename
df = pd.read_csv("multi_scenario_results_" + data + ".csv")
df.rename(columns={
    "GOI": "Gene",
    "hop_order": "Hops",
    "misclassification_percent": "Misclassification"
}, inplace=True)

# Function to compute ANOVA and Kruskal–Wallis per hop
def get_stats_for_hop(hop_value):
    subset = df[df["Hops"] == hop_value]
    grouped = [group["Misclassification"].values for _, group in subset.groupby("Gene")]
    
    if len(grouped) < 2:
        return f"Hop {hop_value}: Insufficient groups"
    
    anova_p = f_oneway(*grouped).pvalue
    kruskal_p = kruskal(*grouped).pvalue
    
    return f"Hop {hop_value}: ANOVA p = {anova_p:.5f}, Kruskal–Wallis p = {kruskal_p:.5f}"

# Compute stats annotation
hop_values = sorted(df["Hops"].unique())
stat_labels = "\n".join(get_stats_for_hop(h) for h in hop_values)
y_max = df["Misclassification"].max()

# Compute means to display
mean_labels_df = (
    df.groupby(["Gene", "Hops"])["Misclassification"]
    .mean()
    .reset_index()
)
mean_labels_df["label"] = mean_labels_df.apply(
    lambda row: f"Mean (Hop = {int(row['Hops'])}): {row['Misclassification']:.1f}", axis=1
)

y_min = df["Misclassification"].min()
# Stagger y-positions based on hop value
mean_labels_df["y"] = mean_labels_df["Hops"].apply(lambda h: y_min - 1.5 if h == 2 else y_min - 0.5)


# Plot
plot = (
    ggplot(df, aes(x='Gene', y='Misclassification', fill='factor(Hops)'))
    + geom_violin(width=0.9, alpha=0.6, color="gray", draw_quantiles=[0.25, 0.5, 0.75])
    + geom_boxplot(width=0.15, outlier_shape='', color='black', fill='white')
    + geom_text(
        aes(x=1.5, y=y_max - 2),
        label=stat_labels,
        size=12,
        ha='center',
        va='bottom',
        inherit_aes=False
    )
    + geom_text(
        data=mean_labels_df,
        mapping=aes(x='Gene', y='y', label='label', group='Hops'),
        size=12,
        color="black",
        va='bottom',
        inherit_aes=False
    )

    + theme_bw()
    + theme(
        figure_size=(12, 7),
        axis_text_x=element_text(rotation=0, size=11),
        axis_text_y=element_text(size=11),
        axis_title=element_text(size=13),
        legend_title=element_text(size=12),
        legend_text=element_text(size=11)
    )
    + labs(
        title="Misclassification Across GOI genes and Hop Lengths in " + data,
        x="Gene(s) of Interest",
        y="Misclassification %",
        fill="Hop Length"
    )
)

plot.save("S5_" + data + ".png", dpi=600)
