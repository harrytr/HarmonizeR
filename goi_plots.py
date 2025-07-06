import pandas as pd
import numpy as np
from plotnine import *
from scipy.stats import ttest_ind, mannwhitneyu

# Load and rename
df = pd.read_csv("multi_scenario_results_TCGA.csv")
df.rename(columns={
    "GOI": "Gene",
    "hop_order": "Hops",
    "misclassification_percent": "Misclassification"
}, inplace=True)

# Extract hop groups
group1 = df[df["Hops"] == 1]["Misclassification"]
group2 = df[df["Hops"] == 2]["Misclassification"]

# Welch's t-test
t_stat, t_p = ttest_ind(group1, group2, equal_var=False)

# Mann–Whitney U test
u_stat, u_p = mannwhitneyu(group1, group2, alternative='two-sided')

# Create precise label
stat_label = (
    f"Welch's t-test p = {t_p:.5f}\n"
    f"Mann–Whitney U p = {u_p:.5f}"
)

# Compute max y for label placement
y_max = df["Misclassification"].max()

# Build plot
plot = (
    ggplot(df, aes(x='Gene', y='Misclassification', fill='factor(Hops)'))
    + geom_violin(width=0.9, alpha=0.6, color="gray", draw_quantiles=[0.25, 0.5, 0.75])
    + geom_boxplot(width=0.15, outlier_shape='', color='black', fill='white')
    + geom_text(
        aes(x=1.5, y=y_max + 2),
        label=stat_label,
        size=9,
        ha='center',
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
        title="Misclassification Across GOI genes and Hop Lengths in TCGA",
        x="Gene(s) of Interest",
        y="Misclassification %",
        fill="Hop Length"
    )
)

# Save to file
plot.save("misclassification_global_test.png", dpi=600)
print("✅ Saved to Desktop with p-values to 5 decimal places")
