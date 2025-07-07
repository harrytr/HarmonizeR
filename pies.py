import pandas as pd
import matplotlib.pyplot as plt

def plot_mutation_distribution(csv_path, output_file):
    # Read the CSV
    df = pd.read_csv(csv_path)

    # Explicitly use 'labels' column if it exists
    if 'labels' in df.columns:
        label_series = df['labels']
    else:
        label_series = df.iloc[:, 1]  # fallback if no header

    # Count occurrences
    label_counts = label_series.value_counts()
    total = label_counts.sum()

    # Labels with percentages
    pie_labels = [f"{label} ({count / total * 100:.1f}%)" for label, count in label_counts.items()]

    # Plot
    fig, ax = plt.subplots(figsize=(10, 8))
    wedges, _ = ax.pie(label_counts, startangle=140)
    ax.set_title("Mutation Type Distribution (CCLE)", fontsize=25)
    ax.axis('equal')

    # External legend
    ax.legend(wedges, pie_labels, title="Mutation Types", loc="center left",
              bbox_to_anchor=(0.95, 0.5), fontsize=15, title_fontsize=16)

    # Reserve space for legend
    plt.subplots_adjust(right=0.65)

    # Save
    plt.savefig(output_file, dpi=600)
    plt.show()

if __name__ == "__main__":
    #plot_mutation_distribution("GNN_labels_TCGA.csv", "mutation_distribution_TCGA.png")
    plot_mutation_distribution("GNN_labels_CCLE.csv", "mutation_distribution_CCLE.png")
