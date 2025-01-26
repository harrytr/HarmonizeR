import os
import pandas as pd
import networkx as nx
from sklearn.linear_model import LinearRegression
import numpy as np
import argparse

def load_grns(grn_dir):
    print("Loading GRNs...")
    grn_dict = {}
    for folder in os.listdir(grn_dir):
        folder_path = os.path.join(grn_dir, folder)
        if os.path.isdir(folder_path):
            graphml_files = [f for f in os.listdir(folder_path) if f.endswith(".graphml")]
            for graphml_file in graphml_files:
                grn_path = os.path.join(folder_path, graphml_file)
                grn_dict[graphml_file] = nx.read_graphml(grn_path)
    return grn_dict


def prune_with_anm(graph, expression_data):
    print("Pruning with Additive Noise Models (ANMs)...")
    genes = list(graph.nodes)
    expression_data = expression_data[genes]

    for u, v in list(graph.edges):
        try:
            X = expression_data[[u]].values
            Y = expression_data[[v]].values
            model = LinearRegression().fit(X, Y)
            residuals = Y - model.predict(X)

            corr = np.corrcoef(residuals.T, X.T)[0, 1]
            if abs(corr) > 0.05:
                graph.remove_edge(u, v)
        except:
            graph.remove_edge(u, v)

    return graph

def main(grn_dir, expression_path, output_dir):
    grn_dict = load_grns(grn_dir)
    print("Loading expression matrix...")
    expression_matrix = pd.read_csv(expression_path, index_col=0)

    all_genes = set(expression_matrix.columns)
    os.makedirs(output_dir, exist_ok=True)

    print("Processing GRNs...")
    for grn_name, grn in grn_dict.items():
        print(f"Processing {grn_name}...")

        grn = grn.subgraph([node for node in grn.nodes if node in all_genes])
        pruned_grn = prune_with_anm(grn, expression_matrix)

        #pruned_name = grn_name.replace(".graphml", "_pruned.graphml")
        output_path = os.path.join(output_dir, grn_name)
        nx.write_graphml(pruned_grn, output_path)

    print(f"Pruned GRNs saved to {output_dir}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prune GRNs using causal disentanglement methods.")
    parser.add_argument("grn_dir", type=str, help="Path to the directory containing folders with GRN GraphML files.")
    parser.add_argument("expression_path", type=str, help="Path to the expression matrix CSV file.")
    parser.add_argument("output_dir", type=str, help="Path to the directory to save pruned GRNs.")
    args = parser.parse_args()

    main(args.grn_dir, args.expression_path, args.output_dir)
