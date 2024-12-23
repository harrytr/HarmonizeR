import os, sys, time
import pandas as pd
import igraph as ig
import torch
from torch_geometric.loader import DataLoader
from torch_geometric.data import Data
from torch.nn import Linear
import torch.nn.functional as F
from torch_geometric.nn import GraphConv, global_mean_pool
from sklearn.metrics import confusion_matrix, roc_auc_score, precision_recall_curve, auc, roc_curve
from sklearn.preprocessing import label_binarize
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from torch_geometric.nn import GINConv, GATConv, GCNConv

from sklearn.model_selection import train_test_split



class GNN(torch.nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim, dropout_prob):
        super(GNN, self).__init__()
        torch.manual_seed(12345)
        self.conv1 = GraphConv(input_dim, hidden_dim)
        self.conv2 = GraphConv(hidden_dim, hidden_dim)
        self.conv3 = GraphConv(hidden_dim, hidden_dim)
        self.lin = Linear(hidden_dim, hidden_dim)

    def forward(self, x, edge_index, batch, dropout_prob):
        # 1. Obtain node embeddings 
        x = self.conv1(x, edge_index)
        x = x.relu()
        x = self.conv2(x, edge_index)
        x = x.relu()
        x = self.conv3(x, edge_index)

        # 2. Readout layer
        x = global_mean_pool(x, batch)  # [batch_size, hidden_channels]

        # 3. Apply a final classifier
        x = F.dropout(x, p=dropout_prob, training=self.training)
        x = self.lin(x)
        
        return x


def main(directory, csv_file, num_epochs, learning_rate, tts, min_obs, bsu, hidden_dim, dropout_prob):
    # Read the CSV file
    label_df = pd.read_csv(csv_file)
    is_multi = label_df['labels'].value_counts()> min_obs
    label_df  = label_df[label_df['labels'].isin(is_multi[is_multi].index)]
    # Create mappings from string labels to integers
    unique_labels = list(label_df["labels"].unique())
    label_mapping = {label: i for i, label in enumerate(unique_labels)}
    reverse_label_mapping = {i: label for label, i in label_mapping.items()}
    
    print("Label Mapping:", label_mapping)

    # Prepare graphs and labels
    graphs = []
    labels = []
    graph_files = []
    graphs_objects = []

    for folder in os.listdir(directory):
        folder_path = os.path.join(directory, folder)
        if os.path.isdir(folder_path):
            for file in os.listdir(folder_path):
                if file.endswith(".graphml"):
                    print(file)
                    graph_path = os.path.join(folder_path, file)
                    graph = ig.Graph.Read_GraphML(graph_path)
                    print("Graph read!")
                    
                    if "fillcolor" in graph.vs.attributes():
                        print("Attributes found... Analysing...")
                        # Compute node features
                        #degrees = graph.degree()
                        betweenness = graph.betweenness()
                        # clustering = graph.transitivity_local_undirected(mode="zero")  # Clustering coefficient
                        #closeness = graph.closeness()
                        infomap_communities = graph.community_infomap()
                        community_ids = infomap_communities.membership
                        # pagerank = graph.pagerank()  # Replace eigenvector centrality
        
                        # Optionally include in-degree and out-degree
                        in_degrees = graph.indegree()
                        out_degrees = graph.outdegree()
                        
                        
                        # Access the edge attribute "d7" and compute binary encoding
                        if "arrowhead" in graph.es.attributes():
                            print("Getting MoR per edge and accummulating...")
                            edge_d7 = graph.es["arrowhead"]  # Get "d7" values
                            edge_binary = [1 if d7_value == '"vee"' else -1 if d7_value == '"tee"' else 0 for d7_value in edge_d7]  # Encode as binary
                            edge_index = graph.get_edgelist()
                            node_binary_sum = [0] * graph.vcount()  # Initialize binary sum for nodes
            
                            for (src, tgt), binary_value in zip(edge_index, edge_binary):
                                node_binary_sum[src] += binary_value
                                node_binary_sum[tgt] += binary_value
                        else:
                            print(f"Edge attribute 'arrowhead' not found in graph: {file}")
                            time.sleep(1)
                            node_binary_sum = [0] * graph.vcount() 
                        
                        # Access the node colors
                        # Convert node colors to {-1, 1}
                        node_colors = graph.vs["fillcolor"]

                        print("Getting node colours as gene state...")
                        node_binary_colors = [
                        1 if color == 'lavender' else -1 if color == 'mistyrose' else 0
                            for color in node_colors
                        ]
 
                        # Create feature matrix
                        feature_matrix = list(zip(node_binary_colors, betweenness,community_ids,in_degrees,out_degrees, node_binary_sum))

                        # Standardize features
                        scaler = StandardScaler()
                        feature_matrix = scaler.fit_transform(feature_matrix)
                        num_features = feature_matrix.shape[1]
                        # Add standardized features to the Data object
                        data = Data(
                            x=torch.tensor(feature_matrix, dtype=torch.float),
                            edge_index=torch.tensor(graph.get_edgelist(), dtype=torch.long).t().contiguous(),
                            file_name=file
                        )
                        
                        # Map and assign label
                        filename = file
                        if filename in label_df["filename"].values:
                            string_label = label_df.loc[label_df["filename"] == filename, "labels"].values[0]
                            numeric_label = label_mapping[string_label]
                            data.y = torch.tensor([numeric_label], dtype=torch.long)
                            graphs.append(data)
                            graphs_objects.append(graph)
                            labels.append(numeric_label)
                            graph_files.append(file)
                        else:
                            print(f"Warning: {filename} not found. Skipping.")
                            
                    else:
                        print("No color attribute found in the graph.")
                        time.sleep(1)
 

    print(f"Total Number of graphs: {len(graphs):.4f}")
    time.sleep(1)

    # Assuming `graphs` is your list of graph data objects and `labels` contains the corresponding labels
    X = graphs  # Features (graphs)
    y = labels  # Labels
    
    # Perform stratified split with class balance and a fixed random seed
    train_graphs, test_graphs, train_labels, test_labels = train_test_split(
        X, y, test_size=1-tts, random_state=42, stratify=y
    )
    
    print(f"Training set size: {len(train_graphs)}")
    print(f"Test set size: {len(test_graphs)}")
    

    train_loader = DataLoader(train_graphs, batch_size=bsu, shuffle=True)
    test_loader = DataLoader(test_graphs, batch_size=bsu, shuffle=False)

    # Initialize model, optimizer, and loss function
    model = GNN(input_dim=num_features, hidden_dim=hidden_dim, output_dim=len(unique_labels), dropout_prob = dropout_prob)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    loss_fn = torch.nn.CrossEntropyLoss()

    # Training loop
    for epoch in range(num_epochs):
        model.train()
        total_loss = 0
        for data in train_loader:
            optimizer.zero_grad()
            out = model(data.x, data.edge_index, data.batch, dropout_prob)
            loss = loss_fn(out, data.y)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm = 1)
            optimizer.step()
            total_loss += loss.item()
        print(f"Epoch {epoch + 1}/{num_epochs}, Loss: {total_loss:.4f}")

    # Validation
    model.eval()
    y_true, y_pred, y_pred_prob, embeddings = [], [], [], []
    aligned_graphs = []
    test_graph_files = [data.file_name for data in test_loader.dataset]
    print(len(test_graph_files))
    time.sleep(1)
    with torch.no_grad():
        for data, graph_file in zip(test_loader, test_graph_files):
            out = model(data.x, data.edge_index, data.batch, dropout_prob)
            embeddings.append(out.cpu().numpy())
            y_true.extend(data.y.tolist())
            y_pred.extend(out.argmax(dim=1).tolist())
            y_pred_prob.append(F.softmax(out, dim=1).cpu())  # Collect probabilities
  
    
    aligned_graphs = [graphs_objects[graph_files.index(f)] for f in test_graph_files]

    # Concatenate probabilities for AUROC
    y_pred_prob = torch.cat(y_pred_prob, dim=0).numpy()


    # Confusion matrix
    y_true_strings = [reverse_label_mapping[y] for y in y_true]
    y_pred_strings = [reverse_label_mapping[y] for y in y_pred]

    cm = confusion_matrix(y_true_strings, y_pred_strings, labels=unique_labels)
    print("Confusion Matrix:\n", cm)
    
    # Plot confusion matrix as a heatmap
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", xticklabels=unique_labels, yticklabels=unique_labels)
    plt.title("Confusion Matrix Heatmap")
    plt.xlabel("Predicted Labels")
    plt.ylabel("True Labels")
    plt.savefig('CM.png', dpi=300,bbox_inches='tight')
    plt.show()
    
    # AUROC

    # Ensure labels are binarized for multiclass AUROC
    y_true_binarized = label_binarize(y_true, classes=list(range(len(unique_labels))))
    y_pred = y_pred_prob.argmax(axis=1)  # Convert probabilities to predicted class labels
    
    # Calculate misclassification
    total_graphs = len(y_true)
    misclassified = (y_pred != y_true).sum()
    misclassification_percentage = (misclassified / total_graphs) * 100
    
    # Print misclassification percentage
    print(f"Total Graphs: {total_graphs}")
    print(f"Misclassified Graphs: {misclassified}")
    print(f"Misclassification Percentage: {misclassification_percentage:.2f}%")

    # Compute ROC curve and AUC for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(len(unique_labels)):
        fpr[i], tpr[i], _ = roc_curve(y_true_binarized[:, i], y_pred_prob[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])
    
    # Plot ROC curve for each class
    plt.figure(figsize=(10, 8))
    for i in range(len(unique_labels)):
        plt.plot(fpr[i], tpr[i], label=f"Class {reverse_label_mapping[i]} (AUC = {roc_auc[i]:.2f})")
    
    # Plot diagonal
    plt.plot([0, 1], [0, 1], "k--", label="Random Guess")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curve (One-vs-Rest)")
    plt.legend(loc="lower right")
    plt.savefig('ROC.png',dpi=300,bbox_inches='tight')
    plt.show()
    
    # PCA Visualization
    def plot_pca(embeddings, true_labels, predicted_labels):
        """
        Plots two PCA visualizations:
        1. By predicted class (string labels).
        2. By correctness of predictions (red = incorrect, blue = correct).
        """
        sns.set(style="whitegrid")

        # Ensure labels and embeddings align
        embeddings = np.array(embeddings)
        true_labels = np.array(true_labels)
        predicted_labels = np.array(predicted_labels)

        # Perform PCA
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(embeddings)
        explained_variance = pca.explained_variance_ratio_ * 100  # Percent variance explained

        # Scatter plot for predicted classes
        plt.figure(figsize=(12, 8))
        unique_labels = np.unique(true_labels)
        palette = sns.color_palette("tab10", len(unique_labels))

        for label, color in zip(unique_labels, palette):
                      idx = np.where(true_labels == label)[0]
                      plt.scatter(
                          pca_result[idx, 0],
                          pca_result[idx, 1],
                          label=f"{label} ({len(idx)} graphs)",
                          color=color,
                          s=100,
                          alpha=0.75,
                          edgecolor="k",
                          linewidth=0.5,
                      )
          
        plt.title("PCA of Graph Embeddings by True Classes", fontsize=16)
        plt.xlabel(f"PC1 ({explained_variance[0]:.2f}% Variance Explained)", fontsize=14)
        plt.ylabel(f"PC2 ({explained_variance[1]:.2f}% Variance Explained)", fontsize=14)
        plt.legend(title="True Classes", fontsize=12, title_fontsize=14, loc="best", markerscale=1.2)
        plt.grid(visible=True, linestyle="--", alpha=0.5)
        plt.tight_layout()
        plt.savefig("PCA_true_classes.png",dpi=300,bbox_inches='tight')
        plt.show()

        plt.figure(figsize=(12, 8))

        for label, color in zip(unique_labels, palette):
            idx = np.where(predicted_labels == label)[0]
            plt.scatter(
                pca_result[idx, 0],
                pca_result[idx, 1],
                label=f"{label} ({len(idx)} graphs)",
                color=color,
                s=100,
                alpha=0.75,
                edgecolor="k",
                linewidth=0.5,
            )

        plt.title("PCA of Graph Embeddings by Predicted Classes", fontsize=16)
        plt.xlabel(f"PC1 ({explained_variance[0]:.2f}% Variance Explained)", fontsize=14)
        plt.ylabel(f"PC2 ({explained_variance[1]:.2f}% Variance Explained)", fontsize=14)
        plt.legend(title="Predicted Classes", fontsize=12, title_fontsize=14, loc="best", markerscale=1.2)
        plt.grid(visible=True, linestyle="--", alpha=0.5)
        plt.tight_layout()
        plt.savefig("PCA_predicted_classes.png",dpi=300,bbox_inches='tight')
        plt.show()

        # Scatter plot for prediction correctness
        plt.figure(figsize=(12, 8))
        correct = true_labels == predicted_labels
        incorrect = ~correct

        plt.scatter(
            pca_result[correct, 0],
            pca_result[correct, 1],
            label=f"Correct ({correct.sum()} graphs)",
            color="blue",
            s=100,
            alpha=0.75,
            edgecolor="k",
            linewidth=0.5,
        )
        plt.scatter(
            pca_result[incorrect, 0],
            pca_result[incorrect, 1],
            label=f"Incorrect ({incorrect.sum()} graphs)",
            color="red",
            s=100,
            alpha=0.75,
            edgecolor="k",
            linewidth=0.5,
        )

        plt.title("PCA of Graph Embeddings by Prediction Accuracy", fontsize=16)
        plt.xlabel(f"PC1 ({explained_variance[0]:.2f}% Variance Explained)", fontsize=14)
        plt.ylabel(f"PC2 ({explained_variance[1]:.2f}% Variance Explained)", fontsize=14)
        plt.legend(title="Prediction Accuracy", fontsize=12, title_fontsize=14, loc="best", markerscale=1.2)
        plt.grid(visible=True, linestyle="--", alpha=0.5)
        plt.tight_layout()
        plt.savefig("PCA_prediction_accuracy.png",dpi=300,bbox_inches='tight')
        plt.show()


    # Update the call to `plot_pca` in `main`:
    plot_pca(np.vstack(embeddings), y_true_strings, y_pred_strings)
    
    def plot_pca_interactive_matplotlib(embeddings, true_labels, predicted_labels, graph_files, graphs):
        """
        Creates an interactive PCA plot using Matplotlib. Clicking on a point displays the corresponding graph.
        Parameters:
            embeddings (np.array): Graph-level embeddings from the model.
            true_labels (list): True class labels (strings).
            predicted_labels (list): Predicted class labels (strings).
            graph_files (list): File names of the graphs.
            graphs (list): List of igraph objects corresponding to the graphs.
        """
        import matplotlib.pyplot as plt
        import networkx as nx
        import numpy as np
    
        # Perform PCA
        from sklearn.decomposition import PCA
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(embeddings)
        explained_variance = pca.explained_variance_ratio_ * 100
    
        # Prepare color mapping for predicted labels
        unique_labels = np.unique(true_labels)
        label_to_color = {label: plt.cm.tab10(i / len(unique_labels)) for i, label in enumerate(unique_labels)}
        point_colors = [label_to_color[label] for label in predicted_labels]
    
        # Create the PCA scatter plot
        fig, ax = plt.subplots(figsize=(12, 8))
        scatter = ax.scatter(
            pca_result[:, 0],
            pca_result[:, 1],
            c=point_colors,
            alpha=0.75,
            edgecolor="k"
        )
    
        # Add axis labels and title
        ax.set_title(f"PCA of Graph Embeddings\nPC1 ({explained_variance[0]:.2f}%), PC2 ({explained_variance[1]:.2f}%)")
        ax.set_xlabel(f"PC1 ({explained_variance[0]:.2f}% Variance Explained)")
        ax.set_ylabel(f"PC2 ({explained_variance[1]:.2f}% Variance Explained)")
    
        # Add legend
        legend_handles = [
            plt.Line2D([0], [0], marker='o', color=color, markersize=10, label=label)
            for label, color in label_to_color.items()
        ]
        ax.legend(handles=legend_handles, title="Predicted Labels", loc="best")
    
        def on_click(event):
            if event.inaxes == ax:
                # Get the closest point
                distances = np.sqrt(
                    (pca_result[:, 0] - event.xdata) ** 2 + (pca_result[:, 1] - event.ydata) ** 2
                )
                idx = np.argmin(distances)
                selected_file = graph_files[idx]
                print(f"Selected graph: {selected_file}")
        
                # Get the selected graph (igraph object)
                graph = graphs[idx]
        
                # Extract the "fillcolor" node attribute
                if "fillcolor" in graph.vs.attributes():
                    node_colors = graph.vs["fillcolor"]  # Use fillcolor for node coloring
                else:
                    node_colors = "lightblue"  # Default color if "fillcolor" is not found
        
                # Get node names (if available) or use indices
                node_labels = graph.vs["id"] if "id" in graph.vs.attributes() else [str(i) for i in range(graph.vcount())]
        
                # Convert the igraph graph to networkx
                nx_graph = nx.DiGraph()
                nx_graph.add_nodes_from(range(graph.vcount()))  # Add nodes
                nx_graph.add_edges_from(graph.get_edgelist())  # Add edges

                # Create a list of positions for networkx to plot
                pos = nx.bfs_layout(nx_graph,0)

                # Display the graph using networkx with node colors from "fillcolor"
                fig, ax_graph = plt.subplots(figsize=(8, 6))
                nx.draw(
                    nx_graph,
                    with_labels=True,
                    pos = pos,
                    labels=dict(enumerate(node_labels)),  # Use node labels
                    node_color=node_colors,  # Use fillcolor for node colors
                    edge_color="gray",
                    node_size=500,
                    font_size=10,
                    arrows = True,
                    ax = ax_graph
                )
                
                # Add title with graph file name
                ax_graph.set_title(f"Graph: {selected_file}", fontsize = 12)
                plt.show()


        # Connect the click event to the callback
        fig.canvas.mpl_connect("button_press_event", on_click)
        plt.savefig("PCA_interactive.png",dpi=300,bbox_inches='tight')
        plt.show()




    
        # Update the call to `plot_pca` in `main`:
    plot_pca_interactive_matplotlib(
      embeddings=np.vstack(embeddings),
      true_labels=y_true_strings,
      predicted_labels=y_pred_strings,
      graph_files=test_graph_files,
      graphs=aligned_graphs
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Train GNN on graph datasets")
    parser.add_argument("directory", type=str, help="Directory containing folders with .graphml files")
    parser.add_argument("csv_file", type=str, help="CSV file with labels")
    parser.add_argument("num_epochs", type=int, help="Number of training epochs")
    parser.add_argument("learning_rate", type=float, help="Learning rate for training")
    parser.add_argument("tts", type=float, help="Train-test split ratio")
    parser.add_argument("min_obs", type=int, help="Number of min obs per class")
    parser.add_argument("bsu", type=int, help="Batch size for training/testing")
    parser.add_argument("hidden_dim", type=int, help="Hidden dimension")
    parser.add_argument("dropout_prob", type=float, help="Dropout Probability")
    args = parser.parse_args()
    #python3 GNN3.py "/Users/harrytriantafyllidis/Desktop/Shiny_app_2/CCLE_TP53..7157./cancers_all_cell_lines/opt" "/Users         /harrytriantafyllidis/Desktop/Shiny_app_2/CCLE_TP53..7157./cancers_all_cell_lines/opt/GNN_labels.csv" 500 0.005 0.75 50
    # vee activating, tee inhibitory
    
    main(args.directory, args.csv_file, args.num_epochs, args.learning_rate, args.tts, args.min_obs, args.bsu, args.hidden_dim, args.dropout_prob)

