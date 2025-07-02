from torch_geometric.loader import DataLoader
import torch.nn.functional as F
from sklearn.metrics import confusion_matrix, roc_auc_score, precision_recall_curve, auc, roc_curve
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.model_selection import train_test_split
import os, sys, time
import pandas as pd
import igraph as ig
import torch
import shutil
from sklearn.utils.class_weight import compute_class_weight
from torch_geometric.data import Data
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from torch_geometric.nn import GATv2Conv, global_mean_pool, BatchNorm
from torch.nn import Linear

from torch_geometric.nn import GATv2Conv
from sklearn.metrics import classification_report
from xgboost import XGBClassifier
import torch.nn as nn

class GNN(torch.nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim, dropout_prob, heads, hop_order):
        super(GNN, self).__init__()
        torch.manual_seed(12345)
        edge_dim=3
        self.hop_order = hop_order

        self.edge_encoder = torch.nn.Linear(3, edge_dim)

        self.gat1 = GATv2Conv(input_dim, hidden_dim, heads=heads, dropout=dropout_prob, edge_dim=edge_dim)
        self.bn1 = BatchNorm(hidden_dim * heads)

        self.gat2 = GATv2Conv(hidden_dim * heads, hidden_dim, heads=heads, dropout=dropout_prob, edge_dim=edge_dim)
        self.bn2 = BatchNorm(hidden_dim * heads)

        self.gat3 = GATv2Conv(hidden_dim * heads, hidden_dim, heads=heads, dropout=dropout_prob, edge_dim=edge_dim)
        self.bn3 = BatchNorm(hidden_dim * heads)

        self.lin = Linear(hidden_dim * heads, output_dim)

        # Neighborhood-aware GOI encoder (small GAT)
        self.goi_gnn = GATv2Conv(
            input_dim, input_dim, heads=heads, concat=True,
            dropout=0.5, edge_dim=edge_dim
        )
        self.goi_proj = nn.Linear(input_dim * heads, input_dim)

        # GOI-specific learnable bias vector
        self.goi_bias = torch.nn.Parameter(torch.zeros(input_dim))
        self.goi_attention_scalar = torch.nn.Parameter(torch.tensor(1.0))

    def forward(self, x, edge_index, edge_attr, batch, dropout_prob, goi_mask=None):
        from torch_geometric.utils import k_hop_subgraph

        if goi_mask is not None:
            goi_mask = goi_mask.to(x.device)
            spotlight_indices = goi_mask.nonzero(as_tuple=False).flatten()
        
            x = x.clone()
        
            for idx in spotlight_indices:
                goi_index = idx.item()
        
                # Extract k-hop subgraph
                subset, sub_edge_index, mapping, sub_edge_mask = k_hop_subgraph(
                    node_idx=goi_index,
                    num_hops=self.hop_order,
                    edge_index=edge_index,
                    relabel_nodes=True,
                    num_nodes=x.size(0)
                )
        
                x_sub = x[subset]
                edge_attr_sub = edge_attr[sub_edge_mask]
        
                # GATv2 with concat=True → output dim = input_dim * heads
                x_sub_updated = self.goi_gnn(x_sub, sub_edge_index, edge_attr_sub)
                x_sub_updated = self.goi_proj(x_sub_updated)  # back to input_dim
        
                # Inject back into full graph
                x[subset] = x_sub_updated
        
                # Bias injection (same for all spotlighted subgraph nodes)
                x[subset] += self.goi_bias * self.goi_attention_scalar

        # 3. Encode edge features
        edge_attr = self.edge_encoder(edge_attr)
    
        # 4. Pass through full GATv2 stack
        x = F.elu(self.bn1(self.gat1(x, edge_index, edge_attr)))
        x = F.elu(self.bn2(self.gat2(x, edge_index, edge_attr) + x))
        x = F.elu(self.bn3(self.gat3(x, edge_index, edge_attr) + x))
        x = global_mean_pool(x, batch)
        x = F.dropout(x, p=0.5, training=self.training)
        return self.lin(x)

# --- BASELINE CLASSIFIER (GRAPH-LEVEL) ---
def run_baseline_classifier(graphs, labels, tts, reverse_label_mapping=None, output_path="baseline_roc_auc.png"):
    features = []
    
    for graph in graphs:
        degrees = graph.degree()
        avg_degree = np.mean(degrees)
        avg_betweenness = np.mean(graph.betweenness())
        node_colors = graph.vs["fillcolor"]
        node_binary = [1 if c == "lavender" else -1 if c == "mistyrose" else 0 for c in node_colors]
        edge_d7 = graph.es["arrowhead"]
        edge_binary = [1 if d7 == '"vee"' else -1 if d7 == '"tee"' else 0 for d7 in edge_d7]

        percent_inhibitory = np.sum(np.array(edge_binary) == -1) / len(edge_binary)
        percent_activating = np.sum(np.array(edge_binary) == 1) / len(edge_binary)

        features.append([
            np.mean(node_binary),
            np.mean(edge_binary),
            avg_betweenness
        ])

    X = np.array(features)
    y = np.array(labels)

    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=1-tts, stratify=y, random_state=42)

    # Class Weights
    classes = np.unique(y_train)
    class_weights = compute_class_weight(class_weight='balanced', classes=classes, y=y_train)
    class_weight_dict = {cls: weight for cls, weight in zip(classes, class_weights)}
    print("Class Weights:", class_weight_dict)

    # Classifier with sample weights
    clf = XGBClassifier(use_label_encoder=False, eval_metric='mlogloss')
    clf.fit(X_train, y_train, sample_weight=np.array([class_weight_dict[cls] for cls in y_train]))

    y_pred = clf.predict(X_test)
    y_pred_proba = clf.predict_proba(X_test)

    misclassified_pct = 100 * (y_pred != y_test).sum() / len(y_test)

    print("\n--- BASELINE CLASSIFIER REPORT ---")
    if reverse_label_mapping:
        y_test_str = [reverse_label_mapping[y] for y in y_test]
        y_pred_str = [reverse_label_mapping[y] for y in y_pred]
        print(classification_report(y_test_str, y_pred_str))
        print("Confusion Matrix:")
        print(confusion_matrix(y_test_str, y_pred_str))
    else:
        print(classification_report(y_test, y_pred))
        print("Confusion Matrix:")
        print(confusion_matrix(y_test, y_pred))

    print(f"Misclassification rate: {misclassified_pct:.2f}%")

    # --- ROC AUC Curve ---
    plt.figure()

    # Multi-class needs binarizing
    y_test_binarized = label_binarize(y_test, classes=classes)
    n_classes = y_test_binarized.shape[1]

    for i in range(n_classes):
        fpr, tpr, _ = roc_curve(y_test_binarized[:, i], y_pred_proba[:, i])
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, label=f'Class {reverse_label_mapping[classes[i]] if reverse_label_mapping else classes[i]} (AUC = {roc_auc:.2f})')

    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC AUC - Baseline Classifier')
    plt.legend(loc="lower right")

    plt.tight_layout()
    plt.savefig('RESULTS/Baseline_ROC.png', dpi=600, bbox_inches='tight')
    plt.close()
    
    print(f"ROC AUC curve saved")
    
def main(directory, csv_file, num_epochs, learning_rate, tts, min_obs, bsu, hidden_dim, dropout_prob, heads_user, hops_user, optuna):
    import multiprocessing
    #torch.autograd.set_detect_anomaly(True)
    os.makedirs("RESULTS", exist_ok=True)
    num_workers = multiprocessing.cpu_count()
    label_df = pd.read_csv(csv_file)
    
    import matplotlib.pyplot as plt
    # Report class distribution before subsampling
    print("\n--- Class Distribution Before Subsampling ---")
    class_counts = label_df["labels"].value_counts()
    print(class_counts.to_frame(name="Sample Count"))
    time.sleep(2)
    
    # Automatically infer max_per_class from the largest class
    ##########################################################
    max_per_class = class_counts.max() # change this to ENFORCE subsampling, e.g:
    #max_per_class = 1000
    ##########################################################
    print(f"\nAuto-inferred max_per_class = {max_per_class} (from dominant class)")
    # Compute original and subsampled class distributions
    original_dist = label_df["labels"].value_counts().sort_index()
    label_df = label_df.groupby("labels", group_keys=False).apply(lambda x: x.sample(n=min(len(x), max_per_class), random_state=42))
    subsampled_dist = label_df["labels"].value_counts().sort_index()
    
    # Labels with counts
    original_labels = [f"{label} ({count})" for label, count in zip(original_dist.index, original_dist.values)]
    subsampled_labels = [f"{label} ({count})" for label, count in zip(subsampled_dist.index, subsampled_dist.values)]
    
    # Create side-by-side pie charts
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    
    # --- Left: Before Subsampling ---
    wedges1, _ = axes[0].pie(original_dist.values, startangle=140)
    axes[0].set_title("Before Subsampling", fontsize=14)
    axes[0].legend(wedges1, original_labels, title="Classes", loc="center left", bbox_to_anchor=(1.05, 0.5), fontsize=10)
    
    # --- Right: After Subsampling ---
    wedges2, _ = axes[1].pie(subsampled_dist.values, startangle=140)
    axes[1].set_title("After Subsampling", fontsize=14)
    axes[1].legend(wedges2, subsampled_labels, title="Classes", loc="center left", bbox_to_anchor=(1.05, 0.5), fontsize=10)
    
    plt.tight_layout()
    plt.savefig("RESULTS/class_balance_pies.png", dpi=600, bbox_inches='tight')
    if optuna == "False":
        plt.show()
    
    is_multi = label_df['labels'].value_counts() > min_obs
    label_df = label_df[label_df['labels'].isin(is_multi[is_multi].index)]

    unique_labels = list(label_df["labels"].unique())
    label_mapping = {label: i for i, label in enumerate(unique_labels)}
    reverse_label_mapping = {i: label for label, i in label_mapping.items()}

    print("Label Mapping:", label_mapping)

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

                        # Get community memberships
                        infomap_communities = graph.community_infomap()
                        community_ids = infomap_communities.membership
                       
                        # Initialize betweenness centrality array
                        community_betweenness = np.zeros(graph.vcount())
                        
                        # Compute betweenness centrality separately for each community
                        for community_id in set(community_ids):
                            # Nodes belonging to this community
                            community_nodes = [idx for idx, cid in enumerate(community_ids) if cid == community_id]
                        
                            # Extract community subgraph
                            subgraph = graph.subgraph(community_nodes)
                        
                            # Betweenness centrality within the community
                            #sub_betweenness =  subgraph.betweenness()
                            #sub_betweenness =  subgraph.hub_score()
                            #sub_betweenness = subgraph.authority_score()
                            sub_betweenness = subgraph.pagerank()
                            
                            # Assign community-specific betweenness to original nodes
                            for idx, original_node_idx in enumerate(community_nodes):
                                community_betweenness[original_node_idx] = sub_betweenness[idx]
                        
                        
                        in_degrees = graph.indegree()
                        out_degrees = graph.outdegree()

                        if "arrowhead" in graph.es.attributes():
                            print("Getting MoR per edge and accummulating...")
                            edge_d7 = graph.es["arrowhead"]
                            # One-hot encode edge types: [activating, inhibitory, neutral]
                            edge_attr_list = []
                            for d7_value in edge_d7:
                                if d7_value == '"vee"':
                                    edge_attr_list.append([1, 0, 0])  # Activating
                                elif d7_value == '"tee"':
                                    edge_attr_list.append([0, 1, 0])  # Inhibitory
                                else:
                                    edge_attr_list.append([0, 0, 1])  # Neutral/Other
                            edge_index = graph.get_edgelist()
                            node_in_mor_sum = [0] * graph.vcount()
                            node_out_mor_sum = [0] * graph.vcount()
                            
                            for (src, tgt), attr in zip(graph.get_edgelist(), edge_attr_list):
                                if attr == [1, 0, 0]:       # Activating
                                    binary_value = 1
                                elif attr == [0, 1, 0]:     # Inhibitory
                                    binary_value = -1
                                else:                       # Neutral
                                    binary_value = 0
                            
                                node_out_mor_sum[src] += binary_value  # what the node does
                                node_in_mor_sum[tgt]  += binary_value  # what the node receives

                        else:
                            print(f"Edge attribute 'arrowhead' not found in graph: {file}")
                            node_in_mor_sum = [0] * graph.vcount()
                            node_out_mor_sum = [0] * graph.vcount()

                        print("Getting node colours as gene state...")
                        node_colors = graph.vs["fillcolor"]

                        # Binary feature: 1 = lavender, -1 = mistyrose, 0 = other
                        node_binary_colors = [1 if color == 'lavender' else -1 if color == 'mistyrose' else 0 for color in graph.vs["fillcolor"]]
                        
                        # Convert features to numpy arrays for consistent shape and weight centrality per node activity
                        community_betweenness = [(cb * color) for cb, color in zip(community_betweenness, node_binary_colors)]
                        community_betweenness = np.array(community_betweenness).reshape(-1, 1)
                        
                        in_degrees = np.array(in_degrees).reshape(-1, 1)
                        out_degrees = np.array(out_degrees).reshape(-1, 1)
                        node_in_mor_sum = np.array(node_in_mor_sum).reshape(-1, 1)
                        node_out_mor_sum = np.array(node_out_mor_sum).reshape(-1, 1)
                        
                        # Standardize selected features (not the binary colors)
                        scaler = StandardScaler()
                        scaled_features = scaler.fit_transform(np.hstack([
                            community_betweenness,
                            in_degrees,
                            out_degrees,
                            node_in_mor_sum,
                            node_out_mor_sum
                        ]))
                        
                        # Combine raw binary color feature with scaled numerical features
                        feature_matrix = np.hstack([
                            np.array(node_binary_colors).reshape(-1, 1),  # Keep binary color raw
                            scaled_features
                        ])
  
                        data = Data(
                            x=torch.tensor(feature_matrix, dtype=torch.float),
                            edge_index=torch.tensor(graph.get_edgelist(), dtype=torch.long).t().contiguous(),
                            edge_attr=torch.tensor(edge_attr_list, dtype=torch.float),
                            file_name=file
                        )
                        #genes_of_interest = ["TP53"]
                        genes_of_interest = ["TP53", "MYC"]
                        goi_mask = torch.tensor(
                            [1 if name in genes_of_interest else 0 for name in graph.vs["name"]],
                            dtype=torch.float
                        )
                        data.goi_mask = goi_mask
                        
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
                        time.sleep(2)


    run_baseline_classifier(graphs_objects, labels, tts, reverse_label_mapping)
    
    print(f"Total Number of graphs: {len(graphs):.4f}")
    time.sleep(1)
    
    if optuna:
      # 3-way split: train, val, test
      print("3-way split: train, val, test:")
      time.sleep(5)
      train_graphs_full, test_graphs, train_labels_full, test_labels = train_test_split(
        graphs, labels, test_size=1-tts, stratify=labels, random_state=42
      )
      train_graphs, val_graphs, train_labels, val_labels = train_test_split(
        train_graphs_full, train_labels_full, test_size=1-tts, stratify=train_labels_full, random_state=42
      )
    else:
    # Standard 2-way split
      train_graphs, test_graphs, train_labels, test_labels = train_test_split(
        graphs, labels, test_size=1 - tts, stratify=labels, random_state=42
      )
    print(f"Training set size: {len(train_graphs)}")
    print(f"Test set size: {len(test_graphs)}")

    # Print the size of each dataset
    print(f"Training set size: {len(train_graphs)}")
    print(f"Test set size: {len(test_graphs)}")


    train_loader = DataLoader(train_graphs, batch_size=bsu, shuffle=True)
    test_loader = DataLoader(test_graphs, batch_size=1, shuffle=False)

    model = GNN(
        input_dim=feature_matrix.shape[1],
        hidden_dim=hidden_dim,
        output_dim=len(unique_labels),
        dropout_prob=dropout_prob,
        heads = heads_user,
        hop_order=hops_user # will encode GOI and all nodes within hop_order hops
    )

    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

    
    from sklearn.utils.class_weight import compute_class_weight
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    # Convert to tensor for torch.unique
    classes = torch.unique(torch.tensor(train_labels)).cpu().numpy()
    
    # Compute class weights
    class_weights = compute_class_weight(class_weight='balanced',
                                         classes=classes,
                                         y=np.array(train_labels))  # <- np.array here
    
    # Convert to tensor for PyTorch loss
    class_weights = torch.tensor(class_weights, dtype=torch.float).to(device)
    loss_fn = torch.nn.CrossEntropyLoss(weight=class_weights)
    
    from tqdm import tqdm  # Import tqdm for progress bar

    # Initialize tqdm progress bar for all epochs
    progress_bar = tqdm(range(num_epochs), desc="Training Progress", unit="epoch")
    
    for epoch in progress_bar:
        model.train()
        total_loss = 0
        for data in train_loader:
            optimizer.zero_grad()
            out = model(data.x, data.edge_index, data.edge_attr, data.batch, dropout_prob,goi_mask=data.goi_mask)
            loss = loss_fn(out, data.y)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1)
            optimizer.step()
            total_loss += loss.item()
        
        #avg_loss = total_loss / len(train_loader)
        
        # Update progress bar with current epoch and loss
        progress_bar.set_description(f"Epoch {epoch + 1}/{num_epochs}")
        progress_bar.set_postfix(loss=total_loss)
        if hasattr(model, 'goi_attention_scalar') and hasattr(model, 'goi_bias'):
          print(f"[Epoch {epoch}] Spotlight attention scalar: {model.goi_attention_scalar.item():.4f}")
          print(f"[Epoch {epoch}] Spotlight bias norm: {torch.norm(model.goi_bias).item():.4f}")
    
    progress_bar.close()
        
    # Visualize learned edge type embeddings
    
    encoder_weights = model.edge_encoder.weight.detach().cpu().numpy()
    encoder_weights = encoder_weights.T  # [3, edge_dim]
    
    edge_type_labels = ['Activating', 'Inhibitory', 'Unknown']
    
    plt.figure(figsize=(10, 3))
    sns.heatmap(encoder_weights, annot=True, cmap="coolwarm",
                xticklabels=[f'dim{i}' for i in range(encoder_weights.shape[1])],
                yticklabels=edge_type_labels)
    plt.title("Learned Edge Type Embeddings")
    plt.xlabel("Edge Feature Dimension")
    plt.ylabel("Edge Type")
    plt.tight_layout()
    plt.savefig("RESULTS/edge_learning.png",dpi=600)
    #plt.show()
    
    y_true, y_pred, y_pred_prob, embeddings = [], [], [], []
    aligned_graphs = []
    test_graph_files = [data.file_name for data in test_loader.dataset]
    print(len(test_graph_files))
    time.sleep(1)
    
    eval_graphs = val_graphs if optuna else test_graphs
    eval_graph_files = [data.file_name for data in eval_graphs]
    eval_loader = DataLoader(eval_graphs, batch_size=1, shuffle=False)
    
    model.eval()
    with torch.no_grad():
        for data, graph_file in zip(eval_loader, eval_graph_files):
            out = model(data.x, data.edge_index, data.edge_attr, data.batch, dropout_prob,goi_mask=data.goi_mask)
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

    print("\n--- GNN CLASSIFIER REPORT ---")
    print(classification_report(y_true_strings, y_pred_strings))


    cm = confusion_matrix(y_true_strings, y_pred_strings, labels=unique_labels)
    print("Confusion Matrix:\n", cm)
    
    # Plot confusion matrix as a heatmap
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", xticklabels=unique_labels, yticklabels=unique_labels)
    plt.title("Confusion Matrix Heatmap")
    plt.xlabel("Predicted Labels")
    plt.ylabel("True Labels")
    plt.savefig('RESULTS/CM.png', dpi=600,bbox_inches='tight')
    #plt.show()
    
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

    if len(unique_labels) == 2:
      # Special handling for binary classification
      fpr, tpr, _ = roc_curve(y_true, y_pred_prob[:, 1])
      roc_auc = auc(fpr, tpr)
      plt.figure(figsize=(10, 8))
      plt.plot(fpr, tpr, label=f"Class {reverse_label_mapping[1]} (AUC = {roc_auc:.2f})")
      plt.plot([0, 1], [0, 1], "k--", label="Random Guess")
      plt.xlim([0.0, 1.0])
      plt.ylim([0.0, 1.05])
      plt.xlabel("False Positive Rate")
      plt.ylabel("True Positive Rate")
      plt.title("ROC Curve (Binary Classification)")
      plt.legend(loc="lower right")
      plt.savefig('RESULTS/ROC_binary.png', dpi=600, bbox_inches='tight')
      #plt.show()
    else:
      # Multiclass case (existing logic)
      fpr = dict()
      tpr = dict()
      roc_auc = dict()
      for i in range(len(unique_labels)):
          fpr[i], tpr[i], _ = roc_curve(y_true_binarized[:, i], y_pred_prob[:, i])
          roc_auc[i] = auc(fpr[i], tpr[i])
      plt.figure(figsize=(10, 8))
      for i in range(len(unique_labels)):
          plt.plot(fpr[i], tpr[i], label=f"Class {reverse_label_mapping[i]} (AUC = {roc_auc[i]:.2f})")
      plt.plot([0, 1], [0, 1], "k--", label="Random Guess")
      plt.xlim([0.0, 1.0])
      plt.ylim([0.0, 1.05])
      plt.xlabel("False Positive Rate")
      plt.ylabel("True Positive Rate")
      plt.title("ROC Curve (One-vs-Rest)")
      plt.legend(loc="lower right")
      plt.savefig('RESULTS/ROC_multiclass.png', dpi=600, bbox_inches='tight')
      #plt.show()

    
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
        plt.savefig("RESULTS/PCA_true_classes.png",dpi=600,bbox_inches='tight')
        #plt.show()

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
        plt.savefig("RESULTS/PCA_predicted_classes.png",dpi=600,bbox_inches='tight')
        #plt.show()

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
        plt.savefig("RESULTS/PCA_prediction_accuracy.png",dpi=600,bbox_inches='tight')
        #plt.show()


    # Update the call
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
                if optuna == "False":
                    plt.show()


        # Connect the click event to the callback
        fig.canvas.mpl_connect("button_press_event", on_click)
        plt.savefig("RESULTS/PCA_interactive.png",dpi=600,bbox_inches='tight')
        if optuna == "False":
           plt.show()


    # Update the call
    plot_pca_interactive_matplotlib(
      embeddings=np.vstack(embeddings),
      true_labels=y_true_strings,
      predicted_labels=y_pred_strings,
      graph_files=test_graph_files,
      graphs=aligned_graphs
    )
    return misclassification_percentage

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
    parser.add_argument("heads_user", type=int, help="Number of Attention Heads")
    parser.add_argument("hops_user", type=int, help="Hop order in Spotlight mini-GNN")
    parser.add_argument("optuna", type=lambda x: x.lower() == "true", help="Enable Optuna 3-way split logic (True/False)")
    args = parser.parse_args()

    
    main(args.directory, args.csv_file, args.num_epochs, args.learning_rate, args.tts, args.min_obs, args.bsu, args.hidden_dim, args.dropout_prob, args.heads_user, args.hops_user)

