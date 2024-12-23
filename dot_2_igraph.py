import igraph as ig
import os, sys
import networkx as nx

def convert_dot_to_igraph(dot_path, output_path, violin_gene):
    """
    Convert a .dot file to an igraph object and save it as an .RDS-compatible file for R.
    
    Parameters:
    - dot_path: str, path to the .dot file
    - output_path: str, path to save the converted igraph object
    """
    try:
        # Step 1: Read the .DOT file using networkx
        dot_graph = nx.nx_pydot.read_dot(dot_path)
        
        # Step 2: Check and process node attributes (e.g., 'color')
        # for node, data in dot_graph.nodes(data=True):
        #     if "fillcolor" in data:
        #         print(f"Node {node} color: {data['fillcolor']}")  # Example output
        #     else:
        #         data["fillcolor"] = "none"  # Add default color if missing
        # 
        # Step 3: Save the graph to .graphml
        nx.write_graphml(dot_graph, output_path)
        print("Graph saved as .graphml with node attributes.")
        
    except Exception as e:
        print(f"Failed to process {dot_path}: {e}")

def scan_and_convert(directory,violin_gene):
    """
    Scan through all folders in a directory, find .dot files, and convert them to igraph objects.
    
    Parameters:
    - directory: str, root directory to scan
    """
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".dot"):
                dot_path = os.path.join(root, file)
                output_path = os.path.join(root, file.replace(".dot", ".graphml"))
                print(f"Processing {dot_path}...")
                convert_dot_to_igraph(dot_path, output_path,violin_gene)

if __name__ == "__main__":
    # Specify the directory to scan
    directory_to_scan = sys.argv[1]
    violin_gene = sys.argv[2]
    print(violin_gene)
    print(directory_to_scan)
    scan_and_convert(directory_to_scan,violin_gene)
