import os
import igraph as ig
import pydot
import sys

def convert_dot_to_igraph(dot_path, output_path, violin_gene):
    """
    Convert a .dot file to an igraph object and save it as an .RDS-compatible file for R.
    
    Parameters:
    - dot_path: str, path to the .dot file
    - output_path: str, path to save the converted igraph object
    """
    try:
        # Parse the .dot file using pydot
        (dot_graph,) = pydot.graph_from_dot_file(dot_path)
        
        # Convert the pydot graph to an igraph object
        igraph_graph = ig.Graph(directed=dot_graph.get_type() == "digraph")
        
        # Add vertices
        
        
        dot_graph.del_node(violin_gene, index = 1)
        nodes = dot_graph.get_nodes()
        
        node_names = [node.get_name() for node in nodes if node.get_name() not in ('node', '', None)]
        
        print(node_names)
        igraph_graph.add_vertices(node_names)
        
        
        # Add edges
        edges = dot_graph.get_edges()
        edge_list = [(edge.get_source(), edge.get_destination()) for edge in edges]
        igraph_graph.add_edges(edge_list)
        
        # Save the igraph object in GraphML format (compatible with R)
        igraph_graph.write_graphml(output_path)
        print(f"Saved igraph object to: {output_path}")
        
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
