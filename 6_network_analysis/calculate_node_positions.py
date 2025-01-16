import argparse
import json
import numpy as np
import networkx as nx
from networkx.readwrite import json_graph

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Calculate node positions.')
parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
parser.add_argument('--input', dest='input', required=True, help='Path to the input JSON file')
parser.add_argument('--output', dest='output', required=True, help='Path to the output JSON file')
parser.add_argument('--k', dest='k_param', type=float, required=True, help='k parameter for spring layout')
parser.add_argument('--iter', dest='iter_param', type=int, required=True, help='Number of iterations for spring layout')
args = parser.parse_args()

num_threads = args.threads
print(f"Using {num_threads} threads")

# Load the JSON data from a file
input_file = args.input
with open(input_file, 'r') as f:
    data = json.load(f)

# Create a new graph from the JSON data
G = json_graph.node_link_graph(data)

# Calculate node positions
pos = nx.spring_layout(G, k=args.k_param, iterations=args.iter_param, scale=1)

# save pos as json
def convert_numpy_types(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj

output_file = f"{args.output}_k{args.k_param}_iter{args.iter_param}.json"
with open(output_file, 'w') as f:
    json.dump(pos, f, default=convert_numpy_types)