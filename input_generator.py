import random

def generate_random_graph(num_vertices, num_edges, output_file):
    if num_edges > num_vertices * (num_vertices - 1) // 2:
        raise ValueError("Too many edges for an undirected simple graph.")

    edges = set()
    while len(edges) < num_edges:
        u = random.randint(0, num_vertices - 1)
        v = random.randint(0, num_vertices - 1)
        if u != v:
            edge = (min(u, v), max(u, v))  # Avoid self-loops and duplicate edges
            edges.add(edge)

    with open(output_file, 'w') as f:
        f.write(f"{num_vertices}\n")
        for u, v in sorted(edges):
            f.write(f"{u} {v}\n")

    print(f"Graph with {num_vertices} vertices and {num_edges} edges written to '{output_file}'.")

# Example usage
if __name__ == "__main__":
    num_vertices = 80
    num_edges = 1600
    for i in range(10):
        output_file = f"data/huge/dense/V{num_vertices}_E{num_edges}({i}).txt"
        generate_random_graph(num_vertices, num_edges, output_file)
