#code modified for vector testing
import os
import numpy as np

class Graph:
    def __init__(self, size):
        self.adj_matrix = [[0] * size for _ in range(size)]
        self.size = size
        self.vertex_data = [''] * size
        self.resid_mapping = {}
        self.vertex_coords = {}  # To store coordinates of each vertex

    def add_edge(self, u, v, weight):
        if 0 <= u < self.size and 0 <= v < self.size:
            self.adj_matrix[u][v] = weight
            self.adj_matrix[v][u] = weight  # For undirected graph

    def add_vertex_data(self, vertex, data, resid, coords):
        if 0 <= vertex < self.size:
            self.vertex_data[vertex] = data
            self.resid_mapping[data] = resid
            self.vertex_coords[data] = coords

    def dijkstra(self, start_vertex_data):
        start_vertex = self.vertex_data.index(start_vertex_data)
        distances = [float('-inf')] * self.size
        distances[start_vertex] = 0
        visited = [False] * self.size
        previous = [None] * self.size

        for _ in range(self.size):
            max_distance = float('-inf')
            u = None
            for i in range(self.size):
                if not visited[i] and distances[i] > max_distance:
                    max_distance = distances[i]
                    u = i

            if u is None:
                break

            visited[u] = True

            for v in range(self.size):
                if self.adj_matrix[u][v] != 0 and not visited[v]:
                    alt = max(distances[u], self.adj_matrix[u][v])
                    if alt > distances[v]:
                        distances[v] = alt
                        previous[v] = u

        return distances, previous

    def find_nth_highest_coupling_path(self, start_vertex_data, n):
        distances, previous = self.dijkstra(start_vertex_data)
        sorted_distances = sorted(((d, i) for i, d in enumerate(distances)), reverse=True)
        
        if n-1 >= len(sorted_distances):
            raise ValueError(f"Requested the {n}th highest value, but only {len(sorted_distances)} values are available.")
        
        nth_coupling_value, nth_coupling_index = sorted_distances[n-1]
        
        path = []
        current = nth_coupling_index
        while current is not None:
            path.insert(0, self.vertex_data[current])
            current = previous[current]

        path_resids = [self.resid_mapping[vertex] for vertex in path]
        return path_resids, nth_coupling_value

def read_average_coupling_values(file_path):
    vertex_mapping = {}
    vertex_index = 1  # Start from 1 since 0 is reserved for the source resid
    source_resid = None
    data = []

    with open(file_path, 'r') as file:
        next(file)  # Skip header
        for line in file:
            parts = line.strip().split()
            if len(parts) != 3:
                continue  # Skip any line that doesn't have exactly 3 elements

            src_resid = int(parts[0])
            tgt_resid = int(parts[1])
            avg_cpl = float(parts[2])  # Convert from eV to meV
            
            if source_resid is None:
                source_resid = src_resid
                vertex_mapping[src_resid] = 0  # Map source resid to 0
            
            if tgt_resid not in vertex_mapping:
                vertex_mapping[tgt_resid] = vertex_index
                vertex_index += 1

            data.append((src_resid, tgt_resid, avg_cpl))
    
    size = len(vertex_mapping)
    graph = Graph(size)
    graph.add_vertex_data(0, 'A', source_resid, (0, 0, 0))  # Add source resid as vertex 0

    # Add the remaining vertex data
    for resid, index in vertex_mapping.items():
        if index != 0:
            label = chr(ord('A') + index)
            coords = (0, 0, 0)  # Placeholder coordinates
            graph.add_vertex_data(index, label, resid, coords)

    # Add edges to the graph
    for src_resid, tgt_resid, avg_cpl in data:
        u = vertex_mapping[src_resid]
        v = vertex_mapping[tgt_resid]
        graph.add_edge(v, u, avg_cpl)

    return graph

def append_to_output_file(output_file, source_resid, final_target_resid, avg_coupling, rank, is_in_direction):
    if not os.path.exists(output_file):
        with open(output_file, 'w') as f:
            f.write("Source_resid\tFinal_target_resid\tAverage_cpl\tRank\tIn_Direction\n")
    
    with open(output_file, 'a') as f:
        f.write(f"{source_resid}\t{final_target_resid}\t{avg_coupling}\t{rank}\t{is_in_direction}\n")

#.........................................................................................................................................
# import os
# import sys

# class Graph:
#     def __init__(self, size):
#         self.adj_matrix = [[0] * size for _ in range(size)]
#         self.size = size
#         self.vertex_data = [''] * size
#         self.resid_mapping = {}

#     def add_edge(self, u, v, weight):
#         if 0 <= u < self.size and 0 <= v < self.size:
#             self.adj_matrix[u][v] = weight
#             self.adj_matrix[v][u] = weight  # For undirected graph

#     def add_vertex_data(self, vertex, data, resid):
#         if 0 <= vertex < self.size:
#             self.vertex_data[vertex] = data
#             self.resid_mapping[data] = resid

#     def dijkstra(self, start_vertex_data):
#         start_vertex = self.vertex_data.index(start_vertex_data)
#         distances = [float('-inf')] * self.size
#         distances[start_vertex] = 0
#         visited = [False] * self.size
#         previous = [None] * self.size

#         for _ in range(self.size):
#             max_distance = float('-inf')
#             u = None
#             for i in range(self.size):
#                 if not visited[i] and distances[i] > max_distance:
#                     max_distance = distances[i]
#                     u = i

#             if u is None:
#                 break

#             visited[u] = True

#             for v in range(self.size):
#                 if self.adj_matrix[u][v] != 0 and not visited[v]:
#                     alt = max(distances[u], self.adj_matrix[u][v])
#                     if alt > distances[v]:
#                         distances[v] = alt
#                         previous[v] = u

#         return distances, previous

#     def find_nth_highest_coupling_path(self, start_vertex_data, n):
#         distances, previous = self.dijkstra(start_vertex_data)
#         sorted_distances = sorted(((d, i) for i, d in enumerate(distances)), reverse=True)
        
#         if n-1 >= len(sorted_distances):
#             raise ValueError(f"Requested the {n}th highest value, but only {len(sorted_distances)} values are available.")
        
#         nth_coupling_value, nth_coupling_index = sorted_distances[n-1]
        
#         path = []
#         current = nth_coupling_index
#         while current is not None:
#             path.insert(0, self.vertex_data[current])
#             current = previous[current]

#         path_resids = [self.resid_mapping[vertex] for vertex in path]
#         return path_resids, nth_coupling_value

# def read_average_coupling_values(file_path):
#     vertex_mapping = {}
#     vertex_index = 1  # Start from 1 since 0 is reserved for the source resid
#     source_resid = None
#     data = []

#     with open(file_path, 'r') as file:
#         next(file)  # Skip header
#         for line in file:
#             parts = line.strip().split()
#             if len(parts) != 3:
#                 continue  # Skip any line that doesn't have exactly 3 elements

#             src_resid = int(parts[0])
#             tgt_resid = int(parts[1])
#             avg_cpl = float(parts[2])  # Convert from eV to meV
            
#             if source_resid is None:
#                 source_resid = src_resid
#                 vertex_mapping[src_resid] = 0  # Map source resid to 0
            
#             if tgt_resid not in vertex_mapping:
#                 vertex_mapping[tgt_resid] = vertex_index
#                 vertex_index += 1

#             data.append((src_resid, tgt_resid, avg_cpl))
    
#     size = len(vertex_mapping)
#     graph = Graph(size)
#     graph.add_vertex_data(0, 'A', source_resid)  # Add source resid as vertex 0

#     # Add the remaining vertex data
#     for resid, index in vertex_mapping.items():
#         if index != 0:
#             label = chr(ord('A') + index)
#             graph.add_vertex_data(index, label, resid)

#     # Add edges to the graph
#     for src_resid, tgt_resid, avg_cpl in data:
#         u = vertex_mapping[src_resid]
#         v = vertex_mapping[tgt_resid]
#         graph.add_edge(v, u, avg_cpl)

#     return graph

# def append_to_output_file(output_file, source_resid, final_target_resid, avg_coupling, n):
#     if not os.path.exists(output_file):
#         with open(output_file, 'w') as f:
#             f.write("Source_resid\tFinal_target_resid\tAverage_cpl\tRank\n")
    
#     with open(output_file, 'a') as f:
#         f.write(f"{source_resid}\t{final_target_resid}\t{avg_coupling}\t{n}\n")

# def main():
#     file_path = 'average_coupling_values.txt'
#     output_file = 'output_path.txt'
#     graph = read_average_coupling_values(file_path)

#     if len(sys.argv) < 2:
#         raise ValueError("Rank argument is required.")
    
#     n = int(sys.argv[1])

#     try:
#         # Find the resid of the path with the nth largest coupling
#         nth_highest_path, nth_highest_value = graph.find_nth_highest_coupling_path('A', n)
#         final_target_resid = nth_highest_path[-1]  # Get the final target resid from the path
#         print("Source Resid:", graph.resid_mapping['A'])
#         print("Final Target Resid:", final_target_resid)
#         print("Average Coupling Value:", nth_highest_value)
        
#         # Append the path and its average coupling value to the output file
#         append_to_output_file(output_file, graph.resid_mapping['A'], final_target_resid, nth_highest_value, n)
#     except ValueError as e:
#         print(e)
    
#     # Remove the input file
#     os.remove(file_path)

# if __name__ == "__main__":
#     main()
