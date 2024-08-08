#simplified version of fallback mechanism code with debugging print statements and same initial vctor and ransom threshold throughout
import sys
import os
import numpy as np
import subprocess
import random
from dj_algo import read_average_coupling_values, append_to_output_file

def calculate_vector(start_coords, end_coords):
    return end_coords - start_coords

def calculate_dot_product(v1, v2):
    return np.dot(v1, v2)

def read_gro_file(gro_file):
    if not os.path.exists(gro_file):
        raise FileNotFoundError(f"GRO file not found: {gro_file}")
    
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    atoms = [(int(line[0:5].strip()), float(line[20:28]), float(line[28:36]), float(line[36:44])) for line in lines[2:-1]]
    box_dimensions = np.array(lines[-1].strip().split(), dtype=float)
    
    return atoms, box_dimensions

def find_valid_molecule(gro_file, source_resid, cutoff_distance, output_file, initial_vector, atoms, end_resid, start_coords, random_threshold):
    subprocess.run(["python3", "extract_molecules_dj.py", gro_file, "extracted.gro", str(source_resid), str(cutoff_distance)])
    
    if not os.path.exists("extracted.gro"):
        return None

    with open("extracted.gro", 'r') as f:
        output_resids = set(int(line.split()[0]) for line in f.readlines()[2:-1] if line.strip())
    
    if not output_resids:
        return None

    subprocess.run(["python3", "weights_source_resid_dj.py", str(source_resid), str(cutoff_distance)])
    
    if not os.path.exists('average_coupling_values.txt'):
        return None
    
    graph = read_average_coupling_values('average_coupling_values.txt')
    sorted_paths = [(path, highest_value) for n in range(1, len(output_resids) + 1) 
                    for path, highest_value in [graph.find_nth_highest_coupling_path('A', n)]]

    for rank, (path, highest_value) in enumerate(sorted_paths, start=1):
        final_target_resid = path[-1]
        target_coords = np.mean([atom[1:] for atom in atoms if atom[0] == final_target_resid], axis=0)
        vector_to_target = calculate_vector(start_coords, target_coords)
        dot_product = calculate_dot_product(initial_vector, vector_to_target)

        print(f"Attempt with cutoff distance {cutoff_distance}:")
        print(f"  Initial vector: {initial_vector}")
        print(f"  Vector to target: {vector_to_target}")
        print(f"  Dot product: {dot_product}")
        print(f"  Random threshold: {random_threshold}")

        if dot_product >= random_threshold:
            append_to_output_file(output_file, source_resid, final_target_resid, highest_value, rank, is_in_direction=True)
            return final_target_resid
    
    return None

def find_closest_residue(atoms, end_resid):
    end_coords = np.mean([atom[1:] for atom in atoms if atom[0] == end_resid], axis=0)
    min_distance = float('inf')
    closest_resid = None

    for resid, x, y, z in atoms:
        if resid == end_resid:
            continue
        coords = np.array([x, y, z])
        distance = np.linalg.norm(end_coords - coords)
        if distance < min_distance:
            min_distance = distance
            closest_resid = resid
    
    return closest_resid

def main():
    if len(sys.argv) != 6:
        print("Usage: python3 control_movement_dj.py <gro_file> <output_file> <source_resid> <end_resid> <cutoff_distance>")
        sys.exit(1)
    
    gro_file, output_file = sys.argv[1], sys.argv[2]
    source_resid, end_resid = int(sys.argv[3]), int(sys.argv[4])
    cutoff_distance = float(sys.argv[5])
    
    atoms, _ = read_gro_file(gro_file)
    start_coords = np.mean([atom[1:] for atom in atoms if atom[0] == source_resid], axis=0)
    end_coords = np.mean([atom[1:] for atom in atoms if atom[0] == end_resid], axis=0)
    initial_vector = calculate_vector(start_coords, end_coords)

    # Calculate random threshold once
    random_threshold = random.uniform(0, np.linalg.norm(initial_vector))

    visited_residues = set()
    visited_residues.add(source_resid)
    
    attempts = 0
    adjustment_value = 0.5  # Adjust cutoff distance by 0.5 units

    while True:
        valid_molecule = find_valid_molecule(gro_file, source_resid, cutoff_distance, output_file, initial_vector, atoms, end_resid, start_coords, random_threshold)
        
        if valid_molecule:
            if valid_molecule == end_resid:
                print(f"End residue {end_resid} found. Terminating script.")
                os.remove('average_coupling_values.txt')
                sys.exit(0)
            
            if valid_molecule in visited_residues:
                print(f"Residue {valid_molecule} has already been visited. Skipping.")
                valid_molecule = None
            else:
                visited_residues.add(valid_molecule)
                source_resid = valid_molecule
                start_coords = np.mean([atom[1:] for atom in atoms if atom[0] == source_resid], axis=0)
                attempts = 0  # Reset attempts if successful
        else:
            if attempts == 0:
                cutoff_distance += adjustment_value  # Increase the cutoff distance
                print(f"Increasing cutoff distance to {cutoff_distance} and retrying.")
            elif attempts == 1:
                cutoff_distance -= 2 * adjustment_value  # Decrease the cutoff distance by 1 unit
                print(f"Decreasing cutoff distance to {cutoff_distance} and retrying.")
            else:
                print("Max attempts reached. Fallback to finding closest residue.")
                closest_resid = find_closest_residue(atoms, end_resid)
                if closest_resid is not None:
                    if closest_resid in visited_residues:
                        print(f"Closest residue {closest_resid} has already been visited. Skipping.")
                        closest_resid = None
                    else:
                        source_resid = closest_resid
                        start_coords = np.mean([atom[1:] for atom in atoms if atom[0] == source_resid], axis=0)
                        visited_residues.add(closest_resid)
                        attempts = 0  # Reset attempts after fallback
                else:
                    print("No valid residues found. Terminating.")
                    sys.exit(1)
            
            attempts += 1
        
        os.remove('average_coupling_values.txt')

if __name__ == "__main__":
    main()
#.............................................................................................................................................
#code with fallback mechanism
# import sys
# import os
# import numpy as np
# import subprocess
# import random
# from dj_algo import read_average_coupling_values, append_to_output_file

# def calculate_vector(start_coords, end_coords):
#     return end_coords - start_coords

# def calculate_dot_product(v1, v2):
#     return np.dot(v1, v2)

# def read_gro_file(gro_file):
#     atoms = []
#     box_dimensions = []
#     if not os.path.exists(gro_file):
#         raise FileNotFoundError(f"GRO file not found: {gro_file}")
    
#     with open(gro_file, 'r') as f:
#         lines = f.readlines()
    
#     box_dimensions = np.array(lines[-1].strip().split(), dtype=float)
#     for line in lines[2:-1]:  # Skip header, number of atoms, and box dimensions lines
#         resid = int(line[0:5].strip())
#         x, y, z = map(float, (line[20:28], line[28:36], line[36:44]))
#         atoms.append((resid, x, y, z))
#     return atoms, box_dimensions

# def find_valid_molecule(gro_file, source_resid, cutoff_distance, output_file, initial_vector, atoms, end_resid, start_coords):
#     subprocess.run(["python3", "extract_molecules_dj.py", gro_file, "extracted.gro", str(source_resid), str(cutoff_distance)])
    
#     if not os.path.exists("extracted.gro"):
#         return None

#     with open("extracted.gro", 'r') as f:
#         lines = f.readlines()
#         output_resids = set(int(line.split()[0]) for line in lines[2:-1] if line.strip())
    
#     if not output_resids:
#         return None

#     subprocess.run(["python3", "weights_source_resid_dj.py", str(source_resid), str(cutoff_distance)])
    
#     if not os.path.exists('average_coupling_values.txt'):
#         return None
    
#     graph = read_average_coupling_values('average_coupling_values.txt')
#     sorted_paths = [(path, highest_value) for path, highest_value in 
#                     (graph.find_nth_highest_coupling_path('A', n) for n in range(1, len(output_resids) + 1))]

#     for rank, (path, highest_value) in enumerate(sorted_paths, start=1):
#         final_target_resid = path[-1]
#         target_coords = np.mean([atom[1:] for atom in atoms if atom[0] == final_target_resid], axis=0)
#         vector_to_target = calculate_vector(start_coords, target_coords)
#         dot_product = calculate_dot_product(initial_vector, vector_to_target)
#         random_threshold = random.uniform(0, np.linalg.norm(initial_vector))

#         print(f"Initial vector: {initial_vector}")
#         print(f"Vector to target: {vector_to_target}")
#         print(f"Dot product: {dot_product}")
#         print(f"Random threshold: {random_threshold}")

#         if dot_product >= random_threshold:
#             append_to_output_file(output_file, source_resid, final_target_resid, highest_value, rank, is_in_direction=True)
#             return final_target_resid
    
#     return None

# def find_closest_residue(atoms, end_resid):
#     end_coords = np.mean([atom[1:] for atom in atoms if atom[0] == end_resid], axis=0)
#     min_distance = float('inf')
#     closest_resid = None

#     for resid, x, y, z in atoms:
#         if resid == end_resid:
#             continue
#         coords = np.array([x, y, z])
#         distance = np.linalg.norm(end_coords - coords)
#         if distance < min_distance:
#             min_distance = distance
#             closest_resid = resid
    
#     return closest_resid

# def main():
#     if len(sys.argv) != 6:
#         print("Usage: python3 control_movement_dj.py <gro_file> <output_file> <source_resid> <end_resid> <cutoff_distance>")
#         sys.exit(1)
    
#     gro_file = sys.argv[1]
#     output_file = sys.argv[2]
#     source_resid = int(sys.argv[3])
#     end_resid = int(sys.argv[4])
#     cutoff_distance = float(sys.argv[5])
    
#     atoms, box_dimensions = read_gro_file(gro_file)

#     start_coords = np.mean([atom[1:] for atom in atoms if atom[0] == source_resid], axis=0)
#     end_coords = np.mean([atom[1:] for atom in atoms if atom[0] == end_resid], axis=0)
#     initial_vector = calculate_vector(start_coords, end_coords)

#     attempts = 0
#     max_attempts = 2  # Two attempts: one with increased and one with decreased cutoff distance
#     adjustment_value = 0.5  # The value to adjust the cutoff distance by

#     while True:
#         valid_molecule = find_valid_molecule(gro_file, source_resid, cutoff_distance, output_file, initial_vector, atoms, end_resid, start_coords)
        
#         if valid_molecule:
#             if valid_molecule == end_resid:
#                 print(f"End residue {end_resid} found. Terminating script.")
#                 os.remove('average_coupling_values.txt')
#                 sys.exit(0)
            
#             source_resid = valid_molecule
#             start_coords = np.mean([atom[1:] for atom in atoms if atom[0] == source_resid], axis=0)
#             initial_vector = calculate_vector(start_coords, end_coords)
#             attempts = 0  # Reset attempts if successful
#         else:
#             if attempts == 0:
#                 cutoff_distance += adjustment_value  # Increase the cutoff distance
#                 print(f"Increasing cutoff distance to {cutoff_distance} and retrying.")
#                 attempts += 1
#             elif attempts == 1:
#                 cutoff_distance -= 2 * adjustment_value  # Decrease the cutoff distance by 1 unit
#                 print(f"Decreasing cutoff distance to {cutoff_distance} and retrying.")
#                 attempts += 1
#             else:
#                 print("Max attempts reached. Fallback to finding closest residue.")
#                 closest_resid = find_closest_residue(atoms, end_resid)
#                 if closest_resid is not None:
#                     source_resid = closest_resid
#                     start_coords = np.mean([atom[1:] for atom in atoms if atom[0] == source_resid], axis=0)
#                     initial_vector = calculate_vector(start_coords, end_coords)
#                     attempts = 0  # Reset attempts after fallback
#                 else:
#                     print("No valid residues found. Terminating.")
#                     sys.exit(1)
        
#         os.remove('average_coupling_values.txt')

# if __name__ == "__main__":
#     main()
#...............................................................................................................................................
# #simplified code with unique resids
# import sys
# import os
# import numpy as np
# import pandas as pd
# import subprocess
# import random
# from dj_algo import read_average_coupling_values, append_to_output_file

# def calculate_vector(start_coords, end_coords):
#     return end_coords - start_coords

# def calculate_dot_product(v1, v2):
#     return np.dot(v1, v2)

# def read_gro_file(gro_file):
#     if not os.path.exists(gro_file):
#         raise FileNotFoundError(f"GRO file not found: {gro_file}")
    
#     with open(gro_file, 'r') as f:
#         lines = f.readlines()
    
#     box_dimensions = np.array(lines[-1].strip().split(), dtype=float)
#     data = [(line[0:5].strip(), line[20:28], line[28:36], line[36:44]) for line in lines[2:-1]]
#     df = pd.DataFrame(data, columns=['resid', 'x', 'y', 'z'], dtype=float)
#     df['resid'] = df['resid'].astype(int)
#     return df, box_dimensions

# def main():
#     if len(sys.argv) != 6:
#         print("Usage: python3 control_movement_dj.py <gro_file> <output_file> <source_resid> <end_resid> <cutoff_distance>")
#         sys.exit(1)
    
#     gro_file, output_file, source_resid, end_resid, cutoff_distance = sys.argv[1:6]
#     source_resid, end_resid, cutoff_distance = int(source_resid), int(end_resid), float(cutoff_distance)
    
#     df, box_dimensions = read_gro_file(gro_file)

#     start_coords = df[df['resid'] == source_resid][['x', 'y', 'z']].mean().to_numpy()
#     end_coords = df[df['resid'] == end_resid][['x', 'y', 'z']].mean().to_numpy()
#     initial_vector = calculate_vector(start_coords, end_coords)
    
#     visited_resids = set([source_resid])

#     while True:
#         subprocess.run(["python3", "extract_molecules_dj.py", gro_file, "extracted.gro", str(source_resid), str(cutoff_distance)])
        
#         if not os.path.exists("extracted.gro"):
#             print("Extracted file not found. Terminating.")
#             break

#         extracted_df = pd.read_csv("extracted.gro", delim_whitespace=True, skiprows=2, header=None)
#         output_resids = set(extracted_df[0].astype(int))

#         output_resids -= visited_resids

#         if not output_resids:
#             print("No new molecules found around the source residue. Terminating.")
#             break

#         subprocess.run(["python3", "weights_source_resid_dj.py", str(source_resid), str(cutoff_distance)])
        
#         if not os.path.exists('average_coupling_values.txt'):
#             print("Average coupling values file not found. Terminating.")
#             break
        
#         graph = read_average_coupling_values('average_coupling_values.txt')
#         sorted_paths = [(path, highest_value) for path, highest_value in 
#                         (graph.find_nth_highest_coupling_path('A', n) for n in range(1, len(output_resids) + 1))]

#         found_valid_molecule = False
        
#         for rank, (path, highest_value) in enumerate(sorted_paths, start=1):
#             final_target_resid = path[-1]
            
#             if final_target_resid in visited_resids:
#                 continue

#             target_coords = df[df['resid'] == final_target_resid][['x', 'y', 'z']].mean().to_numpy()
#             vector_to_target = calculate_vector(start_coords, target_coords)
#             dot_product = calculate_dot_product(initial_vector, vector_to_target)
#             random_threshold = random.uniform(0, np.linalg.norm(initial_vector))

#             print(f"Initial vector: {initial_vector}")
#             print(f"Vector to target: {vector_to_target}")
#             print(f"Dot product: {dot_product}")
#             print(f"Random threshold: {random_threshold}")

#             if dot_product >= random_threshold:
#                 print(f"Valid molecule found: {final_target_resid} with coupling value {highest_value}")
#                 is_in_direction = True
#                 append_to_output_file(output_file, source_resid, final_target_resid, highest_value, rank, is_in_direction)
#                 found_valid_molecule = True
#                 visited_resids.add(final_target_resid)
                
#                 if final_target_resid == end_resid:
#                     print(f"End residue {end_resid} found. Terminating script.")
#                     os.remove('average_coupling_values.txt')
#                     sys.exit(0)

#                 break
        
#         if not found_valid_molecule:
#             print("No valid molecules found that satisfy the condition. Terminating.")
#             break

#         source_resid = final_target_resid
#         start_coords = target_coords
#         initial_vector = calculate_vector(start_coords, end_coords)

#         os.remove('average_coupling_values.txt')

# if __name__ == "__main__":
#     main()


# #code to terminate when it finds the end_resid
# import sys
# import os
# import numpy as np
# import subprocess
# import random
# from dj_algo import read_average_coupling_values, append_to_output_file

# def calculate_vector(start_coords, end_coords):
#     return end_coords - start_coords

# def calculate_dot_product(v1, v2):
#     return np.dot(v1, v2)

# def read_gro_file(gro_file):
#     atoms = []
#     box_dimensions = []
#     if not os.path.exists(gro_file):
#         raise FileNotFoundError(f"GRO file not found: {gro_file}")
    
#     with open(gro_file, 'r') as f:
#         lines = f.readlines()
    
#     box_dimensions = np.array(lines[-1].strip().split(), dtype=float)
#     for line in lines[2:-1]:  # Skip header, number of atoms, and box dimensions lines
#         resid = int(line[0:5].strip())
#         x, y, z = map(float, (line[20:28], line[28:36], line[36:44]))
#         atoms.append((resid, x, y, z))
#     return atoms, box_dimensions

# def main():
#     if len(sys.argv) != 6:
#         print("Usage: python3 control_movement_dj.py <gro_file> <output_file> <source_resid> <end_resid> <cutoff_distance>")
#         sys.exit(1)
    
#     gro_file = sys.argv[1]
#     output_file = sys.argv[2]
#     source_resid = int(sys.argv[3])
#     end_resid = int(sys.argv[4])
#     cutoff_distance = float(sys.argv[5])
    
#     atoms, box_dimensions = read_gro_file(gro_file)

#     start_coords = np.mean([atom[1:] for atom in atoms if atom[0] == source_resid], axis=0)
#     end_coords = np.mean([atom[1:] for atom in atoms if atom[0] == end_resid], axis=0)
#     initial_vector = calculate_vector(start_coords, end_coords)

#     while True:
#         subprocess.run(["python3", "extract_molecules_dj.py", gro_file, "extracted.gro", str(source_resid), str(cutoff_distance)])
        
#         if not os.path.exists("extracted.gro"):
#             print("Extracted file not found. Terminating.")
#             break

#         with open("extracted.gro", 'r') as f:
#             lines = f.readlines()
#             output_resids = set(int(line.split()[0]) for line in lines[2:-1] if line.strip())
        
#         if not output_resids:
#             print("No molecules found around the source residue. Terminating.")
#             break

#         subprocess.run(["python3", "weights_source_resid_dj.py", str(source_resid), str(cutoff_distance)])
        
#         if not os.path.exists('average_coupling_values.txt'):
#             print("Average coupling values file not found. Terminating.")
#             break
        
#         graph = read_average_coupling_values('average_coupling_values.txt')
#         sorted_paths = [(path, highest_value) for path, highest_value in 
#                         (graph.find_nth_highest_coupling_path('A', n) for n in range(1, len(output_resids) + 1))]

#         found_valid_molecule = False
        
#         for rank, (path, highest_value) in enumerate(sorted_paths, start=1):
#             final_target_resid = path[-1]
#             target_coords = np.mean([atom[1:] for atom in atoms if atom[0] == final_target_resid], axis=0)
#             vector_to_target = calculate_vector(start_coords, target_coords)
#             dot_product = calculate_dot_product(initial_vector, vector_to_target)
#             random_threshold = random.uniform(0, np.linalg.norm(initial_vector))

#             print(f"Initial vector: {initial_vector}")
#             print(f"Vector to target: {vector_to_target}")
#             print(f"Dot product: {dot_product}")
#             print(f"Random threshold: {random_threshold}")

#             if dot_product >= random_threshold:
#                 print(f"Valid molecule found: {final_target_resid} with coupling value {highest_value}")
#                 is_in_direction = True
#                 append_to_output_file(output_file, source_resid, final_target_resid, highest_value, rank, is_in_direction)
#                 found_valid_molecule = True
                
#                 # Check if the final target residue is the end residue
#                 if final_target_resid == end_resid:
#                     print(f"End residue {end_resid} found. Terminating script.")
#                     os.remove('average_coupling_values.txt')
#                     sys.exit(0)  # Terminate the script successfully

#                 break
        
#         if not found_valid_molecule:
#             print("No valid molecules found that satisfy the condition. Terminating.")
#             # Optionally clear the output file or add any other cleanup code here
#             break

#         source_resid = final_target_resid
#         start_coords = target_coords
#         initial_vector = calculate_vector(start_coords, end_coords)

#         os.remove('average_coupling_values.txt')

# if __name__ == "__main__":
#     main()

#.......................................................................................................................................
# #code for running over nth highest coupling values
# import sys
# import os
# import numpy as np
# import subprocess
# import random
# from dj_algo import read_average_coupling_values, append_to_output_file

# def calculate_vector(start_coords, end_coords):
#     return end_coords - start_coords

# def calculate_dot_product(v1, v2):
#     return np.dot(v1, v2)

# def read_gro_file(gro_file):
#     atoms = []
#     box_dimensions = []
#     if not os.path.exists(gro_file):
#         raise FileNotFoundError(f"GRO file not found: {gro_file}")
    
#     with open(gro_file, 'r') as f:
#         lines = f.readlines()
    
#     box_dimensions = np.array(lines[-1].strip().split(), dtype=float)
#     for line in lines[2:-1]:  # Skip header, number of atoms, and box dimensions lines
#         resid = int(line[0:5].strip())
#         x, y, z = map(float, (line[20:28], line[28:36], line[36:44]))
#         atoms.append((resid, x, y, z))
#     return atoms, box_dimensions

# def main():
#     if len(sys.argv) != 6:
#         print("Usage: python3 control_movement_dj.py <gro_file> <output_file> <source_resid> <end_resid> <cutoff_distance>")
#         sys.exit(1)
    
#     gro_file = sys.argv[1]
#     output_file = sys.argv[2]
#     source_resid = int(sys.argv[3])
#     end_resid = int(sys.argv[4])
#     cutoff_distance = float(sys.argv[5])
    
#     atoms, box_dimensions = read_gro_file(gro_file)

#     start_coords = np.mean([atom[1:] for atom in atoms if atom[0] == source_resid], axis=0)
#     end_coords = np.mean([atom[1:] for atom in atoms if atom[0] == end_resid], axis=0)
#     initial_vector = calculate_vector(start_coords, end_coords)

#     while True:
#         subprocess.run(["python3", "extract_molecules_dj.py", gro_file, "extracted.gro", str(source_resid), str(cutoff_distance)])
        
#         if not os.path.exists("extracted.gro"):
#             print("Extracted file not found. Terminating.")
#             break

#         with open("extracted.gro", 'r') as f:
#             lines = f.readlines()
#             output_resids = set(int(line.split()[0]) for line in lines[2:-1] if line.strip())
        
#         if not output_resids:
#             print("No molecules found around the source residue. Terminating.")
#             break

#         subprocess.run(["python3", "weights_source_resid_dj.py", str(source_resid), str(cutoff_distance)])
        
#         if not os.path.exists('average_coupling_values.txt'):
#             print("Average coupling values file not found. Terminating.")
#             break
        
#         graph = read_average_coupling_values('average_coupling_values.txt')
#         sorted_paths = [(path, highest_value) for path, highest_value in 
#                         (graph.find_nth_highest_coupling_path('A', n) for n in range(1, len(output_resids) + 1))]

#         found_valid_molecule = False
        
#         for rank, (path, highest_value) in enumerate(sorted_paths, start=1):
#             final_target_resid = path[-1]
#             target_coords = np.mean([atom[1:] for atom in atoms if atom[0] == final_target_resid], axis=0)
#             vector_to_target = calculate_vector(start_coords, target_coords)
#             dot_product = calculate_dot_product(initial_vector, vector_to_target)
#             random_threshold = random.uniform(0, np.linalg.norm(initial_vector))

#             print(f"Initial vector: {initial_vector}")
#             print(f"Vector to target: {vector_to_target}")
#             print(f"Dot product: {dot_product}")
#             print(f"Random threshold: {random_threshold}")

#             if dot_product >= random_threshold:
#                 print(f"Valid molecule found: {final_target_resid} with coupling value {highest_value}")
#                 is_in_direction = True
#                 append_to_output_file(output_file, source_resid, final_target_resid, highest_value, rank, is_in_direction)
#                 found_valid_molecule = True
#                 break
        
#         if not found_valid_molecule:
#             print("No valid molecules found that satisfy the condition. Terminating.")
#             # Optionally clear the output file or add any other cleanup code here
#             break

#         source_resid = final_target_resid
#         start_coords = target_coords
#         initial_vector = calculate_vector(start_coords, end_coords)

#         os.remove('average_coupling_values.txt')

# if __name__ == "__main__":
#     main()

#......................................................................................................................................

#code just to check whether it is moving in correct vector direction
# import sys
# import os
# import numpy as np
# import subprocess
# import random
# from dj_algo import read_average_coupling_values, append_to_output_file

# def calculate_vector(start_coords, end_coords):
#     return end_coords - start_coords

# def calculate_dot_product(v1, v2):
#     return np.dot(v1, v2)

# def read_gro_file(gro_file):
#     atoms = []
#     box_dimensions = []
#     if not os.path.exists(gro_file):
#         raise FileNotFoundError(f"GRO file not found: {gro_file}")
    
#     with open(gro_file, 'r') as f:
#         lines = f.readlines()
    
#     box_dimensions = np.array(lines[-1].strip().split(), dtype=float)
#     for line in lines[2:-1]:  # Skip header, number of atoms, and box dimensions lines
#         resid = int(line[0:5].strip())
#         x, y, z = map(float, (line[20:28], line[28:36], line[36:44]))
#         atoms.append((resid, x, y, z))
#     return atoms, box_dimensions

# def main():
#     if len(sys.argv) != 6:
#         print("Usage: python3 control_movement_dj.py <gro_file> <output_file> <source_resid> <end_resid> <cutoff_distance>")
#         sys.exit(1)
    
#     gro_file = sys.argv[1]
#     output_file = sys.argv[2]
#     source_resid = int(sys.argv[3])
#     end_resid = int(sys.argv[4])
#     cutoff_distance = float(sys.argv[5])
    
#     atoms, box_dimensions = read_gro_file(gro_file)

#     start_coords = np.mean([atom[1:] for atom in atoms if atom[0] == source_resid], axis=0)
#     end_coords = np.mean([atom[1:] for atom in atoms if atom[0] == end_resid], axis=0)
#     initial_vector = calculate_vector(start_coords, end_coords)
    
#     while True:
#         subprocess.run(["python3", "extract_molecules_dj.py", gro_file, "extracted.gro", str(source_resid), str(cutoff_distance)])
        
#         if not os.path.exists("extracted.gro"):
#             print("Extracted file not found. Terminating.")
#             break

#         with open("extracted.gro", 'r') as f:
#             lines = f.readlines()
#             output_resids = set(int(line.split()[0]) for line in lines[2:-1] if line.strip())
        
#         if not output_resids:
#             print("No molecules found around the source residue. Terminating.")
#             break

#         subprocess.run(["python3", "weights_source_resid_dj.py", str(source_resid), str(cutoff_distance)])
        
#         if not os.path.exists('average_coupling_values.txt'):
#             print("Average coupling values file not found. Terminating.")
#             break
        
#         graph = read_average_coupling_values('average_coupling_values.txt')
#         path, highest_value = graph.find_nth_highest_coupling_path('A', 1)
#         final_target_resid = path[-1]

#         target_coords = np.mean([atom[1:] for atom in atoms if atom[0] == final_target_resid], axis=0)
#         vector_to_target = calculate_vector(start_coords, target_coords)

#         dot_product = calculate_dot_product(initial_vector, vector_to_target)
#         random_threshold = random.uniform(0, np.linalg.norm(initial_vector))

#         print(f"Initial vector: {initial_vector}")
#         print(f"Vector to target: {vector_to_target}")
#         print(f"Dot product: {dot_product}")
#         print(f"Random threshold: {random_threshold}")

#         if dot_product < random_threshold:
#             print(f"Terminating: The molecule {final_target_resid} is not in the positive x direction.")
#             in_reference_direction = False
#             append_to_output_file(output_file, source_resid, final_target_resid, highest_value, in_reference_direction)
#             break

#         in_reference_direction = True
#         append_to_output_file(output_file, source_resid, final_target_resid, highest_value, in_reference_direction)

#         source_resid = final_target_resid
#         start_coords = target_coords
#         initial_vector = calculate_vector(start_coords, end_coords)

#         os.remove('average_coupling_values.txt')

# if __name__ == "__main__":
#     main()
