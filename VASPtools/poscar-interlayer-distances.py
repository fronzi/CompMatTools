import numpy as np
from ase.io import read
from ase.geometry import get_distances

def read_poscar(filename):
    return read(filename, format='vasp')

def classify_distances(atoms, layer_tolerance=0.1, distance_tolerance=0.05):
    # Determine unique z-coordinates to identify layers
    z_coords = np.unique(atoms.positions[:, 2].round(decimals=10))
    layers = {}
    for z in z_coords:
        layer_indices = [i for i, atom in enumerate(atoms) if np.isclose(atom.position[2], z, atol=layer_tolerance)]
        layers[z] = layer_indices

    # Calculate and classify distances within each layer
    for z, indices in layers.items():
        print(f"\nDistances within the layer at z = {z:.3f}:")
        distances = []
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                idx_i, idx_j = indices[i], indices[j]
                element_i = atoms[idx_i].symbol
                element_j = atoms[idx_j].symbol
                distance = atoms.get_distance(idx_i, idx_j, mic=True)
                distances.append((f"{element_i} (atom {idx_i})", f"{element_j} (atom {idx_j})", distance))

        # Sort distances and classify them with tolerance
        distances_sorted = sorted(distances, key=lambda x: x[2])
        classified_distances = []
        rank = 1

        if distances_sorted:
            current_distance_group = [distances_sorted[0]]
            for current in distances_sorted[1:]:
                if np.isclose(current[2], current_distance_group[-1][2], atol=distance_tolerance):
                    current_distance_group.append(current)
                else:
                    classified_distances.append((rank, current_distance_group))
                    rank += 1
                    current_distance_group = [current]

            # Append the last group
            classified_distances.append((rank, current_distance_group))

        # Print results
        for rank, group in classified_distances:
            for item in group:
                print(f"{rank}th shortest: Distance between {item[0]} and {item[1]}: {item[2]:.3f} Å")

# Load the POSCAR file
atoms = read_poscar('POSCAR')  # Replace 'POSCAR' with your actual file path

# Print and classify the distances between atoms in the same layer with tolerance
classify_distances(atoms)
