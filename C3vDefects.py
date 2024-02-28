from itertools import permutations
from ase.io import read, write
from ase.neighborlist import NeighborList
import numpy as np

def substitute_and_write(diamond, sub_elements, filename):
      """
    Substitutes specified atoms in the diamond structure and writes to a POSCAR file in a C3v symmetry.
    
    Parameters:
    diamond (Atoms): ASE Atoms object representing the diamond structure.
    sub_elements (list): List of elements for substitution. 'V' denotes a vacancy.
    filename (str): Filename for the output POSCAR file.
    """
    # Define cutoff for carbon-carbon bond in diamond (approximately 1.6 Å)
    cutoff = 1.6

    # Create a NeighborList
    neighbor_list = NeighborList([cutoff/2.0] * len(diamond), self_interaction=False, bothways=True)
    neighbor_list.update(diamond)

    # Find the central carbon atom
    central_atom_index = find_central_atom(diamond)

    # Function to mark for substitution or vacancy
    def mark_substitution_or_vacancy(atom_index, element):
        if element == 'V':
            return (atom_index, True)  # Mark for deletion
        else:
            diamond[atom_index].symbol = element  # Substitute atom
            return (atom_index, False)  # Not marked for deletion

    # List to store indices of atoms marked for vacancy
    vacancies = []

    # Mark the central carbon atom for substitution or vacancy
    central_marked = mark_substitution_or_vacancy(central_atom_index, sub_elements[0])
    if central_marked[1]:
        vacancies.append(central_marked[0])

    # Find the neighbors of the central atom
    neighbors = neighbor_list.get_neighbors(central_atom_index)[0]

    # Substitute or mark for vacancy the neighbors
    for i, neighbor_index in enumerate(neighbors):
        if i == 0:
            marked = mark_substitution_or_vacancy(neighbor_index, sub_elements[1])
        else:
            marked = mark_substitution_or_vacancy(neighbor_index, sub_elements[2])
        if marked[1]:
            vacancies.append(marked[0])

    # Handle the carbon neighbors of the first neighbor
    C2_index = neighbors[0] if len(neighbors) > 0 else None
    if C2_index is not None:
        C2_neighbors = neighbor_list.get_neighbors(C2_index)[0]
        for neighbor_index in C2_neighbors:
            if diamond[neighbor_index].symbol == 'C':
                marked = mark_substitution_or_vacancy(neighbor_index, sub_elements[3])
                if marked[1]:
                    vacancies.append(marked[0])

    # Delete atoms marked for vacancies
    for index in sorted(vacancies, reverse=True):  # Reverse sort to maintain correct indices
        del diamond[index]

    # Write the modified structure to a new POSCAR file
    write(filename, diamond, format='vasp', vasp5=True, direct=True)

def find_central_atom(atoms):
    center = np.array(atoms.get_cell().diagonal()) / 2
    distances = [np.linalg.norm(atom.position - center) for atom in atoms]
    return np.argmin(distances)
