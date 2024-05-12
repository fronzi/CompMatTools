#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 19:28:17 2023

@author: marco
"""

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher

# Define the elements
# elements = ["Si", "C", "Li"]
elements = ["Si","Si", "Li"]
unique_structures=[]
structures=[]


# Define the range of lattice parameters
a_range = [3, 4]
b_range = [3, 4]
c_range = [3, 4]
alpha_range = [90, 120]
beta_range = [90, 120]
gamma_range = [90, 120]

# Generate all possible structures
structures = []
for a in a_range:
    for b in b_range:
        for c in c_range:
            for alpha in alpha_range:
                for beta in beta_range:
                    for gamma in gamma_range:
                        lattice = [[a, 0, 0], [0, b, 0], [0, 0, c]]
                        lattice_parameters = [a, b, c, alpha, beta, gamma]
                        for i in range(len(elements)):
                            for j in range(len(elements)):
                                if j >= i:
                                    structure = Structure(lattice, [elements[i], elements[j]], [[0, 0, 0], [0.5, 0.5, 0.5]])
                                    structures.append(structure)

# Remove duplicate structures
matcher = StructureMatcher()
unique_structures = matcher.group_structures(structures)
i=0
# Analyze the symmetry of each structure and print the space group symbol
for structure in unique_structures:
    print(i)
    i=i+1
    analyzer = SpacegroupAnalyzer(structure)
    print(analyzer.get_space_group_symbol())
    
    
# Analyze the symmetry of each structure and print the space group symbol
for i, structure in enumerate(structures):
    analyzer = SpacegroupAnalyzer(structure)
    space_group = analyzer.get_space_group_symbol()
    formula = structure.formula.replace(" ", "")
    poscar_filename = f"{formula}_{space_group}_cell{i+1}.vasp"
    structure.to(filename=poscar_filename, fmt="poscar")
    print(f"POSCAR file for {formula} with {space_group} symmetry saved as {poscar_filename}")
