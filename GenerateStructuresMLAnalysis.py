#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 07:32:07 2023

@author: marco
"""
import sys
import os
# import itertools
import pandas as pd
from matminer.utils.io import load_dataframe_from_json
# from matminer.featurizers.structure import StructureFeatures
from matminer.featurizers.conversions import StrToComposition
# from matminer.featurizers import MultipleFeaturizer
from matminer.featurizers.composition import ElementProperty
from pymatgen.core import Lattice

from itertools import combinations
from pymatgen.core import Element
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from matplotlib.ticker import StrMethodFormatter

from pymatgen.core.composition import Element, Composition
from pymatgen.core.periodic_table import Specie



from mp_api.client import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PatchedPhaseDiagram, GrandPotentialPhaseDiagram, PDEntry
import matplotlib.pyplot as plt

from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, CompoundPhaseDiagram
# from pymatgen.core.composition import Elements


###############################################################################
###     generate all crystal structures
################################################################################



#  MaterialsProject API 6NN5rK4gahmRJCdlBuhNpGpueNgkOgDd   
path = "/Users/marco/Dropbox/Work/WORKING_DIR_UNIMELB/GraphiteSiBatteries/Data/"

# Initialize MPRester and get the phase diagram
mpr = MPRester("6NN5rK4gahmRJCdlBuhNpGpueNgkOgDd")
# comp = "Si-Li-O"
# comp = "Si-Li-C"
# comp = "Si-O-C"
# comp = "Si-Li-F"
# comp = "Si-Li-F"
# comp = "Si-C-F"
# comp = "Si-O-F"
# comp = "C-Li-O"
# comp = "C-Li-F"
comp = ["Si","Li","C","O"]
# comp = Compositio["Si-Li-C-O"]
comp = Composition({'Si': 1, 'Li': 1})
# comp2 = Elements(comp)
entries = mpr.get_entries_in_chemsys(comp)

formulas = []
space_groups = []
lattice_constants_a = []
lattice_constants_b = []
lattice_constants_c = []
crystal_systems = []
mp_ids=[]

# Generate the POSCAR file for each entry
poscars=[]




for entry in entries:
    mp_id = entry.material_id
    formula = entry.composition.reduced_formula
    # mp_id = mp_id_suffix.split("-")[0] + "-" + mp_id_suffix.split("-")[1]
    structure = entry.structure
    poscar = structure.to(fmt="poscar")
    poscars.append(poscar)
    # Get structure and symmetry data
    mp_ids.append(mp_id)
    structure = entry.structure
    analyzer = SpacegroupAnalyzer(structure)
    space_group = analyzer.get_space_group_symbol()
    crystal_system = analyzer.get_crystal_system()
    lattice_constant_a = structure.lattice.a
    lattice_constant_b = structure.lattice.b
    lattice_constant_c = structure.lattice.c
    

    # Append data to lists
    formulas.append(entry.composition.reduced_formula)
    space_groups.append(space_group)
    lattice_constants_a.append(lattice_constant_a)
    lattice_constants_b.append(lattice_constant_b)
    lattice_constants_c.append(lattice_constant_c)
    crystal_systems.append(crystal_system)

    # Create dataframe with data
    data = {'Formula': formulas,
            'Space Group': space_groups,
            'Lattice Constant a': lattice_constants_a,
            'Lattice Constant b': lattice_constants_b,
            'Lattice Constant c': lattice_constants_c,
            'Crystal System': crystal_systems}

    df = pd.DataFrame(data)
    df_id = pd.DataFrame(mp_ids)

    # Save dataframe as json
    # df.to_json(path+"structure_statistics.json")
    # df = pd.read_json(path+"structure_statistics.json")
    
    
for i, entry in enumerate(entries):
    structure = entry.structure
    poscar = structure.to(fmt="poscar")
    file_path = os.path.join(path, f"poscar_{i}.vasp")
    with open(file_path, "w") as f:
        f.write(poscar)




sys.exit(-1)

df_id.to_csv(path+'mp_ids.csv')


data = []
###creare dataframe for machine learning #####
mpr.materials.search('mp-600071', fields=['formula_pretty', 'compliance_tensor','elastic_tensor', 'elastic_tensor_original'])   

mpr.get_task_ids_associated_with_material_id('mp-600071')
mpr.oxidation_states.search('mp-600071')
mpr.get_structure_by_material_id('mp-600071')

dppp=mpr.materials.search(crystal_system="Cubic", all_fields=False, fields=['material_id', 'formula', 'nsites', 'space_grop','volume','structure','elastic_anisotropy', 'G_Reuss', 'G_VRH', 'G_Voigt','K_Reuss', 'K_VRH', 'K_Voigt', 'poisson_ratio', 'compliance_tensor','elastic_tensor', 'elastic_tensor_original'])   


############################################################

with MPRester("lL4tFjJiSQT8cZ3CGwwBZJJViCpmPzLf") as mpr:
    entries = mpr.get_entries_in_chemsys(["Li", "Si", "C"])



# phase_diagram = GrandPotentialPhaseDiagram(entries)
phase_diagram = PhaseDiagram(entries)
plotter = PDPlotter(phase_diagram, show_unstable=0, markersize=18)
phases=plotter.get_chempot_range_map_plot([Element("C"), Element("Si")])
phases.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # No decimal places
phases.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # 2 decimal places
phases.xticks(fontsize=24)
phases.yticks(fontsize=24)
phases.xlabel('$\mu_C$ (eV)', fontsize=34)
phases.ylabel('$\mu_O$ (eV)', fontsize=34)
# phases.title('Gibbs Free Energy Volume (Li-2atoms-NO-EP) ', fontsize=16)
phases.tight_layout()
# phases.legend(fontsize=4)
phases.grid()
phases.savefig(path+'1+O.png',dpi=600)


stab=phase_diagram.stable_entries

for i, entry in enumerate(stab):
    structure = entry.structure
    poscar = structure.to(fmt="poscar")
    file_path = os.path.join(path+'/POSCARs', f"poscar_{i}.vasp")
    with open(file_path, "w") as f:
        f.write(poscar)




compound_plotter.show()
contour = plotter.get_contour_pd_plot()
plotter = PDPlotter(phase_diagram, show_unstable=0, markersize=20)
cpd = CompoundPhaseDiagram(entries, [Composition("Li"), Composition("SiC")])
phases.savefig(path+'1+O.png',dpi=600)
# df = pd.read_json(path+"structure_statistics.json")



# # Plot the lattice constants distribution
# fig, ax = plt.subplots()
# ax.hist(df['Lattice Constant'],bins=200)
# ax.set_title('Lattice Constants Distribution')
# ax.set_xlabel('Lattice Constant (Angstroms)')
# ax.set_ylabel('Counts')
# ax.tick_params(axis='x', labelsize=18, rotation=0)
# plt.tight_layout()
# plt.savefig(path+'LatticeConst',dpi=600)
# plt.show()

# # Plot the crystal system distribution
# fig, ax = plt.subplots()
# ax.hist(df['Crystal System'])
# ax.set_title('Crystal System Distribution')
# ax.set_xlabel('Crystal System')
# ax.set_ylabel('Counts')
# ax.tick_params(axis='x', labelsize=20, rotation=90)
# plt.tight_layout()
# plt.savefig(path+'CrystalSystem',dpi=600)
# plt.show()

# # Plot the space group distribution
# fig, ax = plt.subplots()
# ax.hist(df['Space Group'], bins=400)
# ax.set_title('Space Group Distribution')
# ax.set_xlabel('Space Group')
# ax.set_ylabel('Counts')
# ax.tick_params(axis='x', labelsize=4, rotation=90)
# plt.tight_layout()
# plt.savefig(path+'SpaceGroup',dpi=600)
# plt.show()


sys.exit(-1)

###############################################################################
###  geberate POSCAR files and diuagram stability DFT 
################################################################################
    

# for i, entry in enumerate(entries):
#     structure = entry.structure
#     poscar = structure.to(fmt="poscar")
#     file_path = os.path.join(path, f"poscar_{i}.vasp")
#     with open(file_path, "w") as f:
#         f.write(poscar)
  



pd = PhaseDiagram(entries)




# Plot the phase diagram
plotter = PDPlotter(pd,show_unstable=False)
stability=plotter.get_contour_pd_plot()
stability.savefig(path+'stability8',dpi=600)
stability.show()

stability2=plotter.get_plot(label_stable=False)
stability2.show()




sys.exit(-1)



# Get the most stable entry and its formula
most_stable_entry = pd.stable_entries
formula = most_stable_entry.composition.reduced_formula

# Get the thermodynamic properties of the most stable entry
data = mpr.get_data(formula)
energy = data[0]['energy']
formation_energy = data[0]['formation_energy_per_atom']

# Plot the stability of the most stable entry on the phase diagram
plt.plot(formation_energy, energy, 'ro')
plt.annotate(formula, (formation_energy, energy))
plt.savefig(path+'Stabilityternary.png', dpi=600)
plt.show()

sys.exit(-1)
###############################################################################
###  Genertate Descriptoirs
################################################################################
    



from matminer.featurizers.structure import SiteStatsFingerprint

from matminer.featurizers.composition import ElementProperty
import json
from pymatgen.io.vasp import Poscar
from matminer.featurizers.structure import StructureComposition, GlobalSymmetryFeatures, SiteStatsFingerprint

from pymatgen.core.structure import Structure
import pandas as pd
import os


# Load the list of POSCAR files from the JSON file
with open(path+"structures.json", "r") as f:
    structures = json.load(f)

# Create empty lists to store the computed descriptors and corresponding structure IDs
descriptors = []
structure_ids = []

structures=[]

# Loop over each POSCAR file and compute the descriptors
for i, poscar in enumerate(structures):
    print(poscar)
    # Load the structure from the POSCAR file
    structures[i] = Poscar.from_string(poscar).structure
    
    # Compute the descriptors
    comp_features = StructureComposition(structures[i])
    # gs_features = GlobalSymmetryFeatures(structures[i])  not working
     
    
    # comp_desc = comp_features.featurize(structure) not working
    # gs_desc = gs_features.featurize(structure)    not working
    # ssf_desc = ssf_features.featurize(structure)  not working


# Loop through all POSCAR files in the folder
all_descriptors = []
for file in os.listdir(path):
    if file.endswith(".vasp"):
        # Read the structure from the POSCAR file
        structure = Structure.from_file(os.path.join(path, file))
        # Generate descriptors for the structure
        struct_descriptors = descriptors.featurize(structure)
        all_descriptors.append(struct_descriptors)

# Convert the list of descriptor arrays to a DataFrame
df = pd.DataFrame(all_descriptors, columns=descriptors.feature_labels())



sys.exit(-1)








###############################################################################
###  ML to be fixed
################################################################################
    
import os
import pandas as pd
from matminer.featurizers.structure import SiteStatsFingerprint
from matminer.utils.io import load_dataframe_from_json
from sklearn.externals import joblib
import joblib

# Define the path to the directory containing the POSCAR files


# Define the featurizer to use
featurizer = SiteStatsFingerprint.from_preset("CoordinationNumber_ward-prb-2017")

# Load the pre-trained models
# elasticity_model = joblib.load("elasticity_model.joblib")
formation_energy_model = joblib.load("formation_energy_model.joblib")

# Define a function to predict the stability of a single structure
def predict_stability(path):
    # Load the structure from the POSCAR file
    structure = Structure.from_file(poscar_path)

    # Compute features for the structure
    features = featurizer.featurize(structure)

    # Convert features to a pandas DataFrame
    df = pd.DataFrame(features).transpose()

    # Use the models to predict the stability of the structure
    # elasticity = elasticity_model.predict(df)[0]
    formation_energy = formation_energy_model.predict(df)[0]

    # return {"elasticity": elasticity, "formation_energy": formation_energy}
    return { "formation_energy": formation_energy}

# Loop over all POSCAR files in the directory and predict their stability
results = {}
for filename in os.listdir(path):
    if filename.endswith(".POSCAR"):
        poscar_path = os.path.join(path, filename)
        results[filename] = predict_stability(poscar_path)

# Print the results
for filename, result in results.items():
    print(filename, result)



################################################################################
################################################################################
################################################################################


