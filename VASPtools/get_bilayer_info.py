#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 15:03:41 2022

@author: marco
"""

from pymatgen.core.sites import PeriodicSite
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.cif import CifParser
from pymatgen.core.structure import Molecule, Structure

import pymatgen.core.structure

from sklearn.cluster import KMeans
import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
km = KMeans(n_clusters=2, random_state=None)

structure = pymatgen.core.Structure.from_file('CONTCAR')

filename = 'CONTCAR'

cut = 0.9 

lines = open(filename).readlines()
c_axis = lines[4].split()
lattice_parameter = lines[1].split()
split_coords = [line.split() for line in lines[8:8+structure.num_sites]]
z_coords = list()
for coord in split_coords:
        z_coord = float(coord[2]) 
        if z_coord > cut:
            z_coord -= 1
        z_coords.append(z_coord)
max_height = max([z_height for z_height in z_coords])
min_height = min([z_height for z_height in z_coords])
spacing = ((1.0 + min_height) - max_height) * abs(float(c_axis[2]))\
	* float(lattice_parameter[0])

slab = abs(float(c_axis[2]))* float(lattice_parameter[0]) - spacing 

KMeans(n_clusters=2)
km.fit(np.array(z_coords).reshape(-1,1))
dists = euclidean_distances(km.cluster_centers_)

InterlayerDistance = abs(float(c_axis[2]))* float(lattice_parameter[0])*abs(dists[1]-dists[0])

InterlayerDistance[0]

print('Vacuum',spacing)
print('Double layer thickness',slab)
print('Interlayer distance',InterlayerDistance[0])
