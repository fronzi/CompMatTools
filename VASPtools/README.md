#  Interlayer Atomic Distance Classifier README

## Overview




This repository contains Python scripts designed for analyzing and processing electronic structure data from computational materials science simulations using VASP outputs.


### 0. oscar-interlayer-distances.py
This script is designed to analyze POSCAR files using the Atomic Simulation Environment (ASE) to determine and classify distances between atoms within the same z-layers of a crystal structure.
It groups these distances based on a specified tolerance, allowing for detailed examination of atomic arrangements and potentially identifying significant structural features.



### 1. BandGap_collect_save.py
**Function:** Extracts band gap information from multiple `vasprun.xml` files located within a specified directory, collating results into a CSV file.
**Input:** Directory path containing `vasprun.xml` files.
**Output:** `band_gaps.csv` containing columns for directory, band gap energy, and whether the band gap is direct.

### 2. C3vDefects.py
**Function:** Modifies a diamond structure by substituting specified atoms to create defects, outputting the result as a VASP POSCAR file.
**Input:** POSCAR file of the original diamond structure, list of elements for substitution.
**Output:** Modified POSCAR file reflecting the introduced defects.

### 3. ChargeCount.py
**Function:** Calculates the total electronic charge from a `CHGCAR` file and analyzes charge distribution across specified interfaces.
**Input:** `CHGCAR` file.
**Output:** Prints total charge and charge distribution details to standard output.

### 4. DiffusionCoeff.py
**Function:** Calculates the diffusion coefficient from velocity autocorrelation function data files using the Green-Kubo relation.
**Input:** Directory containing data files starting with "VACF".
**Output:** `diffusion_coefficients.csv` with diffusion coefficients for each material.

### 5. ElasticProperties.py
**Function:** Reads an elastic tensor from an OUTCAR file and calculates various mechanical properties such as Young's modulus and shear modulus.
**Input:** OUTCAR file.
**Output:** `elastic_properties.csv` detailing the computed elastic properties.

### 6. GenerateStructure.py
**Function:** Generates and writes crystal structures to POSCAR files based on specified lattice parameters and symmetry.
**Input:** Lattice parameters (a, b, c, alpha, beta, gamma), directory to save output.
**Output:** POSCAR files named according to lattice parameters and symmetry in the specified directory.

### 7. GenerateStructuresMLAnalysis.py
**Function:** Fetches structures from the Materials Project API and computes descriptors for machine learning analysis.
**Input:** API key for Materials Project.
**Output:** JSON or CSV files containing the fetched structures and their descriptors.

### 8. UMAP.py
**Function:** Applies UMAP dimensionality reduction to a dataset to visualize high-dimensional data in two dimensions.
**Input:** CSV file containing materials property data.
**Output:** Plots a UMAP projection, potentially saving the plot as an image file.

### 9. VariablesCorrelations.py
**Function:** Analyzes correlations between different materials properties and visualizes these relationships.
**Input:** Data file containing the properties to be analyzed.
**Output:** Plots of the correlation matrices and any significant findings.

### 10. pDOS-StatesOverlap.py
**Function:** Calculates the overlap between projected density of states (pDOS) for different elements and orbitals using integration.
**Input:** `vasprun.xml` file.
**Output:** Outputs the calculated overlaps to standard output or a file.

### 11. pDOS_collect_plot.py
**Function:** Collects pDOS data from multiple VASP calculations, normalizing and plotting results for comparison.
**Input:** Directory containing multiple `vasprun.xml` files.
**Output:** Plots normalized pDOS and saves figures in the directory.

### 12. radial_distr_funct.py
**Function:** Computes radial distribution functions from a POSCAR or CONTCAR file to analyze atomic distributions within a material.
**Input:** POSCAR or CONTCAR file.
**Output:** Plots the radial distribution function and optionally saves the plot as an image file.
