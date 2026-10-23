#  VASP pre and post processing tools

## Overview




This repository contains Python scripts designed for analyzing and processing electronic structure data from computational materials science simulations using VASP outputs.


### 0. poscar-interlayer-distances.py
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

### 13. BandStructuresPlotOverlap.py
**Function:** Overlays band structures from several VASP calculations on one plot, aligning each to its valence-band maximum; handles spin-polarised and non-spin-polarised runs.
**Input:** `vasprun.xml` from each calculation (the script is set up for a PBE and an HSE run in `PBE-prim/` and `HSE-prim/`; edit the paths).
**Output:** Combined band-structure figure.

### 14. plot_bs.py
**Function:** Minimal band-structure plot for a single calculation using pymatgen's `BSPlotter`, with thicker lines for publication.
**Input:** `vasprun.xml`, `KPOINTS`.
**Output:** Band-structure figure.

### 15. pDOS.py
**Function:** Total and projected DOS plot for a single calculation (element- or orbital-resolved) with configurable energy window.
**Input:** `vasprun.xml` (LORBIT = 11).
**Output:** pDOS figure.

### 16. plot_pdos_template.py
**Function:** Publication-quality projected DOS across one or several systems with a built-in `vasprun.xml` parser (pymatgen optional): spin-mirrored panels, colour-blind-safe palette, per-group magnification shown in the legend, antiferromagnetic sublattice separation from OUTCAR moments, and band-centre calculation with an explicit, always-reported energy window. Run `python plot_pdos_template.py -h`; the module docstring documents the conventions.
**Input:** One or more calculation directories (`vasprun.xml`, or `DOSCAR` + `POSCAR`; `OUTCAR` optional).
**Output:** Figure (PDF/PNG), band centres on stdout and as CSV.

### 17. plot_energy_levels.py and plot_defect_levels.py
**Function:** Parse `EIGENVAL` (NELECT, NKPTS, NBANDS from the header) and plot Kohn–Sham eigenvalues per k-point and spin channel; `plot_defect_levels.py` adds occupation-aware styling to highlight in-gap defect states.
**Input:** `EIGENVAL`.
**Output:** Energy-level diagrams.

### 18. level_scheme_template.py
**Function:** Vector-PDF Kohn–Sham level-scheme template (valence/conduction continua as gradients, one column per spin channel, level bars with occupation markers), in the format used for the diamond colour-centre figures.
**Input:** Level energies and occupations specified in the script's configuration block.
**Output:** PDF figure.

### 19. extract_defect_wfc.py
**Function:** Identifies in-gap defect levels for both spin channels from `EIGENVAL` and drives VASPKIT task 511 to extract their real-space pseudo-wavefunctions.
**Input:** `EIGENVAL`, `WAVECAR`; requires `vaspkit` on the PATH.
**Output:** Wavefunction `.vasp` files per defect level.

### 20. freysoldt_correction.py
**Function:** Freysoldt–Neugebauer–Van de Walle charge-correction for charged defects using the `doped` package, with a structure-matcher patch for NV-centre-type cells.
**Input:** Defect and bulk calculation directories (`LOCPOT`, `OUTCAR`).
**Output:** Corrected formation energies and potential-alignment plots.

### 21. Huang-Rhys.py
**Function:** Huang–Rhys factor calculator from VASP or Quantum ESPRESSO outputs supporting the energy method (four single-point energies), the displacement method (ΔQ between ground- and excited-state geometries) and the mode-resolved method (phonon projection). See the module docstring for the three workflows.
**Input:** Relaxed ground- and excited-state structures and energies; optional phonon eigenvectors.
**Output:** Huang–Rhys factor S, effective phonon energy, and partial HR factors per mode.

### 22. multiple_phonopy.py
**Function:** Plots Phonopy thermal properties from several calculation directories on common axes with a clean, Phonopy-like style.
**Input:** `thermal_properties.yaml` in each directory.
**Output:** Thermal-property figure.

### 23. GenerateBilayers.py
**Function:** Finds commensurate (Zur–McGill, via pymatgen's `ZSLGenerator`) or near-commensurate supercells for two 2D layers (default WO₃ // graphene), writes the matched supercell POSCARs and optional stacked bilayer, and logs every candidate to `bilayers_summary.csv`.
**Input:** Two POSCAR files; tolerance settings in the script.
**Output:** `POSCAR_*_SC` files, optional bilayer POSCAR, CSV log.

### 24. get_bilayer_info.py
**Function:** Reads a relaxed bilayer `CONTCAR`, separates the two layers by k-means clustering of the z-coordinates, and prints the interlayer spacing and slab thickness.
**Input:** `CONTCAR`.
**Output:** Interlayer distance and slab thickness to stdout.

### 25. ChargeTransfer.py
**Function:** Integrates `CHGCAR` charge density along z for a metal/MoS₂ heterostructure series to estimate charge transfer between the layers (uses NELECT from `OUTCAR` for normalisation). Edit the `M` list and directory pattern for other systems.
**Input:** `CHGCAR`, `CONTCAR`, `OUTCAR` in one directory per system.
**Output:** Charge above and below the interface and the net charge transfer, printed per system.

### 26. ChemGroupCount.py
**Function:** Counts chemical groups (hydroxyl, CH, CH₂, CH₃, CH₄) in a POSCAR/CONTCAR from neighbour distances, with adjustable O–H and C–H cutoffs. `python ChemGroupCount.py CONTCAR --oh_cutoff 1.2 --ch_cutoff 1.2`.
**Input:** POSCAR-format structure file.
**Output:** Group counts to stdout.

### 27. ElasticTensorML.py
**Function:** Trains a random-forest regressor on the matminer elastic-tensor dataset using composition and structure featurisers, with cross-validation and a residual plot. A machine-learning counterpart to `ElasticProperties.py`.
**Input:** None (dataset is downloaded through matminer).
**Output:** Cross-validation scores and residual plot.

### 28. generate_bilayer.ipynb
Scratch notebook used while developing `GenerateBilayers.py`.
