# CompMatTools

Scripts I use day to day for first-principles materials modelling: preparing
structures, post-processing VASP, ORCA and Gaussian output, analysing defects
and heterostructures, and benchmarking parallel performance on HPC systems.
Everything here was written for my own research and is shared as is; each
folder has a README describing its scripts, inputs and outputs.

| Folder | Contents |
|---|---|
| [`VASPtools/`](VASPtools/) | 28 Python tools for VASP: structure generation (bilayers, supercells, ML datasets), band structures, total and projected DOS, Kohn–Sham level diagrams, defect analysis (in-gap states, wavefunction extraction, Freysoldt correction, Huang–Rhys factors), elastic properties, diffusion coefficients, radial distribution functions, phonon thermal properties and charge-transfer analysis. |
| [`ORCA/`](ORCA/) | IR/Raman spectrum extraction and plotting from ORCA output; automatic HOMO/LUMO cube-file generation with `orca_plot`. |
| [`gaussian16/`](gaussian16/) | HOMO–LUMO gap extraction from Gaussian 16 log files. |
| [`HPCbenchmarking/`](HPCbenchmarking/) | Speed-up and parallel-efficiency analysis for VASP SCF timings across cores and nodes, and OpenMP scaling. |
| [`scripts/`](scripts/) | Shell helpers, e.g. batch rendering of wavefunction isosurfaces with VESTA. |

## Requirements

Python 3.9+ with `numpy`, `pandas`, `matplotlib`, `pymatgen`, `ase`. Some
tools need extras noted in their README entry: `scikit-learn`, `matminer`,
`doped`, `pyyaml`, VASPKIT, Phonopy. Install the common set with

```bash
pip install numpy pandas matplotlib pymatgen ase scikit-learn pyyaml
```

## History

This repository absorbs two earlier ones, `VASPpostproc` and
`POSCAR-2D-utilities`, whose scripts now live in `VASPtools/`.

## Author

Marco Fronzi, School of Physics, The University of Sydney —
marco.fronzi@sydney.edu.au
