# ORCA tools

### IR_Raman_extract_plot.py
**Function:** Parses the `IR SPECTRUM` and `RAMAN SPECTRUM` blocks of an ORCA output file and plots both spectra with Lorentzian broadening.
**Usage:** `python IR_Raman_extract_plot.py molecule.out [--fwhm 10 --xmin 200 --xmax 2000 --out spectra.png]`
**Input:** ORCA `.out` file from a frequency calculation containing the IR and, if computed, Raman spectrum blocks.
**Output:** Figure with IR and Raman panels.

### gen_homo_lumo_auto.sh
**Function:** Reads the orbital occupations from an ORCA output file to find the HOMO and LUMO indices, then drives `orca_plot` non-interactively to write cube files for both orbitals.
**Usage:** `./gen_homo_lumo_auto.sh basename` (expects `basename.gbw` and `basename.out`; `orca_plot` must be on the PATH).
**Output:** Cube files for the HOMO and LUMO written by `orca_plot`.
