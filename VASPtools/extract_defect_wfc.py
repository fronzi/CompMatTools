#!/usr/bin/env python3
"""
extract_defect_wfc.py

Extract real-space pseudo-wavefunctions of in-gap defect levels for both
spin channels, using VASPKIT task 511 (Wave-Function Analysis).

Workflow
--------
1. Reads EIGENVAL to identify:
     - bulk-like VBM band (highest doubly occupied, ~0 exchange splitting)
     - bulk-like CBM band (lowest doubly empty, ~0 exchange splitting)
     - defect bands = every band strictly between VBM and CBM
2. For each defect band x each spin channel (1 = up, 2 = down), drives
   `vaspkit -task 511` via stdin and renames the resulting *.vasp file(s).

Run in a directory that contains: EIGENVAL, WAVECAR, POSCAR.
Outputs:  <name>_band<N>_spin<S>_<orig>.vasp  (e.g. NLi0_GS_band517_spin1_WFN_REAL.vasp)

Notes
-----
- For HSE06 / gamma-only WAVECAR, vaspkit auto-detects via the WAVECAR
  header. No extra flags needed.
- If your vaspkit version uses a different prompt sequence for task 511,
  edit `VASPKIT_INPUT_TEMPLATE` below (a stale menu mapping is the most
  common failure mode across vaspkit releases).
- This script does NOT regenerate WAVECAR: it requires a converged
  WAVECAR from the GS calculation in the cwd.
"""

import os
import re
import sys
import glob
import shutil
import subprocess
from pathlib import Path
import pandas as pd


# --------------------------------------------------------------------------
# EIGENVAL parsing & frontier detection (open-shell aware)
# --------------------------------------------------------------------------
def parse_eigenval(path='EIGENVAL'):
    with open(path) as f:
        lines = f.readlines()
    nelec, nkpt, nbnd = map(int, lines[5].split())
    rows, cur = [], 7
    for kpt in range(1, nkpt + 1):
        cur += 1
        for _ in range(nbnd):
            p = lines[cur].split()
            rows.append([kpt, int(p[0]), float(p[1]), float(p[2]),
                         float(p[3]), float(p[4])])
            cur += 1
    df = pd.DataFrame(rows, columns=[
        'k-point', 'band number',
        'spin-up energy', 'spin-down energy',
        'spin-up occupancy', 'spin-down occupancy'])
    return df, nelec, nbnd


def find_bulk_vbm_band(df, kpt=1, exch_tol=0.02, occ_tol=0.5):
    g = df[df['k-point'] == kpt].copy()
    g['exch'] = (g['spin-up energy'] - g['spin-down energy']).abs()
    mask = ((g['spin-up occupancy']   > occ_tol) &
            (g['spin-down occupancy'] > occ_tol) &
            (g['exch'] < exch_tol))
    if not mask.any():
        raise RuntimeError("No bulk-like VBM found; loosen exch_tol.")
    return int(g.loc[mask, 'band number'].max())


def find_bulk_cbm_band(df, kpt=1, exch_tol=0.05, occ_tol=0.5):
    g = df[df['k-point'] == kpt].copy()
    g['exch'] = (g['spin-up energy'] - g['spin-down energy']).abs()
    mask = ((g['spin-up occupancy']   < occ_tol) &
            (g['spin-down occupancy'] < occ_tol) &
            (g['exch'] < exch_tol))
    if not mask.any():
        raise RuntimeError("No bulk-like CBM found; loosen exch_tol.")
    return int(g.loc[mask, 'band number'].min())


def defect_bands(df, vbm_band, cbm_band):
    """All bands strictly between bulk VBM and bulk CBM."""
    return list(range(vbm_band + 1, cbm_band))


# --------------------------------------------------------------------------
# VASPKIT driver
# --------------------------------------------------------------------------
# Menu sequence for task 511 in vaspkit 1.3+:
#   51   -> Wave-Function Analysis
#   511  -> Real-space wavefunction of selected (spin, kpt, band)
#   <ispin> <ikpt> <iband>
# Trailing 0\n exits cleanly.
VASPKIT_INPUT_TEMPLATE = "51\n511\n{ispin}\n{ikpt}\n{iband}\n0\n"


def run_vaspkit_511(ispin, ikpt, iband, vaspkit_cmd='vaspkit',
                    timeout=900, log_file=None):
    """Drive vaspkit task 511 via stdin. Returns CompletedProcess."""
    stdin = VASPKIT_INPUT_TEMPLATE.format(ispin=ispin, ikpt=ikpt, iband=iband)
    res = subprocess.run([vaspkit_cmd], input=stdin,
                         capture_output=True, text=True, timeout=timeout)
    if log_file is not None:
        with open(log_file, 'a') as fh:
            fh.write(f"\n===== spin={ispin} kpt={ikpt} band={iband} =====\n")
            fh.write(res.stdout)
            if res.stderr:
                fh.write("\n--- STDERR ---\n" + res.stderr)
    return res


def collect_new_vasp_files(before):
    """Set of *.vasp files in cwd minus the snapshot 'before'."""
    return set(glob.glob('*.vasp')) - before


def safe_rename(src, dst):
    """Rename, refusing to overwrite."""
    if Path(dst).exists():
        print(f"  ! skip rename: {dst} already exists")
        return
    os.rename(src, dst)
    print(f"  -> {dst}")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    cwd = Path.cwd()
    name = cwd.name

    for required in ('EIGENVAL', 'WAVECAR', 'POSCAR'):
        if not (cwd / required).exists():
            sys.exit(f"ERROR: {required} not found in {cwd}")

    if shutil.which('vaspkit') is None:
        sys.exit("ERROR: 'vaspkit' not on PATH. Load the module or update PATH.")

    # --- locate frontier and defect bands ---
    df, nelec, nbnd = parse_eigenval('EIGENVAL')
    vbm = find_bulk_vbm_band(df)
    cbm = find_bulk_cbm_band(df)
    bands = defect_bands(df, vbm, cbm)

    print(f"NELECT = {nelec}, NBANDS = {nbnd}")
    print(f"Bulk-like VBM band: {vbm}")
    print(f"Bulk-like CBM band: {cbm}")
    print(f"Defect bands to extract: {bands}")
    if not bands:
        sys.exit("No defect bands found between VBM and CBM. Nothing to do.")

    log_file = cwd / 'vaspkit_511.log'
    if log_file.exists():
        log_file.unlink()

    # --- loop over (band, spin); rename outputs ---
    for iband in bands:
        for ispin in (1, 2):
            print(f"\n>> vaspkit 511 :: band {iband}, spin {ispin}")
            before = set(glob.glob('*.vasp'))
            try:
                res = run_vaspkit_511(ispin=ispin, ikpt=1, iband=iband,
                                       log_file=log_file)
            except subprocess.TimeoutExpired:
                print(f"  ! timeout for band {iband}, spin {ispin}")
                continue
            if res.returncode != 0:
                print(f"  ! vaspkit returned {res.returncode} "
                      f"(see {log_file.name})")
                continue

            new_files = sorted(collect_new_vasp_files(before))
            if not new_files:
                print("  ! no new .vasp files appeared "
                      "(check vaspkit log; menu may differ in your version)")
                continue

            for f in new_files:
                tag = Path(f).stem  # strip .vasp
                target = f"{name}_band{iband}_spin{ispin}_{tag}.vasp"
                safe_rename(f, target)

    print("\nDone. Wavefunction files written to:", cwd)
    print(f"VASPKIT log: {log_file}")


if __name__ == '__main__':
    main()
