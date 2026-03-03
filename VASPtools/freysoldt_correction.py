import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pymatgen.io.vasp.outputs import Vasprun, Locpot

# ---------------------------------------------------------------------------
# Monkey patch for doped structure-matcher issue (NV center compatibility)
# ---------------------------------------------------------------------------
import doped.utils.parsing as _doped_parsing

def _patched_get_species(composition_diff, el_change):
    try:
        return next(el for el, amt in composition_diff.items() if amt == el_change)
    except StopIteration:
        pass
    if el_change < 0:
        return min(composition_diff, key=lambda el: composition_diff[el])
    else:
        return max(composition_diff, key=lambda el: composition_diff[el])

_doped_parsing._get_species_from_composition_diff = _patched_get_species
print("Monkey-patch applied: doped structure-matcher fixed.")

from doped.analysis import DefectParser

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
bulk_path = "/home/562/mxf562/storagejm11/DATA/Diamond_Bulk/444_pristine_spin2/neutral/PBE"
nv0_path  = "/home/562/mxf562/storagejm11/DATA/Diamond_Bulk/444_gamma_0_PBE_HSE_CAM/PBE/LOC"
nvm_path  = "/home/562/mxf562/storagejm11/DATA/Diamond_Bulk/444_gamma_PBE_HSE_CAM/PBE"

DIELECTRIC = 5.7

# ---------------------------------------------------------------------------
# 1. Parse bulk
# ---------------------------------------------------------------------------
bulk_vr = Vasprun(f"{bulk_path}/vasprun.xml")
E_bulk = bulk_vr.final_energy

band_gap, cbm, vbm, is_direct = bulk_vr.eigenvalue_band_properties

print("\n--- Bulk properties ---")
print(f"VBM:  {vbm:.6f} eV")
print(f"CBM:  {cbm:.6f} eV")
print(f"Gap:  {band_gap:.6f} eV")
print(f"Direct gap: {is_direct}")

# ---------------------------------------------------------------------------
# 2. Parse NV0
# ---------------------------------------------------------------------------
print("\nParsing NV0 (q=0)...")
nv0_entry = DefectParser.from_paths(
    defect_path=nv0_path,
    bulk_path=bulk_path,
    dielectric=DIELECTRIC,
    charge_state=0,
    skip_corrections=True,
).defect_entry

E_nv0 = nv0_entry.sc_entry.energy

print(f"NV0 energy: {E_nv0:.6f} eV")

# ---------------------------------------------------------------------------
# 3. Parse NV-
# ---------------------------------------------------------------------------
print("\nParsing NV- (q=-1)...")
nvm_entry = DefectParser.from_paths(
    defect_path=nvm_path,
    bulk_path=bulk_path,
    dielectric=DIELECTRIC,
    charge_state=-1,
    skip_corrections=True,
).defect_entry

E_nvm = nvm_entry.sc_entry.energy

print(f"NV- energy: {E_nvm:.6f} eV")

# ---------------------------------------------------------------------------
# 4. Load LOCPOTs and apply Freysoldt correction
# ---------------------------------------------------------------------------
print("\nLoading LOCPOT files...")
bulk_locpot = Locpot.from_file(f"{bulk_path}/LOCPOT")
nvm_locpot  = Locpot.from_file(f"{nvm_path}/LOCPOT")

print("Applying Freysoldt correction...")
correction, fig = nvm_entry.get_freysoldt_correction(
    dielectric=DIELECTRIC,
    defect_locpot=nvm_locpot,
    bulk_locpot=bulk_locpot,
    plot=True,
)

E_corr = correction.correction_energy

fig.savefig("freysoldt_correction_NVm.png", dpi=150, bbox_inches="tight")

print(f"Freysoldt correction: {E_corr:.6f} eV")

# ---------------------------------------------------------------------------
# 5. Raw energy differences
# ---------------------------------------------------------------------------
print("\n--- Energy differences ---")
print(f"NV0 - bulk: {E_nv0 - E_bulk:.6f} eV")
print(f"NV- - bulk: {E_nvm - E_bulk:.6f} eV")

# 6. Charge transition level ε(0/-)
# First compute absolute Fermi level (relative to VASP internal reference),
# then convert to "above VBM" by subtracting the bulk VBM.

epsilon_abs = (E_nvm + E_corr) - E_nv0          # absolute reference (VASP zero)
epsilon_0m  = epsilon_abs - vbm                 # relative to VBM

print("\n--- Charge Transition Level ε(0/-) ---")
print(f"ε_abs (VASP ref): {epsilon_abs:.6f} eV")
print(f"ε(0/-) above VBM: {epsilon_0m:.6f} eV")
print(f"As % of gap: {epsilon_0m / band_gap * 100:.2f} %")

if 0 < epsilon_0m < band_gap:
    print("✓ ε(0/-) lies inside the band gap")
else:
    print("⚠ ε(0/-) outside band gap → check alignment / convergence")

# ---------------------------------------------------------------------------
# 7. Optional: print KS HOMO-LUMO gap
# ---------------------------------------------------------------------------
print("\n--- Sanity check ---")
print(f"Kohn-Sham gap: {band_gap:.6f} eV")
