import numpy as np
import matplotlib
matplotlib.use("Agg")   # remove if running locally with a display
import matplotlib.pyplot as plt
from pymatgen.io.vasp.outputs import Vasprun, Locpot

# ---------------------------------------------------------------------------
# MONKEY-PATCH: fix doped's _get_species_from_composition_diff so it doesn't
# crash with StopIteration on complex defects like the NV center.
# Must be applied BEFORE importing DefectParser.
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
print("Monkey-patch applied: doped structure-matcher fixed for complex defects.")

from doped.analysis import DefectParser

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
bulk_path = "/home/562/mxf562/storagejm11/DATA/Diamond_Bulk/444_pristine_spin2/neutral/PBE"
nv0_path  = "/home/562/mxf562/storagejm11/DATA/Diamond_Bulk/444_gamma_0_PBE_HSE_CAM/PBE/LOC"
nvm_path  = "/home/562/mxf562/storagejm11/DATA/Diamond_Bulk/444_gamma_PBE_HSE_CAM/PBE/LOC"

DIELECTRIC = 5.7   # isotropic dielectric constant for diamond

# ---------------------------------------------------------------------------
# 1. Parse NV0 (q=0)
# ---------------------------------------------------------------------------
print("\nParsing NV0 (q=0) ...")
nv0_entry = DefectParser.from_paths(
    defect_path=nv0_path,
    bulk_path=bulk_path,
    dielectric=DIELECTRIC,
    charge_state=0,
    skip_corrections=True,
).defect_entry

print(f"  Name:   {nv0_entry.name}")
print(f"  Charge: {nv0_entry.charge_state:+}")
print(f"  Site:   {nv0_entry.defect_supercell_site.frac_coords}")

# ---------------------------------------------------------------------------
# 2. Parse NV- (q=-1)
# ---------------------------------------------------------------------------
print("\nParsing NV- (q=-1) ...")
nvm_entry = DefectParser.from_paths(
    defect_path=nvm_path,
    bulk_path=bulk_path,
    dielectric=DIELECTRIC,
    charge_state=-1,
    skip_corrections=True,
).defect_entry

print(f"  Name:   {nvm_entry.name}")
print(f"  Charge: {nvm_entry.charge_state:+}")
print(f"  Site:   {nvm_entry.defect_supercell_site.frac_coords}")

# ---------------------------------------------------------------------------
# 3. Load LOCPOTs explicitly
#    Required because skip_corrections=True means doped never loaded them.
#    The error message says: pass as argument or store in calculation_metadata.
#    Passing directly is the cleaner approach.
# ---------------------------------------------------------------------------
print("\nLoading LOCPOT files ...")
bulk_locpot = Locpot.from_file(f"{bulk_path}/LOCPOT")
nvm_locpot  = Locpot.from_file(f"{nvm_path}/LOCPOT")
print("  Done.")

# ---------------------------------------------------------------------------
# 4. Apply Freysoldt correction to NV- and generate plot
#    Mirrors the doped docs example:
#      correction, plot = F_O_1_entry.get_freysoldt_correction(plot=True)
#    The only difference is we pass LOCPOTs explicitly since skip_corrections=True
#    meant doped never cached them in calculation_metadata.
# ---------------------------------------------------------------------------
print("\nApplying Freysoldt correction to NV- ...")
correction, fig = nvm_entry.get_freysoldt_correction(
    dielectric=DIELECTRIC,
    defect_locpot=nvm_locpot,
    bulk_locpot=bulk_locpot,
    plot=True,
)

E_corr = correction.correction_energy
print(f"  Correction energy:  {E_corr:.4f} eV")
print(f"  Corrections dict:   {nvm_entry.corrections}")

fig.savefig("freysoldt_correction_NVm.png", dpi=150, bbox_inches="tight")
print("  Plot saved → freysoldt_correction_NVm.png")

# ---------------------------------------------------------------------------
# 5. Total energies
# ---------------------------------------------------------------------------
bulk_vr = Vasprun(f"{bulk_path}/vasprun.xml")
E_bulk  = bulk_vr.final_energy
E_nv0   = nv0_entry.sc_entry.energy
E_nvm   = nvm_entry.sc_entry.energy

print(f"\nBulk energy:     {E_bulk:.4f} eV")
print(f"NV0 energy:      {E_nv0:.4f} eV")
print(f"NV- energy:      {E_nvm:.4f} eV")
print(f"NV0 - bulk:      {E_nv0 - E_bulk:.4f} eV")
print(f"NV-  - bulk:     {E_nvm - E_bulk:.4f} eV")

# ---------------------------------------------------------------------------
# 6. Band edges from bulk
# ---------------------------------------------------------------------------
cbm = bulk_vr.eigenvalue_band_properties[0]
vbm = bulk_vr.eigenvalue_band_properties[1]
gap = bulk_vr.eigenvalue_band_properties[2]

print(f"\nVBM: {vbm:.4f} eV  |  CBM: {cbm:.4f} eV  |  Gap: {gap:.4f} eV")

# ---------------------------------------------------------------------------
# 7. Charge transition level ε(0/-)
#
#   Setting Ef(NV0) = Ef(NV-) and solving for ε_F (above VBM):
#
#     ε(0/-) = (E_nv0 - E_nvm) - E_corr
#
# ---------------------------------------------------------------------------
epsilon_0m = (E_nv0 - E_nvm) - E_corr

print(f"\nCharge transition level ε(0/-): {epsilon_0m:.4f} eV above VBM")
print(f"Gap = {gap:.4f} eV  →  ε(0/-) at {epsilon_0m / gap * 100:.1f}% of the gap")

if 0 < epsilon_0m < gap:
    print("-> ε(0/-) is inside the gap: physically meaningful result")
    if epsilon_0m > gap / 2:
        print("-> Upper half of gap: NV- is stable for n-type / mid-gap Fermi level")
    else:
        print("-> Lower half of gap: NV0 is stable for most Fermi level positions")
else:
    print("-> WARNING: ε(0/-) is outside the gap — check supercell sizes and potential alignment")
