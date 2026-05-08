import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # non-interactive backend (Gadi/Setonix safe)
import matplotlib.pyplot as plt
import os


# ----------------------------------------------------------------------
# I/O
# ----------------------------------------------------------------------
def parse_eigenval_to_dataframe(file_path):
    """
    Parse a spin-polarised VASP EIGENVAL file.

    Columns: band number, spin-up energy, spin-down energy,
             spin-up occupancy, spin-down occupancy.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # metadata is on line 6 (0-indexed line 5): NELECT, NKPTS, NBANDS
    num_electrons, num_kpoints, num_bands = map(int, lines[5].split())

    data = []
    cur = 7  # first k-point header line
    for kpt in range(1, num_kpoints + 1):
        cur += 1  # skip k-point header
        for _ in range(num_bands):
            parts = lines[cur].split()
            data.append([
                kpt,
                int(parts[0]),
                float(parts[1]),
                float(parts[2]),
                float(parts[3]),
                float(parts[4]),
            ])
            cur += 1

    df = pd.DataFrame(data, columns=[
        'k-point', 'band number',
        'spin-up energy', 'spin-down energy',
        'spin-up occupancy', 'spin-down occupancy',
    ])
    return df, num_electrons, num_bands


def get_fermi_energy_from_outcar(file_path='OUTCAR'):
    """Extract E-fermi from OUTCAR (last occurrence)."""
    ef = None
    with open(file_path, 'r') as f:
        for line in f:
            if 'E-fermi' in line:
                ef = float(line.split()[2])
    if ef is None:
        raise ValueError("Fermi energy not found in OUTCAR")
    return ef


# ----------------------------------------------------------------------
# Frontier / VBM detection (open-shell aware)
# ----------------------------------------------------------------------
def find_bulk_vbm_band(df, kpt=1, exch_tol=0.02, occ_tol=0.5):
    """
    Highest band that is bulk-like at the given k-point:
      - doubly occupied (occ_up > occ_tol AND occ_dn > occ_tol)
      - small exchange splitting |eps_up - eps_dn| < exch_tol  (eV)

    For a spin-polarised defect supercell this skips singly-occupied
    in-gap states (e.g. NV-like a1, e) and partially split defect levels,
    returning the highest genuinely bulk-like valence band.
    """
    g = df[df['k-point'] == kpt].copy()
    g['exch'] = (g['spin-up energy'] - g['spin-down energy']).abs()
    mask = ((g['spin-up occupancy']   > occ_tol) &
            (g['spin-down occupancy'] > occ_tol) &
            (g['exch'] < exch_tol))
    if not mask.any():
        raise RuntimeError("No bulk-like VBM band found; loosen exch_tol.")
    return int(g.loc[mask, 'band number'].max())


def find_homo_lumo(df, kpt=1, occ_tol=0.5):
    """
    HOMO/LUMO per spin channel, derived from occupations (open-shell safe).

    Returns (homo_up, lumo_up, homo_dn, lumo_dn) as band indices.
    """
    g = df[df['k-point'] == kpt]
    homo_up = int(g[g['spin-up occupancy']   > occ_tol]['band number'].max())
    lumo_up = int(g[g['spin-up occupancy']   < occ_tol]['band number'].min())
    homo_dn = int(g[g['spin-down occupancy'] > occ_tol]['band number'].max())
    lumo_dn = int(g[g['spin-down occupancy'] < occ_tol]['band number'].min())
    return homo_up, lumo_up, homo_dn, lumo_dn


def report_homo_lumo(df, kpt=1):
    hu, lu, hd, ld = find_homo_lumo(df, kpt=kpt)
    g = df[df['k-point'] == kpt].set_index('band number')
    print(f"HOMO (up)   = band {hu}, {g.loc[hu, 'spin-up energy']:.4f} eV   "
          f"LUMO (up)   = band {lu}, {g.loc[lu, 'spin-up energy']:.4f} eV")
    print(f"HOMO (down) = band {hd}, {g.loc[hd, 'spin-down energy']:.4f} eV   "
          f"LUMO (down) = band {ld}, {g.loc[ld, 'spin-down energy']:.4f} eV")
    return hu, lu, hd, ld


# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------
def plot_energy_levels(df, shift_energy, name,
                       band_window, kpt=1,
                       y_range=(-6, 8), occ_tol=0.5):
    """
    Plot Kohn-Sham levels in the gap region.

    Filled (closed) marker  -> occupied  (occ > occ_tol)
    Hollow (open)   marker  -> unoccupied (occ <= occ_tol)
    Spin-up   = green up-triangle
    Spin-down = magenta down-triangle

    `shift_energy` is subtracted from both spin channels so VBM lies at 0.
    Does NOT mutate the input DataFrame.
    """
    bandINF, bandSUP = band_window
    sub = df[(df['k-point'] == kpt) &
             (df['band number'] >= bandINF) &
             (df['band number'] <= bandSUP)].copy()
    sub['spin-up energy']   = sub['spin-up energy']   - shift_energy
    sub['spin-down energy'] = sub['spin-down energy'] - shift_energy

    plt.figure(figsize=(8, 6))
    seen_up = seen_dn = False

    for _, row in sub.iterrows():
        b   = row['band number']
        eu  = row['spin-up energy']
        ed  = row['spin-down energy']
        ou  = row['spin-up occupancy']   > occ_tol
        od  = row['spin-down occupancy'] > occ_tol

        # spin-up
        if ou:
            plt.scatter(b, eu, color='green', marker='^', s=200,
                        label='spin-up' if not seen_up else "")
        else:
            plt.scatter(b, eu, facecolors='none', edgecolors='green',
                        marker='^', s=200, linewidths=2,
                        label='spin-up' if not seen_up else "")
        plt.text(b, eu - 0.5, f"{eu:.2f}", color='green',
                 fontsize=12, ha='center', va='top')
        seen_up = True

        # spin-down
        if od:
            plt.scatter(b, ed, color='magenta', marker='v', s=200,
                        label='spin-down' if not seen_dn else "")
        else:
            plt.scatter(b, ed, facecolors='none', edgecolors='magenta',
                        marker='v', s=200, linewidths=2,
                        label='spin-down' if not seen_dn else "")
        plt.text(b, ed + 0.5, f"{ed:.2f}", color='magenta',
                 fontsize=12, ha='center', va='bottom')
        seen_dn = True

    plt.axhline(0, color='gray', linestyle='--', label='VBM = 0 eV')
    plt.xlabel('Band Number', fontsize=18)
    plt.ylabel('Energy (eV)', fontsize=18)
    plt.xlim(bandINF, bandSUP)
    plt.xticks(range(bandINF, bandSUP + 1), fontsize=16)
    plt.yticks(fontsize=14)
    plt.ylim(*y_range)
    plt.legend(loc='lower left', fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    out = f"{name}_shifted.png"
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"Wrote {out}")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
if __name__ == "__main__":
    cwd = os.getcwd()
    name = os.path.basename(cwd)

    # --- parse EIGENVAL ---
    df, nelec, nbands = parse_eigenval_to_dataframe('EIGENVAL')
    print(f"NELECT = {nelec}, NBANDS = {nbands}")
    print(df.head())

    # --- Fermi energy (informational only) ---
    try:
        ef = get_fermi_energy_from_outcar('OUTCAR')
        print(f"E-fermi (OUTCAR): {ef:.4f} eV")
    except FileNotFoundError:
        print("OUTCAR not found; skipping Fermi readout.")

    # --- locate the bulk VBM band automatically ---
    vbm_band = find_bulk_vbm_band(df, kpt=1)
    vbm = df[(df['k-point'] == 1) &
             (df['band number'] == vbm_band)]['spin-down energy'].iloc[0]
    print(f"Bulk-like VBM: band {vbm_band}, E = {vbm:.4f} eV (spin-down)")

    # --- HOMO/LUMO per channel (occupation-based) ---
    hu, lu, hd, ld = report_homo_lumo(df, kpt=1)

    # --- plot window: bracket VBM and the gap states ---
    # default: VBM-1 up to LUMO_up+1 (covers all in-gap defect states + CBM)
    bandINF = vbm_band - 1
    bandSUP = max(lu, ld) + 1
    print(f"Plot window: bands {bandINF}-{bandSUP}")

    plot_energy_levels(df, shift_energy=vbm, name=name,
                       band_window=(bandINF, bandSUP),
                       kpt=1, y_range=(-6, 8))
