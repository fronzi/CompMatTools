#!/usr/bin/env python3
"""
orca_spectra.py — Parse and plot IR and Raman spectra from ORCA output files.

Usage
-----
    python orca_spectra.py molecule.out
    python orca_spectra.py molecule.out --fwhm 10 --xmin 200 --xmax 2000
    python orca_spectra.py molecule.out --out spectra.png

ORCA section column layouts
----------------------------
IR SPECTRUM:
    Mode  freq (cm**-1)  T**2  TX  TY  TZ
    → freq_col=1, val_col=2 (T**2 = total dipole derivative squared, km/mol)

RAMAN SPECTRUM:
    Mode  freq (cm**-1)  Activity  Depolar
    → freq_col=1, val_col=2 (Activity, Å⁴/amu)
"""

import argparse
import re
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def parse_orca_spectrum(filepath: str):
    """
    Extract IR and Raman spectral data from an ORCA output file.

    The function locates the 'IR SPECTRUM' and 'RAMAN SPECTRUM' blocks by
    scanning for the section header, counting two dashed-separator lines to
    find where the data rows begin, and collecting lines of the form 'N: ...'
    until an empty line signals the end of the block.

    Parameters
    ----------
    filepath : str
        Path to an ORCA .out or .log file.

    Returns
    -------
    f_ir : ndarray   shape (N_ir,)   — IR mode frequencies (cm**-1)
    ir_int : ndarray shape (N_ir,)   — IR intensities, T**2 (km/mol)
    f_ram : ndarray  shape (N_ram,)  — Raman frequencies (cm**-1)
    ram_act : ndarray shape (N_ram,) — Raman activities (Å⁴/amu)
    """
    with open(filepath, 'r') as fh:
        content = fh.read()

    def _parse_section(header: str, freq_col: int, val_col: int):
        idx = content.find(header)
        if idx == -1:
            print(f"  [WARNING] Section '{header}' not found in {filepath}.",
                  file=sys.stderr)
            return np.array([]), np.array([])

        freqs, vals = [], []
        dash_count = 0
        in_data = False

        for line in content[idx:].split('\n'):
            stripped = line.strip()

            # Dashed separator lines bracket the column-header row;
            # data begins immediately after the second one.
            if re.fullmatch(r'-{10,}', stripped):
                dash_count += 1
                if dash_count == 2:
                    in_data = True
                continue

            if in_data:
                if not stripped:           # blank line = end of block
                    break
                parts = stripped.split()
                # Data rows begin with 'N:' (e.g. '  6:')
                if parts and re.fullmatch(r'\d+:', parts[0]):
                    try:
                        freqs.append(float(parts[freq_col]))
                        vals.append(float(parts[val_col]))
                    except (IndexError, ValueError):
                        continue

        return np.array(freqs), np.array(vals)

    # IR: [mode:] [freq] [T**2] [TX] [TY] [TZ]  → val_col=2 is T**2
    f_ir,  ir_int  = _parse_section('IR SPECTRUM',    freq_col=1, val_col=2)
    # Raman: [mode:] [freq] [Activity] [Depolar]  → val_col=2 is Activity
    f_ram, ram_act = _parse_section('RAMAN SPECTRUM', freq_col=1, val_col=2)

    # Filter out imaginary (negative) frequencies
    if len(f_ir):
        mask = f_ir > 0
        f_ir, ir_int = f_ir[mask], ir_int[mask]
    if len(f_ram):
        mask = f_ram > 0
        f_ram, ram_act = f_ram[mask], ram_act[mask]

    return f_ir, ir_int, f_ram, ram_act


def extract_title(filepath: str) -> str:
    """
    Attempt to extract a meaningful molecule/job label from the ORCA header.
    Falls back to the file stem.
    """
    stem = os.path.splitext(os.path.basename(filepath))[0]
    try:
        with open(filepath, 'r') as fh:
            for line in fh:
                # ORCA prints "The Molecule has ..." or the input block title
                m = re.search(r'^\s*\|\s+(.+?)\s+\|', line)
                if m and 'ORCA' not in m.group(1):
                    candidate = m.group(1).strip()
                    if candidate and len(candidate) < 60:
                        return candidate
    except Exception:
        pass
    return stem


# ---------------------------------------------------------------------------
# Spectrum construction
# ---------------------------------------------------------------------------

def lorentzian(x: np.ndarray, x0: float, fwhm: float) -> np.ndarray:
    gamma = fwhm / 2.0
    return gamma**2 / ((x - x0)**2 + gamma**2)


def build_broadened(freqs, intensities, x, fwhm=15.0):
    if len(freqs) == 0:
        return np.zeros_like(x)
    return np.sum(
        v * lorentzian(x, f, fwhm)
        for f, v in zip(freqs, intensities)
    )


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_spectra(
    f_ir, ir_int, f_ram, ram_act,
    title='ORCA Spectra',
    fwhm=15.0,
    xmin=None, xmax=None,
    outfile=None,
):
    # Auto x-range from parsed data with 100 cm-1 padding
    all_f = np.concatenate([f_ir, f_ram]) if (len(f_ir) or len(f_ram)) else np.array([0.0, 4000.0])
    if xmin is None:
        xmin = max(0.0, all_f.min() - 100.0)
    if xmax is None:
        xmax = all_f.max() + 100.0

    x = np.linspace(xmin, xmax, 3000)
    ir_y  = build_broadened(f_ir,  ir_int,  x, fwhm)
    ram_y = build_broadened(f_ram, ram_act, x, fwhm)

    # Normalised stick heights scaled to match the broadened envelope
    def stick_heights(freqs, vals, envelope):
        if len(freqs) == 0 or envelope.max() == 0:
            return np.zeros_like(vals)
        return vals / vals.max() * envelope.max() if vals.max() > 0 else vals

    ir_sticks  = stick_heights(f_ir,  ir_int,  ir_y)
    ram_sticks = stick_heights(f_ram, ram_act, ram_y)

    # ---- Layout ------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(11, 8), sharex=True,
        gridspec_kw={'hspace': 0.06}
    )
    fig.patch.set_facecolor('white')

    colors = {'ir': '#C0392B', 'ram': '#2471A3'}

    for ax, y, freqs, sticks, ylabel, color in [
        (ax1, ir_y,  f_ir,  ir_sticks,  'IR Intensity (km mol$^{-1}$)',       colors['ir']),
        (ax2, ram_y, f_ram, ram_sticks, 'Raman Activity (Å$^4$ amu$^{-1}$)', colors['ram']),
    ]:
        ax.plot(x, y, color=color, lw=2.0, zorder=3)
        ax.fill_between(x, y, color=color, alpha=0.12, zorder=2)

        # Stick spectrum
        if len(freqs):
            ax.vlines(
                freqs, 0, sticks,
                color=color, lw=1.0, alpha=0.55,
                linestyle='--', zorder=4
            )

        ax.set_ylabel(ylabel, fontsize=12, labelpad=8)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(bottom=0)
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3g'))
        ax.tick_params(axis='both', labelsize=11)
        ax.grid(True, linestyle=':', alpha=0.35, zorder=0)
        for sp in ['top', 'right']:
            ax.spines[sp].set_visible(False)

    ax1.set_title(title, fontsize=14, fontweight='bold', pad=10)
    ax2.set_xlabel('Wavenumber (cm$^{-1}$)', fontsize=12, labelpad=8)

    # Shared legend annotation
    fig.text(
        0.97, 0.02,
        f'Lorentzian FWHM = {fwhm:.0f} cm$^{{-1}}$',
        ha='right', va='bottom', fontsize=9,
        color='gray', style='italic'
    )

    plt.tight_layout(rect=[0, 0.02, 1, 1])

    if outfile:
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
        print(f"Figure saved → {outfile}")
    else:
        plt.show()


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Plot IR and Raman spectra from an ORCA output file.'
    )
    parser.add_argument('orca_file',
                        help='Path to ORCA .out or .log file')
    parser.add_argument('--fwhm', type=float, default=15.0,
                        help='Lorentzian FWHM in cm**-1 (default: 15)')
    parser.add_argument('--xmin', type=float, default=None,
                        help='Lower x-axis bound in cm**-1 (auto if omitted)')
    parser.add_argument('--xmax', type=float, default=None,
                        help='Upper x-axis bound in cm**-1 (auto if omitted)')
    parser.add_argument('--out', type=str, default=None,
                        help='Save figure to file, e.g. spectra.png')
    args = parser.parse_args()

    if not os.path.isfile(args.orca_file):
        sys.exit(f"[ERROR] File not found: {args.orca_file}")

    print(f"Parsing: {args.orca_file}")
    f_ir, ir_int, f_ram, ram_act = parse_orca_spectrum(args.orca_file)

    print(f"  IR modes  : {len(f_ir)}"
          + (f"  (freq range {f_ir.min():.1f}–{f_ir.max():.1f} cm**-1)" if len(f_ir) else ""))
    print(f"  Raman modes: {len(f_ram)}"
          + (f"  (freq range {f_ram.min():.1f}–{f_ram.max():.1f} cm**-1)" if len(f_ram) else ""))

    mol   = extract_title(args.orca_file)
    title = f'IR & Raman Spectra — {mol}  (ORCA, Lorentzian)'

    plot_spectra(
        f_ir, ir_int, f_ram, ram_act,
        title=title, fwhm=args.fwhm,
        xmin=args.xmin, xmax=args.xmax,
        outfile=args.out,
    )


if __name__ == '__main__':
    main()
