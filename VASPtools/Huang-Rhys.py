"""
huang_rhys.py
=============
Huang-Rhys (HR) factor calculator from VASP and Quantum ESPRESSO outputs.

Supports three levels of theory:
  1. Energy method    : S = ΔE_relax / ħω_eff  (requires 4 single-point energies)
  2. Displacement method: S = ω_eff ΔQ² / 2ħ   (requires GS+ES geometries + ω_eff)
  3. Mode-resolved    : S_k = ω_k |⟨e_k|ΔR⟩|² / 2ħ  (requires phonon eigenvectors)

Configuration coordinate diagram:
  - Ground state parabola: E_gs(Q)
  - Excited state parabola: E_es(Q)
  - ZPL = E_es(R_es) - E_gs(R_gs)
  - ΔE_abs = E_es(R_gs) - E_gs(R_gs)  (vertical absorption)
  - ΔE_em  = E_es(R_es) - E_gs(R_es)  (vertical emission)
  - E_rel_es = E_es(R_gs) - E_es(R_es) = S_abs * ħω
  - E_rel_gs = E_gs(R_es) - E_gs(R_gs) = S_em  * ħω

Units:
  - Energies: eV throughout
  - Displacements: Å
  - Masses: amu
  - ΔQ: amu^(1/2)·Å  (internally converted to SI for S)
  - Frequencies: meV, cm⁻¹, THz

Dependencies:
  pip install numpy scipy matplotlib ase

Author: auto-generated for Marco's group
"""

import os
import re
import sys
import warnings
import numpy as np
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, Tuple, Dict, List

# Optional imports — graceful fallback
try:
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    warnings.warn("matplotlib not found — plotting disabled.")

try:
    from ase.io import read as ase_read
    from ase import units as ase_units
    HAS_ASE = True
except ImportError:
    HAS_ASE = False
    warnings.warn("ASE not found — using built-in parsers only.")

# ─────────────────────────────────────────────────────────────────────────────
# Physical constants (CODATA 2018)
# ─────────────────────────────────────────────────────────────────────────────
HBAR_eVs  = 6.582119569e-16   # ħ  [eV·s]
HBAR_Js   = 1.054571817e-34   # ħ  [J·s]
EV_TO_J   = 1.602176634e-19   # 1 eV in Joules
AMU_TO_KG = 1.66053906660e-27 # 1 amu in kg
ANG_TO_M  = 1e-10             # 1 Å in metres
EV_TO_CMi = 8065.544          # 1 eV in cm⁻¹
EV_TO_THZ = 241.7991          # 1 eV in THz (via ħω)
MEV_TO_EV = 1e-3


# ─────────────────────────────────────────────────────────────────────────────
# Data containers
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class Structure:
    """Minimal periodic structure."""
    symbols : List[str]           # atomic symbols
    masses  : np.ndarray          # [N] amu
    positions: np.ndarray         # [N,3] Å (Cartesian)
    lattice  : np.ndarray         # [3,3] Å
    energy   : Optional[float] = None  # eV


@dataclass
class PhononData:
    """Phonon eigensystem at Γ."""
    frequencies : np.ndarray   # [3N] eV (ħω)
    eigenvectors: np.ndarray   # [3N, 3N] mass-weighted, normalised
    # eigenvectors[mode, atom*3+cart]  (convention: e_k · e_k' = δ_kk')


@dataclass
class HRResult:
    """Full Huang-Rhys calculation output."""
    # ── Energy method ──────────────────────────────────────────────────────
    ZPL        : Optional[float] = None   # eV
    E_abs      : Optional[float] = None   # vertical absorption  eV
    E_em       : Optional[float] = None   # vertical emission     eV
    E_rel_es   : Optional[float] = None   # ES relaxation energy  eV
    E_rel_gs   : Optional[float] = None   # GS relaxation energy  eV
    S_abs      : Optional[float] = None   # HR factor (absorption side)
    S_em       : Optional[float] = None   # HR factor (emission side)
    omega_eff_eV: Optional[float] = None  # effective phonon energy eV

    # ── Displacement method ────────────────────────────────────────────────
    delta_Q    : Optional[float] = None   # amu^(1/2)·Å
    S_disp     : Optional[float] = None   # from ΔQ and ω_eff

    # ── Mode-resolved ──────────────────────────────────────────────────────
    S_modes    : Optional[np.ndarray] = None   # [3N] partial HR factors
    omega_modes: Optional[np.ndarray] = None   # [3N] eV

    # ── Derived spectroscopic quantities ──────────────────────────────────
    @property
    def S_total_modes(self):
        return float(np.sum(self.S_modes)) if self.S_modes is not None else None

    @property
    def omega_eff_from_modes(self):
        """⟨ω⟩ = Σ S_k ω_k / S_total  (HR-weighted mean frequency)."""
        if self.S_modes is None: return None
        S = self.S_modes; w = self.omega_modes
        return float(np.dot(S, w) / np.sum(S)) if np.sum(S) > 0 else None

    def summary(self) -> str:
        lines = ["=" * 58, "  Huang-Rhys Factor Summary", "=" * 58]

        def _fmt(label, val, unit=""):
            if val is None: return f"  {label:<30s}  {'N/A':>10s}"
            return f"  {label:<30s}  {val:>10.4f}  {unit}"

        if self.ZPL is not None:
            lines += [
                "",
                "  ── Configuration Coordinate Energies ──────────────",
                _fmt("ZPL",                 self.ZPL,       "eV"),
                _fmt("Vertical absorption", self.E_abs,     "eV"),
                _fmt("Vertical emission",   self.E_em,      "eV"),
                _fmt("ES relax. energy",    self.E_rel_es,  "eV"),
                _fmt("GS relax. energy",    self.E_rel_gs,  "eV"),
                _fmt("Stokes shift",
                     (self.E_rel_es + self.E_rel_gs)
                     if (self.E_rel_es and self.E_rel_gs) else None, "eV"),
            ]
        if self.omega_eff_eV is not None:
            w = self.omega_eff_eV
            lines += [
                "",
                "  ── Effective Phonon Energy ─────────────────────────",
                _fmt("ħω_eff", w, "eV"),
                _fmt("ħω_eff", w * 1e3, "meV"),
                _fmt("ħω_eff", w * EV_TO_CMi, "cm⁻¹"),
                _fmt("ħω_eff", w * EV_TO_THZ, "THz"),
            ]
        lines += [
            "",
            "  ── Huang-Rhys Factors ──────────────────────────────",
        ]
        if self.S_abs is not None:
            lines += [
                _fmt("S (energy, absorption)", self.S_abs),
                _fmt("S (energy, emission)",   self.S_em),
            ]
        if self.delta_Q is not None:
            lines += [_fmt("ΔQ", self.delta_Q, "amu^(1/2)·Å")]
        if self.S_disp is not None:
            lines += [_fmt("S (displacement)",  self.S_disp)]
        if self.S_modes is not None:
            lines += [
                _fmt("S (mode-resolved sum)", self.S_total_modes),
                _fmt("⟨ω⟩_HR", self.omega_eff_from_modes, "eV"),
            ]
        lines.append("=" * 58)
        return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# Parsers
# ─────────────────────────────────────────────────────────────────────────────

class VASPParser:
    """Parse VASP output files."""

    @staticmethod
    def read_structure(poscar_path: str, outcar_path: str = None) -> Structure:
        """
        Read structure from POSCAR/CONTCAR.
        Optionally read final energy from OUTCAR.
        """
        if HAS_ASE:
            atoms = ase_read(poscar_path, format="vasp")
            symbols   = list(atoms.get_chemical_symbols())
            masses    = atoms.get_masses()
            positions = atoms.get_positions()
            lattice   = atoms.get_cell().array
        else:
            symbols, masses, positions, lattice = VASPParser._parse_poscar(poscar_path)

        energy = None
        if outcar_path and os.path.isfile(outcar_path):
            energy = VASPParser._read_last_energy(outcar_path)

        return Structure(symbols=symbols, masses=masses,
                         positions=positions, lattice=lattice, energy=energy)

    @staticmethod
    def _parse_poscar(path: str):
        """Minimal POSCAR/CONTCAR parser (no ASE dependency)."""
        with open(path) as f:
            lines = f.readlines()
        scale   = float(lines[1])
        lattice = scale * np.array([[float(x) for x in l.split()] for l in lines[2:5]])
        species = lines[5].split()
        counts  = list(map(int, lines[6].split()))
        mode    = lines[7].strip()[0].lower()  # S(elective), D(irect), C(artesian)
        if mode == 's':
            mode = lines[8].strip()[0].lower()
            coord_start = 9
        else:
            coord_start = 8

        # Expand species list
        symbols = []
        for sp, n in zip(species, counts):
            symbols += [sp] * n

        # Atomic masses via simple lookup
        mass_table = {
            'H':1.008,'He':4.003,'Li':6.941,'Be':9.012,'B':10.811,'C':12.011,
            'N':14.007,'O':15.999,'F':18.998,'Ne':20.180,'Na':22.990,'Mg':24.305,
            'Al':26.982,'Si':28.086,'P':30.974,'S':32.065,'Cl':35.453,'Ar':39.948,
            'K':39.098,'Ca':40.078,'Sc':44.956,'Ti':47.867,'V':50.942,'Cr':51.996,
            'Mn':54.938,'Fe':55.845,'Co':58.933,'Ni':58.693,'Cu':63.546,'Zn':65.38,
            'Ga':69.723,'Ge':72.630,'As':74.922,'Se':78.971,'Br':79.904,'Kr':83.798,
            'Rb':85.468,'Sr':87.620,'Y':88.906,'Zr':91.224,'Nb':92.906,'Mo':95.960,
            'Tc':98.000,'Ru':101.07,'Rh':102.91,'Pd':106.42,'Ag':107.87,'Cd':112.41,
            'In':114.82,'Sn':118.71,'Sb':121.76,'Te':127.60,'I':126.90,'Xe':131.29,
            'Cs':132.91,'Ba':137.33,'La':138.91,'Ce':140.12,'Hf':178.49,'Ta':180.95,
            'W':183.84,'Re':186.21,'Os':190.23,'Ir':192.22,'Pt':195.08,'Au':196.97,
            'Hg':200.59,'Tl':204.38,'Pb':207.20,'Bi':208.98,'B':10.811,'N':14.007,
        }
        masses = np.array([mass_table.get(s, 1.0) for s in symbols])

        raw_coords = []
        for i in range(len(symbols)):
            vals = list(map(float, lines[coord_start + i].split()[:3]))
            raw_coords.append(vals)
        raw_coords = np.array(raw_coords)

        if mode == 'd':  # fractional → Cartesian
            positions = raw_coords @ lattice
        else:
            positions = raw_coords * scale

        return symbols, masses, positions, lattice

    @staticmethod
    def _read_last_energy(outcar_path: str) -> float:
        """Extract last 'energy(sigma->0)' from OUTCAR."""
        energy = None
        pattern = re.compile(r"energy\(sigma->0\)\s+=\s+([-\d.]+)")
        with open(outcar_path) as f:
            for line in f:
                m = pattern.search(line)
                if m:
                    energy = float(m.group(1))
        if energy is None:
            raise ValueError(f"No energy found in {outcar_path}")
        return energy

    @staticmethod
    def read_phonons(outcar_path: str) -> PhononData:
        """
        Parse Γ-point phonon eigensystem from OUTCAR.
        Requires IBRION = 5, 6, 7, or 8.
        """
        freqs, eigvecs = [], []
        with open(outcar_path) as f:
            content = f.read()

        # Regex to capture each phonon mode block
        mode_pattern = re.compile(
            r"(\d+) f(?:/i)?\s*=\s*([\d.]+) THz.*?meV\n((?:\s+.*\n)+?)\s*\n",
            re.MULTILINE
        )
        # Better: parse line-by-line
        lines = content.split("\n")
        i = 0
        while i < len(lines):
            line = lines[i]
            # Match frequency line
            m = re.match(
                r"\s*\d+\s+f(?:/i)?\s*=\s*([\d.]+)\s+THz\s+"
                r"([\d.]+)\s+2PiTHz\s+([\d.]+)\s+cm-1\s+([\d.]+)\s+meV",
                line
            )
            if m:
                freq_meV = float(m.group(4))
                # Check imaginary
                is_imag = "/i" in line
                freq_eV = (-1 if is_imag else 1) * freq_meV * MEV_TO_EV
                freqs.append(freq_eV)
                # Read eigenvector block: next lines with "X Y Z dx dy dz"
                eig = []
                i += 1
                while i < len(lines):
                    parts = lines[i].split()
                    if len(parts) == 6:
                        try:
                            dx, dy, dz = float(parts[3]), float(parts[4]), float(parts[5])
                            eig.extend([dx, dy, dz])
                            i += 1
                        except ValueError:
                            break
                    else:
                        break
                eigvecs.append(eig)
            else:
                i += 1

        if not freqs:
            raise ValueError(f"No phonon data found in {outcar_path}. "
                             f"Run VASP with IBRION=5/6/7/8.")

        freqs   = np.array(freqs)
        eigvecs = np.array(eigvecs)  # [3N, 3N]

        # VASP eigenvectors are NOT mass-weighted; they are Cartesian displacements
        # normalised per mode. We return them as-is; mass-weighting done in HR calc.
        return PhononData(frequencies=freqs, eigenvectors=eigvecs)


class QEParser:
    """Parse Quantum ESPRESSO output files."""

    @staticmethod
    def read_structure(output_path: str) -> Structure:
        """
        Parse pw.x output (scf or relax).
        Reads final geometry and energy.
        """
        if HAS_ASE:
            atoms = ase_read(output_path, format="espresso-out", index=-1)
            symbols   = list(atoms.get_chemical_symbols())
            masses    = atoms.get_masses()
            positions = atoms.get_positions()
            lattice   = atoms.get_cell().array
            energy    = atoms.get_potential_energy()  # eV
            return Structure(symbols=symbols, masses=masses,
                             positions=positions, lattice=lattice, energy=energy)
        else:
            return QEParser._parse_pw_output(output_path)

    @staticmethod
    def _parse_pw_output(path: str) -> Structure:
        """Fallback QE parser (no ASE)."""
        with open(path) as f:
            lines = f.readlines()

        # --- lattice parameter (alat in bohr) ---
        alat_bohr = None
        for l in lines:
            m = re.search(r"lattice parameter \(alat\)\s*=\s*([\d.]+)\s+a\.u\.", l)
            if m:
                alat_bohr = float(m.group(1))

        BOHR_TO_ANG = 0.529177210903

        # --- lattice vectors ---
        lattice = None
        for i, l in enumerate(lines):
            if "crystal axes:" in l or "a(1)" in l:
                try:
                    vecs = []
                    for j in range(3):
                        m = re.search(r"a\(\d\)\s*=\s*\(\s*([\d.\s-]+)\)", lines[i+j+1])
                        if m:
                            v = list(map(float, m.group(1).split()))
                            vecs.append(v)
                    if len(vecs) == 3:
                        lattice = np.array(vecs) * alat_bohr * BOHR_TO_ANG
                except Exception:
                    pass
            if lattice is not None:
                break

        # --- final atomic positions (angstrom) ---
        positions, symbols = [], []
        for i, l in enumerate(lines):
            if "ATOMIC_POSITIONS" in l and "angstrom" in l.lower():
                j = i + 1
                while j < len(lines) and lines[j].strip() and not lines[j].startswith("End"):
                    parts = lines[j].split()
                    if len(parts) >= 4:
                        try:
                            symbols.append(parts[0])
                            positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
                        except ValueError:
                            break
                    j += 1

        # --- final energy ---
        energy = None
        RY_TO_EV = 13.605693122994
        for l in lines:
            m = re.search(r"!\s+total energy\s*=\s*([-\d.]+)\s+Ry", l)
            if m:
                energy = float(m.group(1)) * RY_TO_EV

        mass_table = {  # abbreviated — extend as needed
            'H':1.008,'C':12.011,'N':14.007,'O':15.999,'F':18.998,'Si':28.086,
            'P':30.974,'S':32.065,'Cl':35.453,'Ga':69.723,'As':74.922,'Ge':72.630,
            'B':10.811,'Al':26.982,'In':114.82,'Sb':121.76,'Zn':65.38,'Se':78.971,
            'Hf':178.49,'W':183.84,'Mo':95.960,'Cr':51.996,'Mn':54.938,'Fe':55.845,
            'Co':58.933,'Ni':58.693,'Cu':63.546,'Mg':24.305,'Ca':40.078,'Ti':47.867,
        }
        masses = np.array([mass_table.get(s, 28.0) for s in symbols])

        return Structure(symbols=symbols, masses=masses,
                         positions=np.array(positions), lattice=lattice, energy=energy)

    @staticmethod
    def read_phonons(dynmat_path: str) -> PhononData:
        """
        Parse phonon data from ph.x dynmat output or matdyn.x output.
        Expects the standard QE dynmat.out format.
        """
        if HAS_ASE:
            # ASE can read QE phonon output for simple cases
            try:
                from ase.io.espresso import read_espresso_ph
                # This may not exist in all ASE versions; fall through to manual parse
            except ImportError:
                pass

        return QEParser._parse_dynmat(dynmat_path)

    @staticmethod
    def _parse_dynmat(path: str) -> PhononData:
        """Parse QE dynmat.out (output of dynmat.x or ph.x at q=0)."""
        with open(path) as f:
            lines = f.readlines()

        CM_TO_EV = 1.0 / EV_TO_CMi  # 1 cm⁻¹ in eV (ħω)

        freqs, eigvecs = [], []
        i = 0
        while i < len(lines):
            # Frequency line: "omega( N) = ... [cm-1]"
            m = re.match(r"\s*omega\(\s*\d+\)\s*=\s*([-\d.]+)\s+\[cm-1\]", lines[i])
            if m:
                freq_cm = float(m.group(1))
                freq_eV = freq_cm * CM_TO_EV  # negative → imaginary
                freqs.append(freq_eV)
                eig = []
                i += 1
                # Read eigenvector: lines like "atom N  e(x) e(y) e(z)"
                while i < len(lines):
                    parts = lines[i].split()
                    # Expected: "atom  N  x  y  z" or "(  x,  y,  z  )" patterns
                    if len(parts) >= 5 and parts[0] == "atom":
                        try:
                            ex, ey, ez = float(parts[2]), float(parts[3]), float(parts[4])
                            eig.extend([ex, ey, ez])
                            i += 1
                        except ValueError:
                            break
                    elif len(parts) == 3:
                        try:
                            eig.extend([float(x) for x in parts])
                            i += 1
                        except ValueError:
                            break
                    else:
                        break
                eigvecs.append(eig)
            else:
                i += 1

        if not freqs:
            raise ValueError(f"No phonon modes found in {path}.")

        return PhononData(frequencies=np.array(freqs),
                          eigenvectors=np.array(eigvecs))


# ─────────────────────────────────────────────────────────────────────────────
# Core Huang-Rhys Calculator
# ─────────────────────────────────────────────────────────────────────────────

class HuangRhysCalculator:
    """
    Compute Huang-Rhys factors for a defect/colour-centre system.

    Notation (Freysoldt et al., RMP 2014; Alkauskas et al., PRB 2014):
      GS = ground state  (state 1)
      ES = excited state (state 2)

    Four energies needed for energy method:
      E_gs(R_gs) : GS energy at GS geometry       [relaxed GS]
      E_gs(R_es) : GS energy at ES geometry       [GS single-point at ES geom]
      E_es(R_gs) : ES energy at GS geometry       [ES single-point at GS geom]
      E_es(R_es) : ES energy at ES geometry       [relaxed ES]
    """

    def __init__(self, verbose: bool = True):
        self.verbose = verbose

    # ── Public API ────────────────────────────────────────────────────────────

    def from_energies(
        self,
        E_gs_Rgs: float,  # eV
        E_gs_Res: float,  # eV
        E_es_Rgs: float,  # eV
        E_es_Res: float,  # eV
        omega_eff_eV: float,  # eV
    ) -> HRResult:
        """
        Energy method (single effective phonon mode approximation).

        Parameters
        ----------
        E_gs_Rgs : total energy of GS at GS geometry
        E_gs_Res : total energy of GS at ES geometry  (GS SP at ES coords)
        E_es_Rgs : total energy of ES at GS geometry  (ES SP at GS coords)
        E_es_Res : total energy of ES at ES geometry
        omega_eff_eV : effective phonon energy ħω_eff  [eV]
        """
        ZPL      = E_es_Res - E_gs_Rgs
        E_abs    = E_es_Rgs - E_gs_Rgs   # vertical absorption energy
        E_em     = E_es_Res - E_gs_Res   # vertical emission energy
        E_rel_es = E_es_Rgs - E_es_Res   # ES relaxation (Franck-Condon energy, abs)
        E_rel_gs = E_gs_Res - E_gs_Rgs   # GS relaxation (Franck-Condon energy, em)

        if E_rel_es < 0 or E_rel_gs < 0:
            warnings.warn(
                "Negative relaxation energy detected — check energy ordering. "
                "Ensure ES calculation used constrained or excited-state method."
            )

        S_abs = E_rel_es / omega_eff_eV
        S_em  = E_rel_gs / omega_eff_eV

        result = HRResult(
            ZPL=ZPL, E_abs=E_abs, E_em=E_em,
            E_rel_es=E_rel_es, E_rel_gs=E_rel_gs,
            S_abs=S_abs, S_em=S_em,
            omega_eff_eV=omega_eff_eV,
        )
        if self.verbose:
            print(result.summary())
        return result

    def from_structures(
        self,
        gs_struct: Structure,
        es_struct: Structure,
        omega_eff_eV: float,
        phonons: Optional[PhononData] = None,
    ) -> HRResult:
        """
        Displacement method (ΔQ) + optional mode-resolved.

        Computes ΔQ = sqrt(Σ_i m_i |Δr_i|²) and
        S_disp = ω_eff ΔQ² / (2ħ).

        If phonons provided, also computes mode-resolved S_k.
        """
        delta_Q = self._compute_delta_Q(gs_struct, es_struct)

        # Convert ΔQ from amu^(1/2)·Å to SI (kg^(1/2)·m)
        delta_Q_SI = delta_Q * np.sqrt(AMU_TO_KG) * ANG_TO_M

        # omega in rad/s from eV
        omega_SI = omega_eff_eV * EV_TO_J / HBAR_Js

        S_disp = 0.5 * omega_SI * delta_Q_SI**2 / HBAR_Js

        # Build result
        result = HRResult(
            omega_eff_eV=omega_eff_eV,
            delta_Q=delta_Q,
            S_disp=S_disp,
        )

        # Absorb energy information if available
        if gs_struct.energy is not None and es_struct.energy is not None:
            result.ZPL = es_struct.energy - gs_struct.energy

        # Mode-resolved
        if phonons is not None:
            S_modes = self._mode_resolved_HR(
                gs_struct, es_struct, phonons
            )
            result.S_modes    = S_modes
            result.omega_modes = phonons.frequencies

        if self.verbose:
            print(result.summary())
        return result

    def full_calculation(
        self,
        gs_struct : Structure,
        es_struct : Structure,
        E_gs_Rgs  : float,
        E_gs_Res  : float,
        E_es_Rgs  : float,
        E_es_Res  : float,
        omega_eff_eV : Optional[float] = None,
        phonons   : Optional[PhononData] = None,
    ) -> HRResult:
        """
        Combined energy + displacement + (optional) mode-resolved HR.

        If omega_eff_eV is not provided, it is estimated self-consistently from
        the energy method: ħω_eff ≈ sqrt(E_rel_es · E_rel_gs) (geometric mean).
        """
        E_rel_es = E_es_Rgs - E_es_Res
        E_rel_gs = E_gs_Res - E_gs_Rgs

        if omega_eff_eV is None:
            omega_eff_eV = np.sqrt(abs(E_rel_es * E_rel_gs))
            if self.verbose:
                print(f"  ω_eff estimated from √(E_rel_es·E_rel_gs) = "
                      f"{omega_eff_eV*1e3:.2f} meV")

        # Energy method
        result = self.from_energies(
            E_gs_Rgs, E_gs_Res, E_es_Rgs, E_es_Res, omega_eff_eV
        )

        # Displacement
        delta_Q   = self._compute_delta_Q(gs_struct, es_struct)
        dQ_SI     = delta_Q * np.sqrt(AMU_TO_KG) * ANG_TO_M
        omega_SI  = omega_eff_eV * EV_TO_J / HBAR_Js
        S_disp    = 0.5 * omega_SI * dQ_SI**2 / HBAR_Js

        result.delta_Q = delta_Q
        result.S_disp  = S_disp

        # Mode-resolved
        if phonons is not None:
            S_modes = self._mode_resolved_HR(gs_struct, es_struct, phonons)
            result.S_modes    = S_modes
            result.omega_modes = phonons.frequencies

        if self.verbose:
            print(f"\n  ΔQ = {delta_Q:.4f}  amu^(1/2)·Å")
            print(f"  S (displacement) = {S_disp:.4f}")
            if phonons:
                print(f"  S (mode-resolved sum) = {result.S_total_modes:.4f}")
        return result

    # ── Private methods ───────────────────────────────────────────────────────

    @staticmethod
    def _compute_delta_Q(gs: Structure, es: Structure) -> float:
        """
        Generalised displacement coordinate:
          ΔQ = sqrt(Σ_i m_i |Δr_i|²)   [amu^(1/2)·Å]

        Wraps periodic images to minimise displacement.
        """
        assert len(gs.masses) == len(es.masses), \
            "Ground and excited state structures must have the same atoms."

        delta_r = es.positions - gs.positions

        # Minimum-image convention (orthorhombic approximation)
        # For non-orthorhombic cells, use fractional coordinates
        if gs.lattice is not None:
            inv_lat = np.linalg.inv(gs.lattice)
            delta_f = delta_r @ inv_lat
            delta_f -= np.round(delta_f)
            delta_r  = delta_f @ gs.lattice

        delta_Q_sq = np.sum(gs.masses[:, None] * delta_r**2)
        return float(np.sqrt(delta_Q_sq))

    @staticmethod
    def _mode_resolved_HR(
        gs: Structure,
        es: Structure,
        phonons: PhononData,
    ) -> np.ndarray:
        """
        Partial HR factor for each phonon mode k:
          S_k = ω_k / (2ħ) · |⟨Δ̃R | e_k⟩|²

        where Δ̃R_iα = sqrt(m_i) Δr_iα  is the mass-weighted displacement
        and e_k are mass-weighted normalised eigenvectors from phonon calc.

        VASP OUTCAR gives Cartesian eigenvectors u_k (NOT mass-weighted).
        We mass-weight: ẽ_k,iα = sqrt(m_i) u_k,iα  then renormalise.
        """
        N = len(gs.masses)

        # Mass-weighted displacement vector [3N]
        delta_r = es.positions - gs.positions
        if gs.lattice is not None:
            inv_lat = np.linalg.inv(gs.lattice)
            delta_f = delta_r @ inv_lat
            delta_f -= np.round(delta_f)
            delta_r  = delta_f @ gs.lattice

        sqrt_m = np.repeat(np.sqrt(gs.masses), 3)   # [3N]
        Delta  = sqrt_m * delta_r.flatten()           # mass-weighted ΔR [3N]  [amu^(1/2)·Å]

        S_modes = np.zeros(len(phonons.frequencies))
        for k, (omega_k, evec_k) in enumerate(
            zip(phonons.frequencies, phonons.eigenvectors)
        ):
            if omega_k <= 0:
                continue  # skip imaginary/acoustic modes

            # Mass-weight and normalise eigenvector
            evec = np.array(evec_k[:3*N])
            evec_mw = sqrt_m * evec
            norm = np.linalg.norm(evec_mw)
            if norm < 1e-12:
                continue
            evec_mw /= norm

            # Projection [amu^(1/2)·Å]
            proj = float(np.dot(Delta, evec_mw))

            # Convert to SI
            proj_SI    = proj * np.sqrt(AMU_TO_KG) * ANG_TO_M
            omega_SI   = omega_k * EV_TO_J / HBAR_Js
            S_modes[k] = 0.5 * omega_SI * proj_SI**2 / HBAR_Js

        return S_modes


# ─────────────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────────────

class HRPlotter:
    """Visualisation utilities."""

    @staticmethod
    def configuration_coordinate_diagram(
        result: HRResult,
        ax=None,
        n_points: int = 300,
        title: str = "Configuration Coordinate Diagram",
    ):
        """
        Plot parabolic GS and ES potential energy surfaces.
        Requires result from from_energies() or full_calculation().
        """
        if not HAS_MPL:
            warnings.warn("matplotlib not available.")
            return

        if result.omega_eff_eV is None or result.ZPL is None:
            warnings.warn("Need energy-method result for CC diagram.")
            return

        show = ax is None
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 5))

        # Parabola curvature:  E(Q) = ½ ω² Q²  (in eV, Q in amu^(1/2)·Å)
        # ω in rad/s → convert curvature to eV/(amu·Å²)
        omega_eV = result.omega_eff_eV
        omega_SI = omega_eV * EV_TO_J / HBAR_Js
        # k = mω² → for generalised Q: curvature = ω²/(2) in SI
        # E [eV] = (ω_SI²/2) * Q_SI² / EV_TO_J
        #        = (ω_eV² / (2ħ²)) * Q_SI² * (HBAR_Js² / EV_TO_J)
        # Simpler: set Q₀ = ΔQ, then E_rel_es = ½ k ΔQ² → k = 2E_rel_es/ΔQ²
        if result.delta_Q is not None and result.delta_Q > 0:
            dQ = result.delta_Q
        else:
            # Recover ΔQ from S_disp and ω
            dQ = np.sqrt(2 * result.S_abs * HBAR_Js /
                         (omega_SI * AMU_TO_KG * ANG_TO_M**2)) if result.S_abs else 1.0

        Q_min = -0.3 * dQ
        Q_max =  1.5 * dQ
        Q = np.linspace(Q_min, Q_max, n_points)

        # GS parabola centred at Q=0, E=0
        k = 2 * result.E_rel_gs / dQ**2 if result.E_rel_gs else omega_SI**2 / 2
        E_gs = 0.5 * k * Q**2

        # ES parabola centred at Q=ΔQ, E=ZPL
        E_es = result.ZPL + 0.5 * k * (Q - dQ)**2

        ax.plot(Q, E_gs, color="#2C6FAC", lw=2.5, label="Ground state")
        ax.plot(Q, E_es, color="#C0392B", lw=2.5, label="Excited state")

        # Transition arrows
        ax.annotate("", xy=(0, result.E_abs), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="->", color="#27AE60", lw=1.8))
        ax.annotate("", xy=(dQ, result.E_em), xytext=(dQ, result.ZPL),
                    arrowprops=dict(arrowstyle="->", color="#E67E22", lw=1.8))
        ax.axhline(result.ZPL,  color="gray", ls="--", lw=1, alpha=0.7)
        ax.axhline(result.E_abs, color="#27AE60", ls=":", lw=1, alpha=0.7)

        ax.text(Q_min + 0.05*(Q_max-Q_min), result.E_abs + 0.02,
                f"$E_{{abs}}$ = {result.E_abs:.3f} eV", color="#27AE60", fontsize=9)
        ax.text(dQ + 0.03*dQ, result.E_em - 0.08,
                f"$E_{{em}}$ = {result.E_em:.3f} eV", color="#E67E22", fontsize=9)
        ax.text(Q_min + 0.05*(Q_max-Q_min), result.ZPL - 0.08,
                f"ZPL = {result.ZPL:.3f} eV", color="gray", fontsize=9)

        ax.set_xlabel(r"Configuration coordinate $Q$  (amu$^{1/2}$·Å)", fontsize=11)
        ax.set_ylabel("Energy (eV)", fontsize=11)
        ax.set_title(title, fontsize=12)
        ax.legend(fontsize=10)
        ax.set_xlim(Q_min, Q_max)
        plt.tight_layout()
        if show:
            plt.show()
        return ax

    @staticmethod
    def spectral_function(
        result: HRResult,
        sigma_meV: float = 5.0,
        ax=None,
        title: str = "Huang-Rhys Spectral Function",
    ):
        """
        Plot A(ω) = Σ_k S_k δ(ω - ω_k) convolved with Gaussian.
        Requires mode-resolved result.
        """
        if not HAS_MPL:
            return
        if result.S_modes is None:
            warnings.warn("Mode-resolved HR factors required for spectral function.")
            return

        show = ax is None
        if ax is None:
            fig, ax = plt.subplots(figsize=(7, 4))

        sigma = sigma_meV * MEV_TO_EV
        omega = result.omega_modes
        S     = result.S_modes

        mask = (omega > 0) & (S > 1e-6)
        w_range = np.linspace(0, omega[mask].max() * 1.1, 1000)
        A = np.zeros_like(w_range)
        for wk, sk in zip(omega[mask], S[mask]):
            A += sk * np.exp(-0.5 * ((w_range - wk) / sigma)**2) \
                    / (sigma * np.sqrt(2 * np.pi))

        ax.fill_between(w_range * 1e3, A, alpha=0.4, color="#2C6FAC")
        ax.plot(w_range * 1e3, A, color="#2C6FAC", lw=1.5)
        ax.vlines(omega[mask] * 1e3, 0, S[mask], color="#C0392B",
                  alpha=0.7, linewidth=0.8, label=r"$S_k$")

        ax.set_xlabel(r"Phonon energy $\hbar\omega$ (meV)", fontsize=11)
        ax.set_ylabel(r"$A(\omega)$ = $\Sigma_k S_k\,\delta(\omega - \omega_k)$", fontsize=11)
        ax.set_title(title, fontsize=12)
        ax.legend(fontsize=10)
        plt.tight_layout()
        if show:
            plt.show()
        return ax

    @staticmethod
    def partial_HR_bar(result: HRResult, top_n: int = 30, ax=None):
        """Bar chart of largest partial HR factors."""
        if not HAS_MPL:
            return
        if result.S_modes is None:
            return

        show = ax is None
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 4))

        mask = (result.omega_modes > 0)
        w = result.omega_modes[mask] * 1e3   # meV
        S = result.S_modes[mask]
        idx = np.argsort(S)[::-1][:top_n]

        ax.bar(w[idx], S[idx], width=1.5, color="#2C6FAC", alpha=0.8)
        ax.set_xlabel(r"$\hbar\omega$ (meV)", fontsize=11)
        ax.set_ylabel(r"$S_k$", fontsize=11)
        ax.set_title(f"Top-{top_n} partial Huang-Rhys factors", fontsize=12)
        plt.tight_layout()
        if show:
            plt.show()
        return ax


# ─────────────────────────────────────────────────────────────────────────────
# Convenience wrappers
# ─────────────────────────────────────────────────────────────────────────────

def hr_from_vasp(
    gs_contcar   : str,
    es_contcar   : str,
    gs_outcar    : str,
    es_outcar    : str,
    gs_sp_outcar : str = None,  # GS energy computed at ES geometry
    es_sp_outcar : str = None,  # ES energy computed at GS geometry
    phonon_outcar: str = None,
    omega_eff_eV : float = None,
) -> HRResult:
    """
    One-call interface for VASP.

    Minimum required: gs_contcar, es_contcar, gs_outcar, es_outcar.
    For energy method: also gs_sp_outcar, es_sp_outcar.
    For mode-resolved: also phonon_outcar.

    Single-point energies:
      gs_sp_outcar  : GS functional, evaluated at ES geometry (static VASP run)
      es_sp_outcar  : ES functional, evaluated at GS geometry (static VASP run)
    """
    parser = VASPParser()
    calc   = HuangRhysCalculator()

    gs = parser.read_structure(gs_contcar, gs_outcar)
    es = parser.read_structure(es_contcar, es_outcar)

    phonons = None
    if phonon_outcar:
        phonons = parser.read_phonons(phonon_outcar)

    # ── Energy method ───────────────────────────────────────────────────────
    if gs_sp_outcar and es_sp_outcar:
        E_gs_Rgs = gs.energy
        E_es_Res = es.energy
        E_gs_Res = parser._read_last_energy(gs_sp_outcar)
        E_es_Rgs = parser._read_last_energy(es_sp_outcar)
        return calc.full_calculation(
            gs, es,
            E_gs_Rgs, E_gs_Res, E_es_Rgs, E_es_Res,
            omega_eff_eV=omega_eff_eV,
            phonons=phonons,
        )
    else:
        # Displacement method only
        if omega_eff_eV is None:
            raise ValueError("Provide omega_eff_eV if single-point energies are not available.")
        return calc.from_structures(gs, es, omega_eff_eV, phonons=phonons)


def hr_from_qe(
    gs_output    : str,
    es_output    : str,
    gs_sp_output : str = None,
    es_sp_output : str = None,
    dynmat_output: str = None,
    omega_eff_eV : float = None,
) -> HRResult:
    """
    One-call interface for Quantum ESPRESSO.

    Parameters mirror hr_from_vasp but for QE pw.x / dynmat.x outputs.
    """
    parser = QEParser()
    calc   = HuangRhysCalculator()

    gs = parser.read_structure(gs_output)
    es = parser.read_structure(es_output)

    phonons = None
    if dynmat_output:
        phonons = parser.read_phonons(dynmat_output)

    if gs_sp_output and es_sp_output:
        E_gs_Rgs = gs.energy
        E_es_Res = es.energy
        E_gs_Res = parser.read_structure(gs_sp_output).energy
        E_es_Rgs = parser.read_structure(es_sp_output).energy
        return calc.full_calculation(
            gs, es,
            E_gs_Rgs, E_gs_Res, E_es_Rgs, E_es_Res,
            omega_eff_eV=omega_eff_eV,
            phonons=phonons,
        )
    else:
        if omega_eff_eV is None:
            raise ValueError("Provide omega_eff_eV if single-point energies not available.")
        return calc.from_structures(gs, es, omega_eff_eV, phonons=phonons)


# ─────────────────────────────────────────────────────────────────────────────
# Demo / self-test
# ─────────────────────────────────────────────────────────────────────────────

def _demo():
    """
    Synthetic demo: NV-centre-like defect in a 215-atom supercell.
    All numbers are illustrative; replace with your VASP/QE outputs.
    """
    print("\n" + "=" * 58)
    print("  DEMO: Synthetic NV-centre-like system")
    print("=" * 58)

    # ── Fake structures ──────────────────────────────────────────────────────
    np.random.seed(42)
    N = 64   # atoms
    symbols  = ["C"] * (N-1) + ["N"]
    masses   = np.array([12.011] * (N-1) + [14.007])
    lattice  = np.eye(3) * 14.0   # 14 Å cubic supercell

    # GS: perfect lattice
    frac_gs = np.random.rand(N, 3)
    pos_gs  = frac_gs @ lattice

    # ES: small distortion around defect site (atom 0)
    pos_es = pos_gs.copy()
    displace = np.zeros((N, 3))
    displace[:5] = np.random.randn(5, 3) * 0.05   # Å displacement near defect

    pos_es += displace

    gs = Structure(symbols=symbols, masses=masses, positions=pos_gs,
                   lattice=lattice, energy=-3450.123)
    es = Structure(symbols=symbols, masses=masses, positions=pos_es,
                   lattice=lattice, energy=-3447.891)

    # ── Four energies (CC diagram) ──────────────────────────────────────────
    omega_eff = 0.065  # eV (~524 cm⁻¹ — typical diamond phonon)
    E_gs_Rgs  = -3450.123
    E_es_Res  = -3447.891
    E_rel_es  = 0.18   # ES relaxation energy
    E_rel_gs  = 0.16   # GS relaxation energy
    E_es_Rgs  = E_es_Res + E_rel_es
    E_gs_Res  = E_gs_Rgs + E_rel_gs

    calc   = HuangRhysCalculator(verbose=False)
    result = calc.full_calculation(
        gs, es,
        E_gs_Rgs, E_gs_Res, E_es_Rgs, E_es_Res,
        omega_eff_eV=omega_eff,
    )

    # ── Fake mode-resolved for spectral function demo ───────────────────────
    n_modes = 3 * N
    freqs   = np.sort(np.abs(np.random.randn(n_modes)) * 0.04 + 0.03)
    freqs[:6] = 1e-5   # acoustic
    evecs   = np.random.randn(n_modes, n_modes)
    for i in range(n_modes):  # rough orthonormalisation
        evecs[i] /= np.linalg.norm(evecs[i]) + 1e-12
    phonons = PhononData(frequencies=freqs, eigenvectors=evecs)

    result2 = calc.from_structures(gs, es, omega_eff, phonons=phonons)

    print(result.summary())
    print(f"\n  [mode-resolved] S_total = {result2.S_total_modes:.4f}")
    print(f"  [mode-resolved] ⟨ω⟩_HR  = "
          f"{result2.omega_eff_from_modes*1e3:.2f} meV")

    # ── Plots ────────────────────────────────────────────────────────────────
    if HAS_MPL:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        fig.suptitle("Huang-Rhys Demo (synthetic NV-centre-like)", fontsize=13)
        result.delta_Q = result2.delta_Q      # copy ΔQ for CC diagram
        HRPlotter.configuration_coordinate_diagram(result, ax=axes[0])
        HRPlotter.spectral_function(result2, ax=axes[1])
        HRPlotter.partial_HR_bar(result2, top_n=20, ax=axes[2])
        plt.tight_layout()
        out = "/mnt/user-data/outputs/hr_demo.png"
        plt.savefig(out, dpi=150)
        print(f"\n  Plot saved → {out}")
        plt.close()


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def _cli():
    import argparse

    parser = argparse.ArgumentParser(
        description="Huang-Rhys factor from VASP / QE calculations.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples (VASP):
  python huang_rhys.py vasp \\
      --gs-poscar gs/CONTCAR --es-poscar es/CONTCAR \\
      --gs-outcar gs/OUTCAR  --es-outcar es/OUTCAR  \\
      --gs-sp-outcar gs_sp/OUTCAR --es-sp-outcar es_sp/OUTCAR \\
      --phonon-outcar phonons/OUTCAR

Examples (QE):
  python huang_rhys.py qe \\
      --gs-output gs.relax.out --es-output es.relax.out \\
      --gs-sp-output gs_sp.scf.out --es-sp-output es_sp.scf.out \\
      --dynmat dynmat.out --omega 0.065
        """,
    )
    sub = parser.add_subparsers(dest="code")

    # ── VASP subcommand ──────────────────────────────────────────────────────
    v = sub.add_parser("vasp")
    v.add_argument("--gs-poscar",     required=True)
    v.add_argument("--es-poscar",     required=True)
    v.add_argument("--gs-outcar",     required=True)
    v.add_argument("--es-outcar",     required=True)
    v.add_argument("--gs-sp-outcar",  default=None)
    v.add_argument("--es-sp-outcar",  default=None)
    v.add_argument("--phonon-outcar", default=None)
    v.add_argument("--omega",         type=float, default=None,
                   help="Effective phonon energy ħω_eff  [eV]")
    v.add_argument("--plot",          action="store_true")

    # ── QE subcommand ────────────────────────────────────────────────────────
    q = sub.add_parser("qe")
    q.add_argument("--gs-output",     required=True)
    q.add_argument("--es-output",     required=True)
    q.add_argument("--gs-sp-output",  default=None)
    q.add_argument("--es-sp-output",  default=None)
    q.add_argument("--dynmat",        default=None)
    q.add_argument("--omega",         type=float, default=None)
    q.add_argument("--plot",          action="store_true")

    # ── demo ─────────────────────────────────────────────────────────────────
    sub.add_parser("demo")

    args = parser.parse_args()

    if args.code == "demo" or args.code is None:
        _demo()
        return

    if args.code == "vasp":
        result = hr_from_vasp(
            gs_contcar    = args.gs_poscar,
            es_contcar    = args.es_poscar,
            gs_outcar     = args.gs_outcar,
            es_outcar     = args.es_outcar,
            gs_sp_outcar  = args.gs_sp_outcar,
            es_sp_outcar  = args.es_sp_outcar,
            phonon_outcar = args.phonon_outcar,
            omega_eff_eV  = args.omega,
        )
    elif args.code == "qe":
        result = hr_from_qe(
            gs_output     = args.gs_output,
            es_output     = args.es_output,
            gs_sp_output  = args.gs_sp_output,
            es_sp_output  = args.es_sp_output,
            dynmat_output = args.dynmat,
            omega_eff_eV  = args.omega,
        )

    print(result.summary())
    if args.plot and HAS_MPL:
        HRPlotter.configuration_coordinate_diagram(result)
        if result.S_modes is not None:
            HRPlotter.spectral_function(result)


if __name__ == "__main__":
    _cli()
