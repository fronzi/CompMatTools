#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GenerateBilayer3.2.py
- Finds commensurate (ZSL) or near-commensurate (fallback) supercells for WO3 // graphene.
- Writes POSCAR_graphene_SC / POSCAR_WO3_SC (+ optional bilayer).
- Logs a CSV row per run/candidate/pick into ./bilayers_summary.csv

Tip: we disable numba by default to avoid llvmlite/libffi issues.
"""

import os
import sys
import csv
import json
from datetime import datetime
import numpy as np
from pymatgen.core import Structure, Lattice

# --- Robust imports / env ---
os.environ.setdefault("PMG_USE_NUMBA", "0")
os.environ.setdefault("NUMBA_DISABLE_JIT", "1")

try:
    from pymatgen.analysis.interfaces.zsl import ZSLGenerator
except Exception:
    print("ERROR: Could not import ZSLGenerator (pymatgen). Install/update from conda-forge.")
    sys.exit(1)

try:
    from ase.io import read as ase_read, write as ase_write
    ASE_OK = True
except Exception:
    ASE_OK = False

# ======= INPUT =======
POSCAR_WO3      = "POSCAR_WO3"
POSCAR_GRAPHENE = "POSCAR_graphene"

# ======= ZSL search parameters (relaxed) =======
MAX_MISMATCH   = 0.05
MAX_ANGLE_DIFF = 3.0
MAX_AREA       = 1200.0
ROTATE_SUBSTR  = True
ROT_STEP_DEG   = 0.5

# ======= Near-commensurate fallback =======
FALLBACK_MAX_MULTI_FILM = 18
FALLBACK_MAX_MULTI_SUB  = 12
FALLBACK_MAX_CELLS      = 3000

# Auto-pick thresholds for fallback
AUTOPICK_MAX_ABS_STRAIN = 2.5    # % (max |εa|, |εb|)
AUTOPICK_MAX_AREA       = 900.0  # Å^2

# ======= Quick bilayer output =======
MAKE_BILAYER = True
GAP_A        = 3.3
VACUUM_C     = 20.0

# ======= CSV/logging controls =======
ALWAYS_CREATE_CSV     = True
LOG_ALL_FALLBACK_CAND = True
RUN_ID_TAG            = None  # e.g. "WO3_graphene_test1"

CSV_FILE = "bilayers_summary.csv"
CSV_COLUMNS = [
    "timestamp","run_id","mode","status",
    "score","area","angle_diff_deg",
    "film_strain_a_pct","film_strain_b_pct","film_dgamma_deg",
    "Mf","Ms"
]

def _now():
    return datetime.now().isoformat(timespec="seconds")

def ensure_csv_header(filename: str = CSV_FILE):
    if not os.path.isfile(filename):
        with open(filename, "w", newline="") as f:
            csv.DictWriter(f, fieldnames=CSV_COLUMNS).writeheader()

def append_to_csv(row: dict, filename: str = CSV_FILE):
    out = {k: row.get(k, "") for k in CSV_COLUMNS}
    out["timestamp"] = out["timestamp"] or _now()
    out["run_id"] = out["run_id"] or RUN_ID_TAG or ""
    for key in ("Mf","Ms"):
        if isinstance(out.get(key), (list, tuple, np.ndarray)):
            out[key] = json.dumps(np.array(out[key]).tolist())
    if not os.path.isfile(filename):
        ensure_csv_header(filename)
    with open(filename, "a", newline="") as f:
        csv.DictWriter(f, fieldnames=CSV_COLUMNS).writerow(out)

# ======= Utilities =======
def ab_gamma(latt2x3):
    a, b, _ = latt2x3
    aL, bL = np.linalg.norm(a), np.linalg.norm(b)
    cosg = np.dot(a, b) / (aL * bL)
    cosg = np.clip(cosg, -1, 1)
    gamma = np.degrees(np.arccos(cosg))
    return aL, bL, gamma

def to_3x3(supercell_2x2):
    M = np.eye(3, dtype=int)
    M[:2, :2] = np.array(supercell_2x2, dtype=int)
    return M

def build_supercell(struct: Structure, mat):
    M = to_3x3(mat) if np.array(mat).shape == (2, 2) else np.array(mat, dtype=int)
    sc = struct.copy()
    sc.make_supercell(M)  # in-place
    return sc

def replace_ab(struct: Structure, new_a, new_b):
    old = struct.lattice.matrix
    new_cell = np.vstack([new_a, new_b, old[2]])
    return Structure(Lattice(new_cell), struct.species, struct.frac_coords)

def area_ab(cell):
    a, b, _ = cell
    return np.linalg.norm(np.cross(a, b))

def print_top_matches(scored, n=3):
    k = min(n, len(scored))
    print(f"Found {len(scored)} ZSL matches. Top {k}:")
    for i, (score, area, angle_diff, m) in enumerate(scored[:k], 1):
        print(f"  #{i}: score≈{score:.3f}, area≈{area:.1f} Å², Δγ≈{angle_diff:.2f}°")

# ======= ZSL wrapper =======
def zsl_search(film: Structure, sub: Structure):
    use_init_params = True
    try:
        gen = ZSLGenerator(
            max_area=MAX_AREA,
            max_mismatch=MAX_MISMATCH,
            max_angle_diff=MAX_ANGLE_DIFF,
            rotation_step=(ROT_STEP_DEG if ROTATE_SUBSTR else 0.0),
        )
    except TypeError:
        gen = ZSLGenerator()
        use_init_params = False

    matches = None
    for attr in ["generate", "get_matches", "get_interfaces"]:
        if hasattr(gen, attr):
            try:
                if use_init_params:
                    matches = getattr(gen, attr)(film, sub)
                else:
                    matches = getattr(gen, attr)(
                        film, sub,
                        max_area=MAX_AREA,
                        max_mismatch=MAX_MISMATCH,
                        max_angle_diff=MAX_ANGLE_DIFF,
                        rotation_step=(ROT_STEP_DEG if ROTATE_SUBSTR else 0.0),
                    )
                if matches:
                    break
            except TypeError:
                try:
                    if use_init_params:
                        matches = getattr(gen, attr)(substrate=sub, film=film)
                    else:
                        matches = getattr(gen, attr)(
                            substrate=sub, film=film,
                            max_area=MAX_AREA,
                            max_mismatch=MAX_MISMATCH,
                            max_angle_diff=MAX_ANGLE_DIFF,
                            rotation_step=(ROT_STEP_DEG if ROTATE_SUBSTR else 0.0),
                        )
                    if matches:
                        break
                except Exception:
                    continue
    return matches or []

def score_matches(matches):
    scored = []
    for m in matches:
        film_scM = m.get("film_supercell_matrix", m.get("film_sl_vecs"))
        sub_scM  = m.get("substrate_supercell_matrix", m.get("substrate_sl_vecs"))
        film_strain = m.get("film_strain", None)
        sub_strain  = m.get("substrate_strain", None)
        area        = float(m.get("match_area", m.get("area", -1.0)))
        angle_diff  = float(m.get("angle", m.get("angle_diff", 0.0)))
        def s_mean(x):
            if x is None: return 999.0
            x = np.array(x, dtype=float)
            return float(np.mean(np.abs(x)))
        score = 0.5 * (s_mean(film_strain) + s_mean(sub_strain))
        scored.append((score, area, angle_diff, m))
    scored.sort(key=lambda t: (t[0], t[1]))
    return scored

# ======= Near-commensurate fallback =======
def best_near_commensurate(film_str: Structure, sub_str: Structure,
                           max_multi_film=18, max_multi_sub=12,
                           max_cells=3000):
    def mats(limit):
        for h11 in range(1, limit + 1):
            for h22 in range(1, limit + 1):
                for h12 in range(-min(h11, h22), min(h11, h22) + 1):
                    yield np.array([[h11, h12],
                                    [0,   h22]], dtype=int)

    def embed2(M2):
        M = np.eye(3, dtype=int); M[:2,:2] = M2; return M

    def make_sc(S, M2):
        sc = S.copy(); sc.make_supercell(embed2(M2)); return sc

    def abg(cell):
        a,b,_ = cell
        aL,bL = np.linalg.norm(a), np.linalg.norm(b)
        cosg = np.dot(a,b)/(aL*bL); cosg = np.clip(cosg, -1, 1)
        g = np.degrees(np.arccos(cosg))
        return aL,bL,g

    cand, count = [], 0
    for Mf in mats(max_multi_film):
        for Ms in mats(max_multi_sub):
            count += 1
            if count > max_cells:
                break
            F = make_sc(film_str, Mf)
            S = make_sc(sub_str,  Ms)
            newF = replace_ab(F, S.lattice.matrix[0], S.lattice.matrix[1])
            aF0,bF0,gF0 = abg(F.lattice.matrix)
            aF1,bF1,gF1 = abg(newF.lattice.matrix)
            ea = (aF1 - aF0) / aF0 * 100.0
            eb = (bF1 - bF0) / bF0 * 100.0
            dga = (gF1 - gF0)
            score = 0.67 * (abs(ea) + abs(eb)) / 2.0 + 0.33 * abs(dga)
            areaS = area_ab(S.lattice.matrix)
            cand.append(dict(score=score, area=areaS, Mf=Mf, Ms=Ms,
                             film_strain=(ea, eb), angle_change=dga))
        if count > max_cells:
            break
    cand.sort(key=lambda d: (d["score"], d["area"]))
    return cand[:10]

def write_quick_bilayer():
    if not ASE_OK:
        print("[INFO] ASE not available: skipping quick bilayer write.")
        return
    g = ase_read("POSCAR_graphene_SC", format="vasp")
    w = ase_read("POSCAR_WO3_SC",      format="vasp")
    g.positions[:, 2] -= g.positions[:, 2].min()
    w.positions[:, 2] -= w.positions[:, 2].min()
    t_g = float(g.positions[:, 2].max())
    t_w = float(w.positions[:, 2].max())
    cell_final = w.cell.array.copy()
    cell_final[2] = [0.0, 0.0, VACUUM_C]
    z0 = (VACUUM_C - (t_w + GAP_A + t_g)) * 0.5
    w.positions[:, 2] += z0
    g.positions[:, 2] += z0 + t_w + GAP_A
    w.set_cell(cell_final, scale_atoms=False)
    g.set_cell(cell_final, scale_atoms=False)
    bil = w.copy(); bil.extend(g); bil.set_pbc([True, True, False])
    ase_write("POSCAR_bilayer_check", bil, format="vasp")
    print(f"[OK] Quick bilayer → POSCAR_bilayer_check (gap≈{GAP_A} Å, c≈{VACUUM_C} Å)")

# ======= Main =======
def main():
    if ALWAYS_CREATE_CSV:
        ensure_csv_header()
        append_to_csv({"mode": "RUN", "status": "start"})

    sub  = Structure.from_file(POSCAR_WO3)
    film = Structure.from_file(POSCAR_GRAPHENE)

    # --- 1) ZSL search ---
    matches = zsl_search(film, sub)
    if matches:
        scored = score_matches(matches)
        print_top_matches(scored, n=3)
        best = scored[0][3]
        film_scM = np.array(best.get("film_supercell_matrix", best.get("film_sl_vecs")))
        sub_scM  = np.array(best.get("substrate_supercell_matrix", best.get("substrate_sl_vecs")))
        film_sc = build_supercell(film, film_scM)
        sub_sc  = build_supercell(sub,  sub_scM)
        sub_cell  = sub_sc.lattice.matrix.copy()
        film_cell = film_sc.lattice.matrix.copy()
        film_sc   = replace_ab(film_sc, sub_cell[0], sub_cell[1])
        aF0,bF0,gF0 = ab_gamma(film_cell)
        aF1,bF1,gF1 = ab_gamma(film_sc.lattice.matrix)
        aS0,bS0,gS0 = ab_gamma(sub_sc.lattice.matrix)
        print("\n=== STRAIN (film = Graphene, after align) ===")
        print(f"  a: {(aF1-aF0)/aF0*100:+.2f}%   b: {(bF1-bF0)/bF0*100:+.2f}%   Δγ: {gF1-gF0:+.2f}°")
        print("=== SUBSTRATE (WO3) metrics (unchanged) ===")
        print(f"  a: {aS0:.4f} Å   b: {bS0:.4f} Å   γ: {gS0:.2f}°")
        film_sc.to(fmt="poscar", filename="POSCAR_graphene_SC")
        sub_sc.to(fmt="poscar",  filename="POSCAR_WO3_SC")
        print("\n[OK] Wrote POSCAR_graphene_SC and POSCAR_WO3_SC (ZSL best match)")
        append_to_csv({
            "mode": "ZSL", "status": "picked",
            "score": scored[0][0], "area": scored[0][1], "angle_diff_deg": scored[0][2],
            "film_strain_a_pct": (aF1-aF0)/aF0*100.0,
            "film_strain_b_pct": (bF1-bF0)/bF0*100.0,
            "film_dgamma_deg": (gF1-gF0),
            "Mf": film_scM.tolist(), "Ms": sub_scM.tolist(),
        })
        if MAKE_BILAYER:
            write_quick_bilayer()
        append_to_csv({"mode": "RUN", "status": "end"})
        return

    # --- 2) Fallback: near-commensurate ---
    print("[INFO] ZSL found no matches under current limits.")
    print("[INFO] Running near-commensurate brute-force search…")
    near = best_near_commensurate(
        film, sub,
        max_multi_film=FALLBACK_MAX_MULTI_FILM,
        max_multi_sub=FALLBACK_MAX_MULTI_SUB,
        max_cells=FALLBACK_MAX_CELLS,
    )
    if not near:
        append_to_csv({"mode": "Fallback", "status": "no_candidates"})
        raise RuntimeError("No near-commensurate candidates found. Increase search limits.")

    print("Top near-commensurate candidates:")
    for i, d in enumerate(near, 1):
        ea, eb = d["film_strain"]
        print(f" #{i}: score={d['score']:.3f}, area~{d['area']:.1f} Å², "
              f"film εa={ea:+.2f}%, εb={eb:+.2f}%, Δγ={d['angle_change']:+.2f}°")
        print(f"    Mf (film):\n{d['Mf']}\n    Ms (substrate):\n{d['Ms']}\n")
        if LOG_ALL_FALLBACK_CAND:
            append_to_csv({
                "mode": "Fallback", "status": "candidate",
                "score": d["score"], "area": d["area"],
                "angle_diff_deg": d["angle_change"],
                "film_strain_a_pct": ea, "film_strain_b_pct": eb,
                "film_dgamma_deg": d["angle_change"],
                "Mf": d["Mf"].tolist(), "Ms": d["Ms"].tolist(),
            })

    # Auto-pick
    chosen = None
    if AUTOPICK_MAX_ABS_STRAIN is not None and AUTOPICK_MAX_AREA is not None:
        for d in near:
            ea, eb = d["film_strain"]
            if max(abs(ea), abs(eb)) <= AUTOPICK_MAX_ABS_STRAIN and d["area"] <= AUTOPICK_MAX_AREA:
                chosen = d
                break

    if chosen is None:
        append_to_csv({"mode": "Fallback", "status": "no_pick_under_thresholds"})
        print("No candidate met the auto-pick thresholds. Adjust thresholds or pick manually.")
        append_to_csv({"mode": "RUN", "status": "end"})
        sys.exit(0)

    print("[OK] Auto-picked near-commensurate candidate under thresholds:")
    ea, eb = chosen["film_strain"]
    print(f"     score={chosen['score']:.3f}, area~{chosen['area']:.1f} Å², "
          f"film εa={ea:+.2f}%, εb={eb:+.2f}%, Δγ={chosen['angle_change']:+.2f}°")

    film_sc = build_supercell(film, chosen["Mf"])
    sub_sc  = build_supercell(sub,  chosen["Ms"])
    sub_cell = sub_sc.lattice.matrix.copy()
    film_sc  = replace_ab(film_sc, sub_cell[0], sub_cell[1])
    film_sc.to(fmt="poscar", filename="POSCAR_graphene_SC")
    sub_sc.to(fmt="poscar",  filename="POSCAR_WO3_SC")
    print("[OK] Wrote POSCAR_graphene_SC and POSCAR_WO3_SC (near-commensurate auto-pick)")

    append_to_csv({
        "mode": "Fallback", "status": "picked",
        "score": chosen["score"], "area": chosen["area"],
        "angle_diff_deg": chosen["angle_change"],
        "film_strain_a_pct": chosen["film_strain"][0],
        "film_strain_b_pct": chosen["film_strain"][1],
        "film_dgamma_deg": chosen["angle_change"],
        "Mf": chosen["Mf"].tolist(), "Ms": chosen["Ms"].tolist(),
    })

    if MAKE_BILAYER:
        write_quick_bilayer()

    append_to_csv({"mode": "RUN", "status": "end"})

if __name__ == "__main__":
    main()
