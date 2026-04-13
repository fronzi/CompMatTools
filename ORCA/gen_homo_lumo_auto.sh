#!/bin/bash
# orca_homo_lumo_cube.sh
#
# Automatically extracts HOMO/LUMO indices from an ORCA output file,
# then generates cube files for both orbitals via orca_plot.
#
# Usage: ./orca_homo_lumo_cube.sh [basename]
#   basename : filename stem without extension (default: carbon_fluoride)
#
# Requires: orca_plot on PATH; ${BASENAME}.gbw and ${BASENAME}.out present.
#
# orca_plot interactive menu (version-stable sequence used here):
#   1  → Plot type : molecular orbital
#   1  → MO type   : alpha (use 2 for beta in open-shell)
#   2  → Enter MO index
#   N  → MO number (0-based, matching ORCA orbital energies block)
#   4  → Grid quality
#   120→ Grid points per dimension
#   11 → (confirm / write cube)
#   12 → Quit
# ---------------------------------------------------------------------------

set -euo pipefail

# --- CONFIGURATION ---
BASENAME="${1:-carbon_fluoride}"
GBWFILE="${BASENAME}.gbw"
OUTFILE="${BASENAME}.out"

# Fallback indices if auto-detection fails
HOMO_FALLBACK=20
LUMO_FALLBACK=21

# orca_plot grid resolution (increase for publication quality)
GRID_PTS=120

# ---------------------------------------------------------------------------
# Validate input files
# ---------------------------------------------------------------------------
if [[ ! -f "$GBWFILE" ]]; then
    echo "[ERROR] GBW file not found: ${GBWFILE}" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Auto-detect HOMO/LUMO from ORCA output
# ---------------------------------------------------------------------------
detect_homo_lumo() {
    local outfile="$1"

    # The ORCA "ORBITAL ENERGIES" block has the structure:
    #
    #   ORBITAL ENERGIES
    #   ----------------
    #
    #     NO   OCC          E(Eh)            E(eV)
    #      0   2.0000     -26.405...
    #      ...
    #     20   2.0000      -0.456...   ← HOMO (last line with OCC > 0)
    #     21   0.0000       0.123...   ← LUMO
    #
    # Index column ($1) is 0-based and matches orca_plot's MO numbering.

    awk '
        /^ORBITAL ENERGIES[[:space:]]*$/ { in_block=1; next }
        in_block && /NO[[:space:]]+OCC/  { header_seen=1; next }
        in_block && header_seen && /^[[:space:]]*[0-9]+[[:space:]]+[0-9]/ {
            occ = $2 + 0
            if (occ > 1e-6) homo = $1 + 0
        }
        in_block && header_seen && /^[[:space:]]*$/ { exit }
        END {
            if (homo == "") { exit 1 }
            print homo
        }
    ' "$outfile"
}

if [[ -f "$OUTFILE" ]]; then
    echo "Parsing orbital occupations from: ${OUTFILE}"

    if HOMO=$(detect_homo_lumo "$OUTFILE" 2>/dev/null); then
        LUMO=$(( HOMO + 1 ))
        echo "  → HOMO : orbital ${HOMO}"
        echo "  → LUMO : orbital ${LUMO}"
    else
        echo "[WARNING] Could not parse orbital energies block." >&2
        echo "          Falling back to HOMO=${HOMO_FALLBACK}, LUMO=${LUMO_FALLBACK}" >&2
        HOMO=$HOMO_FALLBACK
        LUMO=$LUMO_FALLBACK
    fi
else
    echo "[WARNING] Output file not found: ${OUTFILE}" >&2
    echo "          Falling back to HOMO=${HOMO_FALLBACK}, LUMO=${LUMO_FALLBACK}" >&2
    HOMO=$HOMO_FALLBACK
    LUMO=$LUMO_FALLBACK
fi

# ---------------------------------------------------------------------------
# Sanity check: ensure LUMO index doesn't exceed the basis set dimension.
# orca_plot will error out gracefully, but we can warn early.
# ---------------------------------------------------------------------------
if [[ -f "$OUTFILE" ]]; then
    N_MOS=$(awk '/^Number of basis functions/{print $NF; exit}' "$OUTFILE" 2>/dev/null || echo "")
    if [[ -n "$N_MOS" && "$LUMO" -ge "$N_MOS" ]]; then
        echo "[WARNING] LUMO index ${LUMO} >= number of MOs (${N_MOS}). Check your system." >&2
    fi
fi

# ---------------------------------------------------------------------------
# Helper: run orca_plot interactively
# ---------------------------------------------------------------------------
run_orca_plot() {
    local label="$1"    # HOMO or LUMO
    local mo_idx="$2"   # 0-based MO index

    echo ""
    echo "Generating ${label} (orbital ${mo_idx})..."

    orca_plot "${GBWFILE}" -i <<EOF
1
1
2
${mo_idx}
4
${GRID_PTS}
11
12
EOF

    # orca_plot names the cube file as <basename>.mo<N>a.cube (alpha)
    # Rename to a more descriptive name for clarity
    # bash 3-safe lowercase: avoids ${label,,} which requires bash 4+
    local label_lower
    label_lower=$(echo "$label" | tr '[:upper:]' '[:lower:]')

    local default_cube="${BASENAME}.mo${mo_idx}a.cube"
    local target_cube="${BASENAME}_${label_lower}.cube"

    if [[ -f "$default_cube" ]]; then
        mv "$default_cube" "$target_cube"
        echo "  → Cube file: ${target_cube}"
    else
        echo "  [WARNING] Expected cube file not found: ${default_cube}" >&2
        echo "            orca_plot may have used a different naming convention." >&2
        ls -1 *.cube 2>/dev/null | tail -5 || true
    fi
}

# ---------------------------------------------------------------------------
# Generate cube files
# ---------------------------------------------------------------------------
run_orca_plot "HOMO" "$HOMO"
run_orca_plot "LUMO" "$LUMO"

echo ""
echo "Done."
echo "  HOMO cube : ${BASENAME}_homo.cube"
echo "  LUMO cube : ${BASENAME}_lumo.cube"
echo ""
echo "Render with VESTA, VMD, or Avogadro."
echo "For VMD isosurface: isovalue ±0.04 (adjust to ~80% electron density capture)."
