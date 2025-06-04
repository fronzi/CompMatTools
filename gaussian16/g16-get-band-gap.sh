#!/bin/bash

# Check if a file name is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <Gaussian log file>"
    exit 1
fi

LOGFILE=$1

# Function to extract HOMO and LUMO energies


extract_homo_lumo() {
    local orbital_type=$1
    local __homo_var=$2
    local __lumo_var=$3

    local last_scf_start
    last_scf_start=$(grep -n "SCF Done" "$LOGFILE" | tail -1 | cut -d':' -f1)

    local homo_line
    homo_line=$(awk "NR>${last_scf_start}" "$LOGFILE" \
                | grep "${orbital_type}  occ. eigenvalues" \
                | tail -1)
    local homo=""
    if [[ -n "$homo_line" ]]; then
        homo=$(echo "$homo_line" | awk '{print $NF}')
    fi

    local lumo_line
    lumo_line=$(awk "NR>${last_scf_start}" "$LOGFILE" \
                | grep "${orbital_type} virt. eigenvalues" \
                | head -1)
    # suppress “field -4” errors
    local lumo=$(echo "$lumo_line" | awk '{print $(NF-4)}' 2>/dev/null)

    if [[ -z "$homo" || -z "$lumo" ]]; then
        eval "$__homo_var=''"
        eval "$__lumo_var=''"
    else
        echo "${orbital_type} orbitals (in Hartree):"
        echo "  HOMO: $homo"
        echo "  LUMO: $lumo"
        eval "$__homo_var=$homo"
        eval "$__lumo_var=$lumo"
    fi
}

# Extract HOMO and LUMO energies for alpha and beta orbitals
extract_homo_lumo "Alpha" HOMO_ALPHA LUMO_ALPHA
extract_homo_lumo "Beta"  HOMO_BETA  LUMO_BETA

# Check alpha extraction succeeded
if [[ -z "$HOMO_ALPHA" || -z "$LUMO_ALPHA" ]]; then
    echo "Error: Alpha orbitals not found in log file. Cannot proceed."
    exit 1
fi

# Convert HOMO/LUMO from Hartree to eV for alpha
HOMO_ALPHA_eV=$(echo "$HOMO_ALPHA * 27.2114" | bc -l)
LUMO_ALPHA_eV=$(echo "$LUMO_ALPHA * 27.2114" | bc -l)

# Calculate alpha band gap (in Hartree and eV)
ALPHA_GAP_HA=$(echo "$LUMO_ALPHA - $HOMO_ALPHA" | bc -l)
ALPHA_GAP_eV=$(echo "$ALPHA_GAP_HA * 27.2114" | bc -l)

echo ""
echo "Alpha Band Gap:"
echo "  HOMO (Ha): $HOMO_ALPHA"
echo "  LUMO (Ha): $LUMO_ALPHA"
echo "  Gap  (Ha): $ALPHA_GAP_HA"
echo "  HOMO (eV): $HOMO_ALPHA_eV"
echo "  LUMO (eV): $LUMO_ALPHA_eV"
echo "  Gap  (eV): $ALPHA_GAP_eV"

# If beta orbitals not found, skip beta and spin-flip calculations
if [[ -z "$HOMO_BETA" || -z "$LUMO_BETA" ]]; then
    echo ""
    echo "No Beta orbitals found. Skipping Beta and spin-flip calculations."
    exit 0
fi

# Convert HOMO/LUMO from Hartree to eV for beta
HOMO_BETA_eV=$(echo "$HOMO_BETA * 27.2114" | bc -l)
LUMO_BETA_eV=$(echo "$LUMO_BETA * 27.2114" | bc -l)

# Calculate beta band gap (in Hartree and eV)
BETA_GAP_HA=$(echo "$LUMO_BETA - $HOMO_BETA" | bc -l)
BETA_GAP_eV=$(echo "$BETA_GAP_HA * 27.2114" | bc -l)

echo ""
echo "Beta Band Gap:"
echo "  HOMO (Ha): $HOMO_BETA"
echo "  LUMO (Ha): $LUMO_BETA"
echo "  Gap  (Ha): $BETA_GAP_HA"
echo "  HOMO (eV): $HOMO_BETA_eV"
echo "  LUMO (eV): $LUMO_BETA_eV"
echo "  Gap  (eV): $BETA_GAP_eV"

# Calculate and display the overall HOMO-LUMO gap (smallest LUMO - largest HOMO)
highest_homo=$(echo "$HOMO_ALPHA $HOMO_BETA" | awk '{if ($1>$2) print $1; else print $2}')
lowest_lumo=$(echo "$LUMO_ALPHA $LUMO_BETA" | awk '{if ($1<$2) print $1; else print $2}')
OVERALL_GAP_eV=$(echo "($lowest_lumo - $highest_homo) * 27.2114" | bc -l)

echo ""
echo "Overall HOMO-LUMO Gap:"
echo "  Largest HOMO (Ha): $highest_homo"
echo "  Smallest LUMO (Ha): $lowest_lumo"
echo "  Gap      (eV): $OVERALL_GAP_eV"

# Calculate and display the spin-flip band gaps
#  - Beta LUMO minus Alpha HOMO
SPIN_FLIP_BETA_ALPHA_HA=$(echo "$LUMO_BETA - $HOMO_ALPHA" | bc -l)
SPIN_FLIP_BETA_ALPHA_eV=$(echo "$SPIN_FLIP_BETA_ALPHA_HA * 27.2114" | bc -l)

#  - Alpha LUMO minus Beta HOMO
SPIN_FLIP_ALPHA_BETA_HA=$(echo "$LUMO_ALPHA - $HOMO_BETA" | bc -l)
SPIN_FLIP_ALPHA_BETA_eV=$(echo "$SPIN_FLIP_ALPHA_BETA_HA * 27.2114" | bc -l)

echo ""
echo "Spin-Flip Band Gaps:"
echo "  Beta LUMO - Alpha HOMO:"
echo "    Gap (Ha): $SPIN_FLIP_BETA_ALPHA_HA"
echo "    Gap (eV): $SPIN_FLIP_BETA_ALPHA_eV"
echo ""
echo "  Alpha LUMO - Beta HOMO:"
echo "    Gap (Ha): $SPIN_FLIP_ALPHA_BETA_HA"
echo "    Gap (eV): $SPIN_FLIP_ALPHA_BETA_eV"
