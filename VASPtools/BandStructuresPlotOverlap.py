from pymatgen.io.vasp import BSVasprun
import matplotlib.pyplot as plt
import numpy as np

# Load the VASP band structure data
vasprun1 = BSVasprun("./PBE-prim/vasprun.xml", parse_projected_eigen=True)
vasprun2 = BSVasprun("./HSE-prim/vasprun.xml", parse_projected_eigen=True)
vasprun3 = BSVasprun("./CAM-prim/vasprun.xml", parse_projected_eigen=True)

# Extract band structures
bs1 = vasprun1.get_band_structure(kpoints_filename="./KPOINTS", line_mode=True)
bs2 = vasprun2.get_band_structure(kpoints_filename="./KPOINTS", line_mode=True)
bs3 = vasprun3.get_band_structure(kpoints_filename="./KPOINTS", line_mode=True)

# Extract spin channel (handles spin and non-spin cases)
spin = list(bs1.bands.keys())[0]  

# Get the Valence Band Maximum (VBM) for each calculation
VBM1 = bs1.get_vbm()['energy']
VBM2 = bs2.get_vbm()['energy']
VBM3 = bs3.get_vbm()['energy']

# Create figure and axis
fig, ax = plt.subplots(figsize=(8, 6))

# Function to extract band data, shift VBM to 0 eV, and plot
def plot_band_structure(bs, VBM, ax, color, label):
    kpoints = bs.distance  # Get k-point distances
    for band in bs.bands[spin]:  
        shifted_band = np.array(band) - VBM  # Shift energy by VBM
        ax.plot(kpoints, shifted_band, color=color, linewidth=1.5, alpha=0.8, label=label if band is bs.bands[spin][0] else "")

# Plot each band structure with different colors and correct VBM shift
plot_band_structure(bs1, VBM1, ax, "black", "PBE")
plot_band_structure(bs2, VBM2, ax, "blue", "HSE06")
plot_band_structure(bs3, VBM3, ax, "red", "CAM-B3LYP")

# ---- Correctly Add k-Point Labels ----
kpt_positions = []
kpt_labels = []
for i, kpt in enumerate(bs1.kpoints):
    if kpt.label is not None:  # Only use labeled k-points (Γ, X, L, etc.)
        kpt_labels.append(kpt.label)
        kpt_positions.append(bs1.distance[i])  # Corresponding k-point distance

ax.set_xticks(kpt_positions)
ax.set_xticklabels(kpt_labels, fontsize=14)

# Beautify the plot
ax.set_title("Overlapping Band Structures (VBM at 0 eV)", fontsize=16)
ax.set_xlabel("Wave Vector", fontsize=18)
ax.set_ylabel("Energy (eV)", fontsize=18)
ax.axhline(0, color="black", linestyle="--", linewidth=1, label="VBM (0 eV)")
ax.legend()
ax.grid(alpha=0.3)
ax.tick_params(axis='both', which='major', labelsize=18)
ax.tick_params(axis='y', which='minor', labelsize=12)

# Set x-axis and y-axis limits
ax.set_ylim(-18, 18)  # Adjust energy range
ax.set_xlim(0, max(bs1.distance))  # Fit the plot within the box

ax.plot([], [], color="black", linewidth=2,linestyle=":",  label="PBE")  # Empty plot for legend
ax.plot([], [], color="blue", linewidth=2, label="HSE")
ax.plot([], [], color="red", linewidth=2, label="CAM-B3LYP")
ax.legend(loc="lower left", fontsize=12)  # Adjust location and font size



# Save and Show
plt.savefig('band_structure_shifted.png', dpi=600, bbox_inches='tight')
plt.show()
