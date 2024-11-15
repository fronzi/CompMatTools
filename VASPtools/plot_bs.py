from pymatgen.io.vasp import BSVasprun
from pymatgen.electronic_structure.plotter import BSPlotter
import matplotlib.pyplot as plt

# Load the VASP band structure data
vasprun = BSVasprun("./vasprun.xml", parse_projected_eigen=True)
bs = vasprun.get_band_structure(kpoints_filename="./KPOINTS", line_mode=True)

# Create the plot
plotter = BSPlotter(bs)
ax = plotter.get_plot(vbm_cbm_marker=False, ylim=(-4, 8))  # No line_width here

# Increase line thickness by updating each line in the Axes object
for line in ax.get_lines():
    line.set_linewidth(2.5)  # Set desired line thickness

# Beautify the plot
ax.set_title("Electronic Band Structure", fontsize=16)
ax.set_xlabel("Wave Vector", fontsize=14)
ax.set_ylabel("Energy (eV)", fontsize=14)
ax.grid(alpha=0.3)  # Add light grid lines
ax.tick_params(axis='both', which='major', labelsize=12)

# Save the figure
plt.savefig('band_structure.png', dpi=600, bbox_inches='tight')

# Show the plot
plt.show()
