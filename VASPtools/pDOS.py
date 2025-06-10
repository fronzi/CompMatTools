import os
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.plotter import DosPlotter

# === Configurable Parameters ===
energy_range = (-10, 10)   # eV range around Fermi level
dos_sigma = 0.05           # Gaussian smearing
stack_dos = False          # Line plots (set to True for stacked area)
plot_title = "Element-Resolved Density of States"

# === Load VASP data ===
vasprun = Vasprun("vasprun.xml")
complete_dos = vasprun.complete_dos

# === Initialize plotter ===
dos_plotter = DosPlotter(sigma=dos_sigma, stack=stack_dos)

# Add total DOS
dos_plotter.add_dos("Total DOS", complete_dos)

# Add element-projected DOS
element_dos_dict = complete_dos.get_element_dos()
for el in element_dos_dict:
    dos_plotter.add_dos(el, element_dos_dict[el])

# === Generate and style the plot ===
plot_data = dos_plotter.get_plot(xlim=energy_range, ylim=(0, None))

# Plot styling
plt.title(plot_title, fontsize=16, fontweight='bold')
plt.xlabel("Energy (eV)", fontsize=14)
plt.ylabel("DOS (states/eV)", fontsize=14)
plt.axvline(x=0, color='k', linestyle='--', linewidth=1)  # Fermi level
plt.grid(True, linestyle=':', alpha=0.6)

# === Remove duplicate legend entries ===
handles, labels = plt.gca().get_legend_handles_labels()
unique = dict(zip(labels, handles))
plt.legend(unique.values(), unique.keys(), fontsize=10)

plt.tight_layout()

# === Save figure with current directory name ===
current_dir_name = os.path.basename(os.getcwd())
output_file = f"{current_dir_name}_pdos.png"
plt.savefig(output_file, dpi=300)

# === Show the plot ===
plt.show()
