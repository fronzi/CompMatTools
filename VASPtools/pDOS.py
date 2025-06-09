import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.plotter import DosPlotter

# === Load data from vasprun.xml ===
vasprun = Vasprun("vasprun.xml", parse_projected=True)
complete_dos = vasprun.complete_dos

# === Total and Elemental DOS ===
dos_plotter = DosPlotter(sigma=0.05, stack=False)

# Add total DOS
dos_plotter.add_dos("Total DOS", complete_dos)

# Add element-specific projected DOS (change elements as needed)
elements_of_interest = ["Cu", "C"]
for el in elements_of_interest:
    dos_plotter.add_dos(el, complete_dos.get_element_dos(el))

# === Plot with Custom Styling ===
plot_data = dos_plotter.get_plot(xlim=(-10, 10), ylim=(0, None), plt=plt)

# Custom plot style
plt.title("Projected Density of States", fontsize=16, fontweight='bold')
plt.xlabel("Energy (eV)", fontsize=14)
plt.ylabel("DOS (states/eV)", fontsize=14)
plt.axvline(x=0, color='k', linestyle='--', linewidth=1)  # Fermi level
plt.legend(fontsize=10)
plt.grid(True, linestyle=':', alpha=0.7)
plt.tight_layout()

# === Save or Show ===
plt.savefig("pdos_beautiful.png", dpi=300)
plt.show()
