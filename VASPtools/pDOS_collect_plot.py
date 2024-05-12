import os
import matplotlib.pyplot as plt
from pymatgen.io.vasp.outputs import Vasprun
from cycler import cycler

def plot_pdos(vasprun_path, plot_path):
    vasprun = Vasprun(vasprun_path)
    pdos = vasprun.complete_dos.get_element_dos()

    # Setting plot parameters for aesthetics
    plt.figure(figsize=(12, 8))
    plt.rc('font', size=20)  # controls default text size
    plt.rc('axes', titlesize=20)  # fontsize of the title
    plt.rc('axes', labelsize=20)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=18)  # fontsize of the x tick labels
    plt.rc('ytick', labelsize=18)  # fontsize of the y tick labels
    plt.rc('legend', fontsize=18)  # fontsize of the legend

    # Setting a color cycle for multiple lines
    plt.gca().set_prop_cycle(cycler('color', ['b', 'g', 'r', 'c', 'm', 'y', 'k']))

#    for element, dos in pdos.items():
#        plt.plot(dos.energies - vasprun.efermi, dos.get_densities(), label=str(element), linewidth=2)


    x_min, x_max = -5, 5  # Energy range


    # First, find the maximum density within the specified range
    max_density = 0
    for element, dos in pdos.items():
        energies = dos.energies - vasprun.efermi
        densities = dos.get_densities()
        within_range = (energies >= x_min) & (energies <= x_max)
        if within_range.any():  # Check if there are points within the range
            max_density =  max(densities[within_range])
	
    print(max_density)

    # Plot normalized densities after finding max density
    for element, dos in pdos.items():
	

        energies = dos.energies - vasprun.efermi
        densities = dos.get_densities() / max_density  # Normalize
        within_range = (energies >= x_min) & (energies <= x_max)
        plt.plot(energies[within_range], densities[within_range], label=str(element), linewidth=2)
    print(element)
    plt.xlim(x_min, x_max)
    plt.ylim(0, 1)
    plt.xlabel("Energy (eV)")
    plt.ylabel("Normalized Density of States")
    plt.title(f"Normalized Projected Density of States {subdir_name}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(plot_path, dpi=300)
    plt.close()

# Iterate over subdirectories and plot normalized pDOS
for subdir in next(os.walk('.'))[1]:
    vasprun_path = os.path.join(subdir, 'vasprun.xml')
    subdir_name = os.path.basename(subdir)
    plot_filename = f'normalized_pdos_{subdir_name}.png'
    plot_path = os.path.join('.', plot_filename)

    if os.path.exists(vasprun_path):
        plot_pdos(vasprun_path, plot_path)
    else:
        print(f"vasprun.xml not found in {subdir}")
