#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 11:55:44 2023

@author: marco
"""



import os
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.core import OrbitalType
from pymatgen.electronic_structure.plotter import DosPlotter
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import itertools
import csv
import warnings
from pymatgen.io.vasp import Vasprun
import traceback

from scipy.integrate import simps

def overlap_integral(x, f, g):
    """
    Compute the overlap integral of two functions, f(x) and g(x),
    over the domain specified by x.
    """
    return simps(f * g, x)

def calculate_overlap(dos1, dos2, e_min=-6, e_max=6):
    """
    Calculate the overlap of two DOS curves.
    
    Args:
        dos1, dos2: DOS objects for the elements/orbitals of interest.
        e_min, e_max: Energy range of interest.
    
    Returns:
        Overlap value.
    """
    energies = dos1.energies - dos1.efermi  # Shift by Fermi energy
    
    # Filter based on the energy range
    indices = (energies >= e_min) & (energies <= e_max)
    
    key1 = list(dos1.densities.keys())[0]
    key2 = list(dos2.densities.keys())[0]

    # Get the densities for each DOS and filter them based on the energy range
    densities1 = dos1.densities[key1][indices]
    densities2 = dos2.densities[key2][indices]
    
    # Return the overlap integral for the filtered densities and energies
    return overlap_integral(energies[indices], densities1, densities2)





    
    # Return the overlap integral for the filtered densities and energies
    return overlap_integral(energies[indices], densities1, densities2)

# Example of usage:
# Assume we have dos_object_element1 and dos_object_element2 as DOS objects corresponding to two different elements.
# overlap = calculate_overlap(dos_object_element1, dos_object_element2)
# print(f"Overlap: {overlap:.4f}")

def calculate_all_overlaps(cdos, e_min=-6, e_max=6):
    """
    Calculate overlaps for all possible element-orbital pairs.
    
    Args:
        cdos: CompleteDOS object.
        e_min, e_max: Energy range of interest.
    
    Returns:
        Dictionary with overlaps.
    """
    overlaps = {}
    elements = cdos.structure.composition.elements
    orbitals = [OrbitalType.s, OrbitalType.p, OrbitalType.d]  # you can add more if needed
    
    # Generate all possible element-orbital pairs
    pairs = list(itertools.product(elements, orbitals))
    
    for i, pair1 in enumerate(pairs):
        for j, pair2 in enumerate(pairs):
            if j > i:  # To avoid calculating overlap for the same pair twice
                key = f"{pair1[0]}-{pair1[1]} vs {pair2[0]}-{pair2[1]}"
                dos1 = cdos.get_element_spd_dos(str(pair1[0]))[pair1[1]]
                dos2 = cdos.get_element_spd_dos(str(pair2[0]))[pair2[1]]
                overlap = calculate_overlap(dos1, dos2, e_min, e_max)
                overlaps[key] = overlap
                print(f"{key}: Overlap = {overlap:.4f}")
                
    return overlaps




def plot_version1(v, path):
    try:
        # Read the VASP output
        v = Vasprun(path + '/vasprun.xml')
      

        # Get the complete density of states
        cdos = v.complete_dos

        # Get the composition of the structure
        composition = cdos.structure.composition

        # Extract the list of elements in the structure
        elements = composition.elements

        # Initialize the DosPlotter
        plotter = DosPlotter()

        # Initialize variable to keep track of the maximum pDOS value
        max_pDOS = 0

        # Define energy limits
        e_min, e_max = -6, 6

        # Define a list of elements you want to plot; leave it empty to plot all elements
        selected_elements = []  # Add the elements you want to plot here

        # Generate the plot with modified figure size
        fig, ax = plt.subplots(figsize=(8, 6))

        # Loop through each element and add its normalized s, p, d DOS to the plotter
        for element in elements:
            element_symbol = str(element)
            
            # Check if the element is in the selected list or if all elements should be plotted
            if not selected_elements or element_symbol in selected_elements:
                element_spd = cdos.get_element_spd_dos(element_symbol)
                num_atoms = composition[element_symbol]
                
                for orbital_type in [OrbitalType.s, OrbitalType.p, OrbitalType.d]:
                    if orbital_type in element_spd:
                        densities = element_spd[orbital_type].densities
                        energies = element_spd[orbital_type].energies - cdos.efermi
                        
                        for spin in densities.keys():
                            densities[spin] /= num_atoms
                            
                            filtered_densities = [density for energy, density in zip(energies, densities[spin]) if e_min <= energy <= e_max]
                            
                            if filtered_densities:
                                max_pDOS = max(max_pDOS, max(filtered_densities))
                        
                        plotter.add_dos(f"{element_symbol} {orbital_type.name}", element_spd[orbital_type])

        # # Generate the plot
        # plt = plotter.get_plot(xlim=[e_min, e_max], ylim=[0, max_pDOS])
        # Generate the plot
        plot_obj = plotter.get_plot(xlim=[e_min, e_max], ylim=[0, max_pDOS])
        # Get the current axes again after generating the plot
        ax = plot_obj.gca()

        # Get the current axes again after generating the plot
        ax = plt.gca()

        # Modify tick size and labels
        ax.tick_params(axis='both', which='major', labelsize=30)
        ax.tick_params(axis='both', which='minor', labelsize=16)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.set_xlabel('Energy (eV)', fontsize=30)
        ax.set_ylabel('pDOS (arb units)', fontsize=30)
        ax.legend(loc='upper right', fontsize=26)

        # Change the thickness of the plot box
        for spine in ['top', 'bottom', 'left', 'right']:
            ax.spines[spine].set_linewidth(1)

        ax.grid(True)
        plot_obj.tight_layout()
        # directory = os.path.dirname(vasprun_path)
 # Instead of directory, use the path
        save_path = os.path.join(path, 'Element_dos_normalized.png')
        plot_obj.savefig(save_path, dpi=600)
        
        plot_obj.savefig(path + '/Element_spd_dos_normalized.png', dpi=600)
        plt.close()  # Corrected from plot_obj.close()

        
    except Exception as e:
        print(f"An error occurred in directory {path}: {e}")

def plot_version2(vasprun_path):
    v = Vasprun(vasprun_path)
    
    
  # Get the complete density of states
    cdos = v.complete_dos

  

    # Get the composition of the structure
    composition = cdos.structure.composition

    # Extract the list of elements in the structure
    elements = composition.elements

    # Initialize the DosPlotter
    plotter = DosPlotter()

    # Initialize variable to keep track of the maximum pDOS value
    max_pDOS = 0

    # Define energy limits
    e_min, e_max = -6, 6

    # Define a list of elements you want to plot; leave it empty to plot all elements
    selected_elements = []  # Add the elements you want to plot here

    # Loop through each element and add its normalized DOS to the plotter
    for element in elements:
        element_symbol = str(element)
        
        # Check if the element is in the selected list or if all elements should be plotted
        if not selected_elements or element_symbol in selected_elements:
            element_dos = cdos.get_element_dos()[element]
            num_atoms = composition[element_symbol]
            
            # Normalize the DOS by the number of atoms of the element
            densities = element_dos.densities
            energies = element_dos.energies - cdos.efermi  # Shift by Fermi energy
            
            for spin in densities.keys():
                densities[spin] /= num_atoms
                
                # Filter densities based on energy limits
                filtered_densities = [density for energy, density in zip(energies, densities[spin]) if e_min <= energy <= e_max]
                
                if filtered_densities:
                    max_pDOS = max(max_pDOS, max(filtered_densities))  # Update max_pDOS
            
            plotter.add_dos(f"{element_symbol}", element_dos)

    # Generate the plot with modified figure size
    fig, ax = plt.subplots(figsize=(8, 6))

    directory = os.path.dirname(vasprun_path)  # Extract the directory from the vasprun_path
    # Generate the plot
    plot_obj = plotter.get_plot(xlim=[e_min, e_max], ylim=[0, max_pDOS])
    
    # Get the current axes again after generating the plot
    ax = plot_obj.gca()

     # Apply your matplotlib customizations
    ax.tick_params(axis='both', which='major', labelsize=30)
    ax.tick_params(axis='both', which='minor', labelsize=16)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.set_xlabel('Energy (eV)', fontsize=30)
    ax.set_ylabel('pDOS (arb units)', fontsize=30)


    # Change the thickness of the plot box
    for spine in ['top', 'bottom', 'left', 'right']:
        ax.spines[spine].set_linewidth(1)
    
    ax.grid(True)
    plot_obj.tight_layout()

  # Save the figure
    plot_obj.savefig(os.path.join(directory, 'Element_spd_dos_normalized.png'), dpi=600)
    plt.close()


def generate_plot(path):
    try:
        # Print the exact path we are trying to read
        vasprun_path = os.path.join(path, 'vasprun.xml')
        print(f"Trying to read {vasprun_path}")
        
        # Read the VASP output
        v = Vasprun(vasprun_path)

        # Extract the complete density of states
        cdos = v.complete_dos

        # Generate different versions of the plot
        plot_version1(v, path)
        plot_version2(vasprun_path)
        overlaps = calculate_all_overlaps(cdos)

        if overlaps:  # Check if there's any overlap data
            # Saving overlaps to a CSV file
            with open(os.path.join(path, "overlaps.csv"), "w", newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(["Pair", "Overlap"])  # header row
                for key, value in overlaps.items():
                    writer.writerow([key, value])
            
            print(f"Overlaps saved to {os.path.join(path, 'overlaps.csv')}")

        else:
            print(f"No overlaps data for directory {path}")

    except Exception as e:
        print(f"An error occurred in directory {path}: {e}")
        traceback.print_exc()  # This will print the full traceback

        
        
        
def crawl_directory(root_path):
    for root, dirs, files in os.walk(root_path):
        if 'dos' in dirs:
            dos_path = os.path.join(root, 'dos')
            print(f"Generating plot for {dos_path}")
            generate_plot(dos_path)



with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    v = Vasprun('/Users/marco/Dropbox/Work/WORKING_DIR_UTS/InSe/test/InSe-1/dos/vasprun.xml')

# Define the root path to start the search
root_path = '/Users/marco/Dropbox/Work/WORKING_DIR_UTS/InSe/test/'

# Start crawling from the root directory
crawl_directory(root_path)
