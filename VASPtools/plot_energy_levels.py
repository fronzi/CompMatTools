import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os


def parse_eigenval_to_dataframe(file_path):
    """
    Parses the EIGENVAL file formatted with columns: band number, spin-up energy, spin-down energy, 
    spin-up occupancy, spin-down occupancy, and stores the data in a Pandas DataFrame.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Extract metadata from line 6
    metadata_line = lines[5]
    num_electrons, num_kpoints, num_bands = map(int, metadata_line.split())

    # Initialize an empty list to hold data for the DataFrame
    data = []
    
    # Skip the header and parse the k-point and band data
    current_line = 7  # Start after the header
    for kpoint_index in range(1, num_kpoints + 1):
        current_line += 1  # Skip the k-point header line
        for band_index in range(1, num_bands + 1):
            split_line = lines[current_line].split()
            band_number = int(split_line[0])
            spin_up_energy = float(split_line[1])
            spin_down_energy = float(split_line[2])
            spin_up_occupancy = float(split_line[3])
            spin_down_occupancy = float(split_line[4])
            data.append([
                kpoint_index,
                band_number,
                spin_up_energy,
                spin_down_energy,
                spin_up_occupancy,
                spin_down_occupancy
            ])
            current_line += 1

    # Create a Pandas DataFrame
    df = pd.DataFrame(
        data,
        columns=[
            'k-point',
            'band number',
            'spin-up energy',
            'spin-down energy',
            'spin-up occupancy',
            'spin-down occupancy'
        ]
    )
    return df, num_electrons, num_bands

def get_fermi_energy_from_outcar(file_path='OUTCAR'):
    """
    Extracts the Fermi energy from the OUTCAR file.
    """
    with open(file_path, 'r') as file:
        for line in file:
            if 'E-fermi' in line:
                return float(line.split()[2])  # Extract the third value
    raise ValueError("Fermi energy not found in OUTCAR")

def find_homo_lumo(df, num_electrons):
    """
    Finds the HOMO and LUMO indices and their energies from the spin-up channel.
    """
    homo_index = num_electrons // 2 - 1  # HOMO is the last occupied state
    lumo_index = homo_index + 1          # LUMO is the first unoccupied state
    
    spin_up_homo_energy = df[df['band number'] == homo_index + 1]['spin-up energy'].iloc[0]
    spin_up_lumo_energy = df[df['band number'] == lumo_index + 1]['spin-up energy'].iloc[0]
    
    print(f"HOMO index: {homo_index}, LUMO index: {lumo_index}")
    print(f"HOMO energy (spin-up): {spin_up_homo_energy:.4f}, LUMO energy (spin-up): {spin_up_lumo_energy:.4f}")
    
    return homo_index, lumo_index, spin_up_homo_energy, spin_up_lumo_energy



def plot_energy_levels(df, fermi_energy, name):
    """
    Plots the energy levels for the spin-up and spin-down channels, referenced to the Fermi energy.
    Annotates the differences between the Fermi energy and relative energy levels.
    """
    # Adjust energies relative to Fermi energy
    df['spin-up energy'] -= fermi_energy
    df['spin-down energy'] -= fermi_energy

    # Group by k-points
    grouped = df.groupby('k-point')

    plt.figure(figsize=(8, 6))
    for kpoint, group in grouped:
        # Scatter plot for spin-up
        plt.scatter(
            group['band number'],
            group['spin-up energy'],
            label=f'Spin-up',
            alpha=0.7,
            marker='^',
            color='blue',  # Set the color of the markers
            s=120          # Set the size of the markers
        )
        # Annotate differences for spin-up
        for _, row in group.iterrows():
            diff = row['spin-up energy']  # Difference relative to Fermi energy
            plt.text(
                row['band number'], row['spin-up energy'] - 0.6, f'{diff:.2f}',
                fontsize=10, color='blue', ha='left', va='bottom'
            )

        # Scatter plot for spin-down
        plt.scatter(
            group['band number'],
            group['spin-down energy'],
            label=f'Spin-down',
            alpha=0.7,
            marker='v',
            color='red',    # Set the color of the markers
            s=120           # Set the size of the markers
        )
        # Annotate differences for spin-down
        for _, row in group.iterrows():
            diff = row['spin-down energy']  # Difference relative to Fermi energy
            plt.text(
                row['band number'], row['spin-down energy'] + 0.6, f'{diff:.2f}',
                fontsize=10, color='red', ha='left', va='top'
            )

    plt.axhline(0, color='gray', linestyle='--', label='Fermi Level')
    plt.xlabel('Band Number', fontsize=16)
    plt.ylabel('Energy (eV)', fontsize=16)
    plt.xlim(1021, 1025)
    plt.xticks(ticks=range(1021, 1026), fontsize=12)
    plt.yticks(fontsize=14)
#    plt.ylim(-6, 6)
    plt.legend(loc='upper left', fontsize=14)
    plt.grid()
    plt.savefig(f"{name}_vbm")
    plt.show()


# Main execution

# Get the current working directory (full path)
current_path = os.getcwd()

# Extract only the name of the current directory
current_dir_name = os.path.basename(current_path)

name=current_dir_name

eigenval_path = 'EIGENVAL'  # Replace with the path to your EIGENVAL file
outcar_path = 'OUTCAR'      # Replace with the path to your OUTCAR file

# Parse the EIGENVAL file into a DataFrame
df, num_electrons, num_bands = parse_eigenval_to_dataframe(eigenval_path)

# Debug: Display the first few rows of the DataFrame
print(df.head())

# Extract Fermi energy
fermi_energy = get_fermi_energy_from_outcar(outcar_path)
fermi_energy = -1.46
print(f"Fermi Energy: {fermi_energy:.4f} eV")

# Find HOMO and LUMO
homo_index, lumo_index, homo_energy, lumo_energy = find_homo_lumo(df, num_electrons)

# Plot energy levels
plot_energy_levels(df, fermi_energy, name)
