from pymatgen.io.vasp import Chgcar, Poscar
import numpy as np

# Load the CHGCAR file
chgcar = Chgcar.from_file("CHGCAR")

# Get the charge density data (returns a 3D numpy array)
charge_density = chgcar.data['total']

# Calculate the total charge
total_charge = np.sum(charge_density)
print(f"Total charge: {total_charge}")

# Load the POSCAR file
poscar = Poscar.from_file("POSCAR")

# Find the uppermost Au atom and the lowest S atom
au_z_coordinates = [atom[2] for atom in poscar.structure.frac_coords if atom.species_string == "Au"]
s_z_coordinates = [atom[2] for atom in poscar.structure.frac_coords if atom.species_string == "S"]
uppermost_au_z = max(au_z_coordinates)
lowest_s_z = min(s_z_coordinates)

# Take the average of the uppermost Au atom and the lowest S atom to find the interface
interface_z = (uppermost_au_z + lowest_s_z) / 2

# Convert the interface z-coordinate from fractional to Cartesian coordinates
interface_z_cart = interface_z * poscar.structure.lattice.matrix[2][2]

# Convert the interface z-coordinate from Cartesian to grid coordinates
interface_index = int(interface_z_cart / chgcar.structure.lattice.matrix[2][2] * charge_density.shape[-1])

# Calculate the charge above and below the interface
charge_below_interface = np.sum(charge_density[:, :, :interface_index])
charge_above_interface = np.sum(charge_density[:, :, interface_index:])
print(f"Charge below interface: {charge_below_interface}")
print(f"Charge above interface: {charge_above_interface}")
