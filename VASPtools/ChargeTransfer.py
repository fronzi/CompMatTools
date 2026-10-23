from pymatgen.io.vasp import Chgcar, Poscar, Outcar
import numpy as np

    


M=['Au','Ni','Cu','Ag','Pd','Al']


for m in M:
          
    # Load the OUTCAR file
    outcar = Outcar("./"+m+"-MoS2/OUTCAR")
    # Load the POSCAR file
    poscar = Poscar.from_file("./"+m+"-MoS2/CONTCAR")
    volume = poscar.structure.volume
    
    # Extract the NELECT value
    nelect = outcar.nelect
    print(f"Crystal: {m}-MoS2")
    print(f"NELECT: {nelect}")

    # Load the CHGCAR file
    chgcar = Chgcar.from_file("./"+m+"-MoS2/CHGCAR")
    chgcar_metal = Chgcar.from_file("./"+m+"-MoS2/CHGCAR."+m)
    chgcar_mos2 = Chgcar.from_file("./"+m+"-MoS2/CHGCAR.MoS2")
    
    
    
    
    # Get the charge density data (returns a 3D numpy array)
    charge_density = chgcar.data['total']
    charge_density_metal = chgcar_metal.data['total']
    charge_density_mos2 = chgcar_mos2.data['total']
    
    
    # Calculate the total charge
    total_charge = np.sum(charge_density)/volume
    total_charge_metal = np.sum(charge_density_metal)/volume
    total_charge_mos2 = np.sum(charge_density_mos2)/volume
    
    print(f"Total charge: {total_charge}")
    
    
    
    # Find the uppermost Au atom and the lowest S atom
    au_z_coordinates = [site.z for site in poscar.structure if site.species_string == m]
    s_z_coordinates = [site.z for site in poscar.structure if site.species_string == "S"]
    uppermost_au_z = max(au_z_coordinates)
    lowest_s_z = min(s_z_coordinates)
    
    # Take the average of the uppermost Au atom and the lowest S atom to find the interface
    interface_z = (uppermost_au_z + lowest_s_z) / 2
    
    # Convert the interface z-coordinate from fractional to Cartesian coordinates
    interface_z_cart = interface_z# * poscar.structure.lattice.matrix[2][2]
    
    # Convert the interface z-coordinate from Cartesian to grid coordinates
    interface_index = int(interface_z_cart / chgcar.structure.lattice.matrix[2][2] * charge_density.shape[-1])
    
    
    
    # Calculate the charge above and below the interface
    charge_below_interface = np.sum(charge_density[:, :, :interface_index])
    charge_above_interface = np.sum(charge_density[:, :, interface_index:])
    
    # Normaqlize the charge above and below the interface
    
    charge_below_interface =  charge_below_interface/volume
    charge_above_interface =  charge_above_interface/volume
    
    print("---------------------------------")
    print(f"Charge below interface {m}: {charge_below_interface}")
    print(f"Charge above interface {m}: {charge_above_interface}")
   
    print("---------------------------------")

    total_charge_loss = total_charge - (total_charge_metal + total_charge_mos2)  
    
    chgtran_to_metal = charge_below_interface - total_charge_metal  
    chgtran_to_mos2 = charge_above_interface  - total_charge_mos2 
    print(f"Charge to {m} : {chgtran_to_metal}")
    print(f"Charge to mos2 : {chgtran_to_mos2}")
    
    tran1= 
    
    
    print("---------------------------------")

    print(f"Charge loss upon interface formation: {total_charge_loss}")

    print("---------------------------------")
    print("---------------------------------")
    print("---------------------------------")
    
    
    
    
  
    
    
    
    
    
    