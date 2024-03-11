import os
import csv
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import BSPlotter

# Function to extract band gap

def extract_band_gap(vasprun_path):
    try:
        vasprun = Vasprun(vasprun_path)
        band_structure = vasprun.get_band_structure()
        band_gap = band_structure.get_band_gap()
        return band_gap
    except Exception as e:
        print(f"Error parsing {vasprun_path}: {e}")
        return None  # or appropriate default value/error indicator




# Directory containing vasprun.xml files
working_dir = '.'  # current directory; modify as needed

# CSV file to write the band gaps
csv_filename = 'band_gaps.csv'

# Find vasprun.xml files and extract band gaps
band_gaps = []


# In the loop where band gaps are extracted:
for subdir, dirs, files in os.walk(working_dir):
    for file in files:
        if file == 'vasprun.xml':
            file_path = os.path.join(subdir, file)
            band_gap_info = extract_band_gap(file_path)
            if band_gap_info is not None:
                band_gaps.append([subdir, band_gap_info['energy'], band_gap_info['direct']])



# Write band gaps to CSV file
with open(csv_filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Directory', 'Band Gap (eV)', 'Direct'])
    writer.writerows(band_gaps)

print(f"Band gap data written to {csv_filename}")

