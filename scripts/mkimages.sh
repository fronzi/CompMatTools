#!/bin/bash

# Script to batch render VASP wavefunction isosurfaces using VESTA
# Assumes VESTA is installed at /Applications/VESTA.app on macOS

# Loop over all wavefunction files matching pattern
for file in *.vasp; do
    echo "here"
    echo $file 
    echo $PWD
        echo $file
        # Extract base name (without extension)
        echo $file
        echo $PWD
        base="${file%.vasp}"
        echo "${base}"

        echo "→ Rendering $file ..."
        # Run VESTA in new instance with export image option
        open -n -a VESTA --args -open "$PWD/$file" -rotate_x 90 -export_img "${base}.png"

        # Wait briefly to allow VESTA to start up (tune if needed)
        sleep 2
done

echo "✓ All wavefunction images rendered."
