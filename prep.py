## Imports and Install Notes
# conda install -c conda-forge pdbfixer
from pdbfixer import PDBFixer

# conda install -c conda-forge openmm-setup
from openmm.app import PDBFile
import warnings
import sys
import os
warnings.simplefilter('ignore')

# Parse command line arguments
if len(sys.argv) < 2:
    print("Usage: python prep.py <file_name>")
    sys.exit(1)

file_name = sys.argv[1]

print("Working on File:\n" + file_name)

#Load the PBD file
fixer = PDBFixer(filename=file_name)

# Add missing residues and atoms
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(pH=7.0)

# Adjust File Name to Account for Cleaning
base, ext = os.path.splitext(file_name)
cleaned_file_name = base + ".clean.pdb"

# Save
# Save cleaned PDB
with open(cleaned_file_name, 'w') as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)