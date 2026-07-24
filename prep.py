from pdbfixer import PDBFixer
from openmm.app import PDBFile
import tempfile
import warnings
import sys
import os

warnings.simplefilter('ignore')


def strip_hetatm_records(file_name):
    removed = 0
    tmp = tempfile.NamedTemporaryFile('w', suffix='.pdb', delete=False)
    try:
        with open(file_name, 'r') as src:
            for line in src:
                if line.startswith('HETATM'):
                    removed += 1
                    continue
                tmp.write(line)
        tmp.flush()
    finally:
        tmp.close()

    return tmp.name, removed

# Parse command line arguments
if len(sys.argv) < 2:
    print("Usage: python prep.py <file_name>")
    sys.exit(1)

file_name = sys.argv[1]

print("Working on File:\n" + file_name)

if file_name.endswith('.clean.pdb'):
    cleaned_file_name = file_name
else:
    base, ext = os.path.splitext(file_name)
    cleaned_file_name = base + ".clean.pdb"

stripped_file_name, removed_count = strip_hetatm_records(file_name)

try:
    fixer = PDBFixer(filename=stripped_file_name)

    # Add missing residues and atoms
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    tmp_output = cleaned_file_name + ".tmp"
    with open(tmp_output, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    os.replace(tmp_output, cleaned_file_name)
finally:
    if os.path.exists(stripped_file_name):
        os.remove(stripped_file_name)

print(f"Removed {removed_count} HETATM line(s)")
print(f"Wrote cleaned file: {cleaned_file_name}")