"""
This script converts phosphorylated residues in a PDB file to their canonical forms.
Phosphorylated residues are represented by specific three-letter codes (SEP, TPO, PTR)
and are converted to their canonical forms (SER, THR, TYR). The script also removes
the phosphate atoms associated with these residues.
-> This means that we can simply hydrogenate the residues and use their L-configurations in Opsrey
"""


from Bio.PDB import PDBParser, PDBIO, Select

# Mapping phospho to canonical residue names
PHOSPHO_TO_CANONICAL = {
    "SEP": "SER",
    "TPO": "THR",
    "PTR": "TYR"
}

# Atoms associated with phosphate groups
PHOSPHO_ATOMS = {"P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P"}

class PhosphoCleaner(Select):
    def accept_residue(self, residue):
        return True  # Keep all residues for now

    def accept_atom(self, atom):
        resname = atom.get_parent().get_resname()
        atom_name = atom.get_name()
        # Remove phosphate atoms from any residue that could have them
        if (
            resname in PHOSPHO_TO_CANONICAL.keys() or
            resname in PHOSPHO_TO_CANONICAL.values()
        ) and atom_name in PHOSPHO_ATOMS:
            return False
        return True

def convert_phospho_residues(structure):
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                if resname in PHOSPHO_TO_CANONICAL:
                    # Set new canonical residue name
                    residue.resname = PHOSPHO_TO_CANONICAL[resname]
                    # Change HETATM to ATOM by setting id[0] to " "
                    residue.id = (" ", residue.id[1], residue.id[2])

def renumber_atoms(structure):
    serial = 1
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.serial_number = serial
                    serial += 1


def main(input_pdb, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)

    convert_phospho_residues(structure)
    renumber_atoms(structure)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, select=PhosphoCleaner())

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python phosphorylated_to_canonical.py input.pdb output.pdb")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])