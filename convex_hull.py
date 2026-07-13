import sys
import os
from Find_Doublets import SCOPE, rank_flex_overlap

# Prefer adding the CCKStar package directory directly so internal imports like
# `from Make_Convex_Hull import ...` work when importing modules.
# some interactive environments (REPL, Python -c) don't define __file__;
# fall back to the current working directory in that case.

"""
try:
    base_dir = os.path.dirname(__file__)
except NameError:
    base_dir = os.getcwd()
repo_root = os.path.abspath(os.path.join(base_dir, ".."))
cck_dir = os.path.join(repo_root, "src", "main", "python", "CCKStar")
py_src_dir = os.path.join(repo_root, "src", "main", "python")

if cck_dir not in sys.path:
    sys.path.insert(0, cck_dir)

try:
    # import the module as the package author expects (local imports)
    from Find_Doublets import SCOPE
except ModuleNotFoundError:
    # fallback: add the broader python path and try package import
    if py_src_dir not in sys.path:
        sys.path.insert(0, py_src_dir)
    from CCKStar.Find_Doublets import SCOPE
"""

if len(sys.argv) < 2:
    print("Usage: python3 convex_hull.py <file_name>")
    sys.exit(1)

file_name = sys.argv[1]

# ensure output folder exists so SCOPE can write hull PDBs
outfolder = "pdb_hulls"
os.makedirs(outfolder, exist_ok=True)

three_letter_amino_acid_codes = [
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL"
]

## We want to iterate over each of the amino acids in the list 
## and call SCOPE for each one, saving the results to a file.

contacts, interchain = SCOPE(
    file_name,
    "pdb_hulls",
    "B",                      # design chain ID
    ['TRP'],                  # design_AA_type (example)
    False,                      # savePDB
    "L",                      # design_chirality
    []                         # fixed_identity
)
print(contacts)

print(interchain)

# # optional: order the flexible residues by volume overlap with design chain hulls
# # this is useful for prioritizing flexible residues over a large search space
# # returns: {design res, target res : cubic angstrom overlap}
flex_order = rank_flex_overlap('B', 
                               'A', 
                               contacts, 
                               interchain, 
                               'pdb_hulls')

print("Printing Flex Order:")
print(flex_order)