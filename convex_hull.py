import sys
import os

# Prefer adding the CCKStar package directory directly so internal imports like
# `from Make_Convex_Hull import ...` work when importing modules.
# some interactive environments (REPL, Python -c) don't define __file__;
# fall back to the current working directory in that case.
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


# ensure output folder exists so SCOPE can write hull PDBs
outfolder = "pdb_hulls"
os.makedirs(outfolder, exist_ok=True)

contacts, interchain = SCOPE(
    "2LOB_Model_One.clean.pdb",
    "pdb_hulls",
    "C",                      # design chain ID
    ['TRP'],          # design_AA_type (example)
    True,                      # savePDB
    "L",                      # design_chirality
    []                         # fixed_identity
)
print(contacts, interchain)