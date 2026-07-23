import os
import sys

from Find_Doublets import SCOPE, rank_flex_overlap


def run_convex_hull(file_name, outfolder="pdb_hulls"):
    """Run the convex-hull contact scan and return contacts + interchain."""
    os.makedirs(outfolder, exist_ok=True)

    contacts, interchain = SCOPE(
        file_name,
        outfolder,
        "B",                      # design chain ID
        ["TRP"],                 # design_AA_type (example)
        True,                     # savePDB
        "L",                      # design_chirality
        []                        # fixed_identity
    )
    return contacts, interchain


def build_experiment_rows(interchain):
    """Translate the interchain residue lists into IAS experiment rows."""
    rows = [{
        "mode": "baseline",
        "protein_flex_residue": "",
        "ligand_flex_residue": "",
    }]

    for ligand_idx, nearby_proteins in enumerate(interchain):
        ligand_residue = f"B{ligand_idx + 1}"
        for protein_idx in nearby_proteins:
            protein_residue = f"A{protein_idx + 1}"
            rows.append({
                "mode": "single-protein-flex",
                "protein_flex_residue": protein_residue,
                "ligand_flex_residue": ligand_residue,
            })

    return rows


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 convex_hull.py <file_name>")
        sys.exit(1)

    file_name = sys.argv[1]
    contacts, interchain = run_convex_hull(file_name)

    print(contacts)
    print(interchain)

    for obj in interchain:
        print(obj)

    # # optional: order the flexible residues by volume overlap with design chain hulls
    # # this is useful for prioritizing flexible residues over a large search space
    # # returns: {design res, target res : cubic angstrom overlap}
    # flex_order = rank_flex_overlap(
    #     "B",
    #     "A",
    #     contacts,
    #     interchain,
    #     "pdb_hulls",
    # )
    # print("Printing Flex Order:")
    # print(flex_order)


if __name__ == "__main__":
    main()

