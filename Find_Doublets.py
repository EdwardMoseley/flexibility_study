import sys
from Bio.PDB import PDBParser, Select, PDBIO
import warnings
import numpy as np
from collections import OrderedDict
import os
from scipy.spatial import ConvexHull as ScipyConvexHull
from scipy.spatial import Delaunay
import pyvista as pv
import vtk


from Make_Convex_Hull import make_convex_hull

# hide warnings - biopython throws for a lot of atom names, but is ok
warnings.simplefilter('ignore')

# vtk throws warnings for non-overlapping hulls, which we will be working on
vtk.vtkObject.GlobalWarningDisplayOff()


# class for saving backbone-only PDB files
class BBAtoms(Select):
    def accept_atom(self, atom):
        if atom.name in ["C", "CA", "HA", "N"]:
            return True
        else:
            return False


# function that strips all sidechains from a PDB
# only keeps C, CA, HA, N
def make_backbone_PDB(pdb_filename: str, designChirality: str):
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("complex", pdb_filename)
    model = structure[0]

    # make sure backbone is protonated (need for alignment) and adjust GLY atom naming
    done_target = False
    for chain in model:
        for residue in chain:
            resname = residue.get_resname()
            atoms = residue.get_atoms()
            atom_names = [a.name for a in atoms]
            if resname == "GLY":
                if 'HA2' not in atom_names or 'HA3' not in atom_names:
                    print("ERROR! Residue %s %s is missing a backbone proton" % (resname, residue.id[1]),
                          file=sys.stderr)
                    exit()
                for atom in residue:
                    if designChirality == 'D' and done_target:
                        if atom.fullname == ' HA2':
                            atom.fullname = ' HA '
                            atom.name = "HA"
                    elif not done_target or designChirality == 'L':
                        if atom.fullname == ' HA3':
                            atom.fullname = ' HA '
                            atom.name = "HA"
            else:
                if 'HA' not in atom_names:
                    print("ERROR! Residue %s %s is missing a backbone proton" % (resname, residue.id[1]),
                          file=sys.stderr)
                    exit()
        done_target = True

    # strip sidechains
    io = PDBIO()
    io.set_structure(structure)
    new_filename = pdb_filename.split(".")[0] + "-bbonly.pdb"
    io.save(new_filename, BBAtoms())
    print("Saved backbone only PDB to %s" % new_filename)

    return new_filename


def origin_anchors(target_N: list, target_CA: list, target_C: list, target_HA: list):
    x_diff = 0 - target_CA[0]
    y_diff = 0 - target_CA[1]
    z_diff = 0 - target_CA[2]

    target_N[0] = target_N[0] + x_diff
    target_N[1] = target_N[1] + y_diff
    target_N[2] = target_N[2] + z_diff

    target_C[0] = target_C[0] + x_diff
    target_C[1] = target_C[1] + y_diff
    target_C[2] = target_C[2] + z_diff

    target_HA[0] = target_HA[0] + x_diff
    target_HA[1] = target_HA[1] + y_diff
    target_HA[2] = target_HA[2] + z_diff

    return target_N, target_C, target_HA


# class of CH w/ anchor coords (arbitrarily made relative to VALp rotamer of respective chirality)
class ConvexHull:
    def __init__(self, resid, N_back, CA_back, C_back, HA_back, hull_coords):
        self.resid = resid
        self.N_back = N_back
        self.CA_back = CA_back
        self.C_back = C_back
        self.HA_back = HA_back
        self.hull_coords = hull_coords

    def translateCH(self, target_CA: list):
        x_diff = target_CA[0] - self.CA_back[0]
        y_diff = target_CA[1] - self.CA_back[1]
        z_diff = target_CA[2] - self.CA_back[2]

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom != "hull_coords":
                    x_new = x_diff + xyz[0]
                    y_new = y_diff + xyz[1]
                    z_new = z_diff + xyz[2]
                    new_coord = [x_new, y_new, z_new]
                    self.__setattr__(atom, new_coord)
                else:
                    new_hull = []
                    for item in xyz:
                        x_new = x_diff + item[0]
                        y_new = y_diff + item[1]
                        z_new = z_diff + item[2]
                        new_hull.append([x_new, y_new, z_new])
                    self.__setattr__(atom, new_hull)

    def moveCH(self, target_N: list, target_CA: list, target_C: list, target_HA: list):

        # we must be at the origin (CH and target) to find the rotation matrix
        self.translateCH([0, 0, 0])

        origin_target_N, origin_target_C, origin_target_HA = origin_anchors(target_N, target_CA, target_C, target_HA)

        # find rot matrix
        V_matrix = np.array([self.N_back, self.C_back, self.HA_back])
        V_matrix_formatted = np.transpose(V_matrix)
        V_matrix_inv = np.linalg.inv(V_matrix_formatted)

        V_prime_matrix = np.array([origin_target_N, origin_target_C, origin_target_HA])
        V_prime_formatted = np.transpose(V_prime_matrix)
        rotation_matrix = np.matmul(V_prime_formatted, V_matrix_inv)

        # rotate each atom
        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom != "hull_coords":
                    xyz_matrix = np.array(xyz)
                    new_coords = np.matmul(rotation_matrix, xyz_matrix).tolist()
                    self.__setattr__(atom, new_coords)
                else:
                    new_hull = []
                    for item in xyz:
                        xyz_matrix = np.array(item)
                        new_coords = np.matmul(rotation_matrix, xyz_matrix).tolist()
                        new_hull.append(new_coords)
                    self.__setattr__(atom, new_hull)

        # move back to target CA
        self.translateCH(target_CA)

    def print_pdb(self, filename, chainID):
        f = open(filename, "w")
        atom_num = 1
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "CHU".ljust(3)
        pdb_info[4] = chainID.rjust(1)
        pdb_info[5] = "1".rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom in ["CA_back", "N_back", "C_back", "HA_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                    pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                    pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                    pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                    f.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                             pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                             pdb_info[8]))
                    atom_num += 1
                    pdb_info[1] = str(atom_num).rjust(5)
                else:
                    for item in xyz:
                        pdb_info[2] = "C".center(4)
                        pdb_info[6] = str('%8.3f' % (item[0])).rjust(8)
                        pdb_info[7] = str('%8.3f' % (item[1])).rjust(8)
                        pdb_info[8] = str('%8.3f' % (item[2])).rjust(8)
                        f.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                                 pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                                 pdb_info[8]))
                        atom_num += 1
                        pdb_info[1] = str(atom_num).rjust(5)
        f.write("TER\n")
        f.close()


# function for placing CH for each residue on the specified chain
def insert_hulls(bb_pdb: str, target_chain: str, AA_types: list, designChirality: str, FixedIdentity: list):
    all_hulls = []

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("complex", bb_pdb)
    model = structure[0]

    for chain in model:
        if chain.id == target_chain:
            for residue in chain:
                res_name = residue.get_resname()
                res_num = residue.id[1]
                N_target = []
                CA_target = []
                C_target = []
                HA_target = []
                for atom in residue:
                    if atom.name == "N":
                        N_target = list(atom.coord)
                    elif atom.name == "CA":
                        CA_target = list(atom.coord)
                    elif atom.name == "C":
                        C_target = list(atom.coord)
                    elif atom.name == "HA":
                        HA_target = list(atom.coord)
                if not N_target:
                    print("ERROR! Didn't find N anchor for chain %s residue %s %s" % (target_chain,
                                                                                      res_name,
                                                                                      res_num), file=sys.stderr)
                    exit()
                if not CA_target:
                    print("ERROR! Didn't find CA anchor for chain %s residue %s %s" % (target_chain,
                                                                                       res_name,
                                                                                       res_num), file=sys.stderr)
                    exit()
                if not C_target:
                    print("ERROR! Didn't find C anchor for chain %s residue %s %s" % (target_chain,
                                                                                      res_name,
                                                                                      res_num), file=sys.stderr)
                    exit()
                if not HA_target:
                    print("ERROR! Didn't find HA anchor for chain %s residue %s %s" % (target_chain,
                                                                                       res_name,
                                                                                       res_num), file=sys.stderr)
                    exit()

                if AA_types[0] == 'wt' or res_num in FixedIdentity:
                    if res_num in FixedIdentity:
                        print("Making WT hull for fixed residue %s%s" % (res_name, res_num))
                    if res_name in ['GLY', 'PRO', 'ALA']:
                        # approximate residues with too few atoms for CH with VAL
                        convex_hull_coords = make_convex_hull(["VAL"], "", False, designChirality)
                    elif res_name == 'HIS':
                        convex_hull_coords = make_convex_hull(["HID", "HIE", "HIP"], "", False, designChirality)
                    elif res_name not in ['GLY', 'PRO', 'ALA']:
                        convex_hull_coords = make_convex_hull([res_name], "", False, designChirality)
                else:
                    convex_hull_coords = make_convex_hull(AA_types, "", False, designChirality)

                if designChirality == 'L':
                    newCH = ConvexHull(res_num, [0.605000000000004, 0.8560000000000016, 1.016], [0, 0, 0],
                                   [-0.7420000000000044, 0.8540000000000028, -1.0280000000000005],
                                   [0.8059999999999974, -0.5640000000000001, -0.4930000000000003],
                                   convex_hull_coords)
                elif designChirality == 'D':
                    newCH = ConvexHull(res_num, [0.605000000000004, 0.8560000000000016, -1.016], [0, 0, 0],
                                       [-0.7420000000000044, 0.8540000000000028, 1.0280000000000005],
                                       [0.8059999999999974, -0.5640000000000001, 0.4930000000000003],
                                       convex_hull_coords)
                newCH.moveCH(N_target, CA_target, C_target, HA_target)
                all_hulls.append(newCH)

    return all_hulls


def print_hulls(all_hulls: list, outfolder, target_chain):

    if len(all_hulls) <= 52:
        chain_id = ord('A')
        filenames = []
        for h in all_hulls:
            res_id = h.resid
            filename = ("%s/Chain%sRes%s.pdb" % (outfolder, target_chain, res_id))
            filenames.append(filename)
            h.print_pdb(filename, chr(chain_id))
            chain_id += 1

        all_pdb_hull = os.path.join(outfolder, ("Chain%s_all_hulls.pdb" % target_chain))
        print("Saving all CH to %s" % all_pdb_hull)
        with open(all_pdb_hull, 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

    elif len(all_hulls) > 52:
        print("WARNING: Number of residues exceeds number of unique chain IDs, so assigning all as chain A")
        chain_id = "A"
        filenames = []
        for h in all_hulls:
            res_id = h.resid
            filename = ("%s/Chain%sRes%s.pdb" % (outfolder, target_chain, res_id))
            filenames.append(filename)
            h.print_pdb(filename, chain_id)

        all_pdb_hull = os.path.join(outfolder, ("Chain%s_all_hulls.pdb" % target_chain))
        print("Saving all CH to %s" % all_pdb_hull)
        with open(all_pdb_hull, 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)


def polydata_from_hull(pts):
    hull = ScipyConvexHull(pts)
    faces_idx = hull.simplices  # (M, 3) triangle indices

    # PyVista wants a flat [3, i, j, k, 3, i, j, k, ...]
    faces_pv = np.hstack(np.c_[np.full((len(faces_idx), 1), 3, int), faces_idx]).astype(np.int64).ravel()
    mesh = pv.PolyData(pts, faces_pv).triangulate().clean()
    return mesh


# compute volume overlap of two hulls
def find_volume_overlap(hull1, hull2):

    hull1_mesh = polydata_from_hull(hull1)
    hull2_mesh = polydata_from_hull(hull2)

    inter = hull1_mesh.boolean_intersection(hull2_mesh)

    if inter.n_cells > 0:
        return inter.volume

    return 0.0


def find_intrachain_intersects(hulls: list):
    all_doublets = list()
    have_doublet = []

    for h1 in range(0, len(hulls)):
        res1num = h1 + 1
        res1hull = hulls[h1].hull_coords
        for h2 in range(h1 + 1, len(hulls)):
            res2num = h2 + 1
            res2hull = hulls[h2].hull_coords

            overlap = find_volume_overlap(res1hull, res2hull)

            if overlap > 0.0:
                doublet = {res1num, res2num}
                all_doublets.append(doublet)
                have_doublet.append(res1num)
                have_doublet.append(res2num)

    print("Found Doublets:")
    print(all_doublets)

    print("Found Islands:")
    have_doublet.sort()
    filter_doublet = list(set(have_doublet))
    islands = [{i} for i in range(1, len(hulls)+1) if i not in filter_doublet]

    if len(islands) == 0:
        print("None")
    elif len(islands) != 0:
        print(islands)

    return all_doublets, islands


def find_interchain_intersects(hull1: list, hull2: list, chainIDs):
    all_intersects = []

    # for each hull in h1, find all intersects across all hulls in h2
    for h1 in hull1:
        intersects = []
        for h2 in hull2:

            volume_overlap = find_volume_overlap(h1.hull_coords, h2.hull_coords)
            if volume_overlap > 0.0:
                intersects.append(h2.resid)

        print("Chain %s Residue %s intersects with Chain %s residue(s) %s" % (chainIDs[0], h1.resid,
                                                                              chainIDs[1], intersects))

        if intersects:
            all_intersects.append(intersects)

    return all_intersects


# find flexible residues defined by CH intersects between two chains
def SCOPE(pdb_name, outfolder, designID: str, design_AA_type: list, savePDB: bool, design_chirality: str,
          fixed_identity: list):
    print("--Now finding intra and inter chain contacts for singlechain design space--")
    bbpdb = make_backbone_PDB(pdb_name, design_chirality)

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("complex", bbpdb)
    model = structure[0]

    chain_names = []

    for chain in model:
        chain_names.append(chain.id)

    if chain_names[0] != designID:
        chain_names = [chain_names[1], chain_names[0]]

    print("--Making %s-space CH for design Chain %s with %s mutants: %s--" % (design_chirality, designID,
                                                                             len(design_AA_type), design_AA_type))
    design_hulls = insert_hulls(bbpdb, chain_names[0], design_AA_type, design_chirality, fixed_identity)

    print("--Making CH for target Chain--")
    target_hulls = insert_hulls(bbpdb, chain_names[1], ['wt'], 'L', [])

    print("--Finding intrachain contacts for design Chain %s--" % designID)
    lig_doublets, lig_islands = find_intrachain_intersects(design_hulls)
    all_contacts = lig_doublets + lig_islands

    if savePDB:
        print_hulls(design_hulls, outfolder, chain_names[0])
        print_hulls(target_hulls, outfolder, chain_names[1])

    print("--Finding interchain contacts relative to design Chain %s--" % chain_names[0])
    interchain_intersect = find_interchain_intersects(design_hulls, target_hulls, chain_names)

    return all_contacts, interchain_intersect

def pdb_to_coords(pdb_name: str):

    xyz_coords = []

    f = open(pdb_name, "r")
    for line in f:
        if "TER" in line:
            continue
        coords = line.split('     ')[2]
        coords_format = ""
        counter = -1
        for i in coords:
            if i == ' ' and counter == 0:
                coords_format += ","
                counter += 1
            if i in ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ".", "-"]:
                counter = 0
                coords_format += i

        all_coords = coords_format.split(',')
        atom_coords = []
        for i in all_coords:
            atom_coords.append(float(i))

        xyz_coords.append(atom_coords)

    return xyz_coords


def rank_flex_overlap(designID: str, targetID: str, intrachain_contacts: list, interchain_contacts, hull_folder: str):

    ordered_flex = {}
    hulls_loc = hull_folder + "/Chain%sRes%s.pdb"
    design_hulls = {}
    target_hulls = {}

    print("Ordering interchain contacts by volume overlap")

    # get the CH coords for all design chain residues
    intra_residues = {x for pair in intrachain_contacts for x in pair}
    for res in intra_residues:
        hull_loc = hulls_loc % (designID, res)
        coords = pdb_to_coords(hull_loc)
        design_hulls[res] = coords

    # get the CH for only flex target residues
    inter_residues = {x for pair in interchain_contacts for x in pair}
    for res in inter_residues:
        hull_loc = hulls_loc % (targetID, res)
        coords = pdb_to_coords(hull_loc)
        target_hulls[res] = coords

    # compute mesh overlap using vtk boolean
    for tar_res, tar_hull in target_hulls.items():
        for des_res, des_hull in design_hulls.items():
            volume_overlap = find_volume_overlap(tar_hull, des_hull)
            if volume_overlap > 0.0:
                overlap_pair = (des_res, tar_res)
                ordered_flex[overlap_pair] = volume_overlap

    # order by cubic angstrom overlap
    ordered_flex = dict(sorted(ordered_flex.items(), key=lambda item: item[1], reverse=True))
    return ordered_flex



