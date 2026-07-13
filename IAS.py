import sys
from pdbfixer import PDBFixer
import osprey

osprey.start()

if len(sys.argv) < 2:
    print("Usage: python IAS.py <file_name> [protein_flexible_residue] [ligand_residue]")
    sys.exit(1)

def normalize_residue(value, default="", empty_values=None):
    if value is None:
        return default
    value = str(value).strip()
    if value == "" or value.upper() in (empty_values or set()):
        return default
    return value


file_name = sys.argv[1]
protein_flexible_residue = normalize_residue(
    sys.argv[2] if len(sys.argv) > 2 else "",
    default="",
    empty_values={"NONE", "NULL", "N/A", "A0"},
)
target_res = normalize_residue(
    sys.argv[3] if len(sys.argv) > 3 else "",
    default="",
    empty_values={"NONE", "NULL", "N/A", "B0"},
)

print("Using pdb file:", file_name)
print("Protein flexible residue:", protein_flexible_residue or "<none>")
print("Ligand residue:", target_res or "<none>")

fixer = PDBFixer(filename=file_name)

residue_count_chain_a = 0
for chain in fixer.topology.chains():
    if chain.id == 'A':
        for residue in chain.residues():
            residue_count_chain_a += 1

residue_count_chain_b = 0
for chain in fixer.topology.chains():
    if chain.id == 'B':
        for residue in chain.residues():
            residue_count_chain_b += 1

print(f"Number of residues in Chain A: {residue_count_chain_a}")
print(f"Number of residues in Chain B: {residue_count_chain_b}")

chain_a_residues = ["A" + str(i) for i in range(1, residue_count_chain_a + 1)]
chain_b_residues = ["B" + str(i) for i in range(1, residue_count_chain_b + 1)]

prot_flex = [protein_flexible_residue] if protein_flexible_residue else []
lig_flex = [target_res] if target_res else []
mut_residues = [target_res] if target_res else []

mutated_residue_label = mut_residues[0] if mut_residues else "<none>"
print("Mutated residue:", mutated_residue_label)
print("Protein Flexibilities:")
print(prot_flex)
print("Ligand Flexibilities:")
print(lig_flex)

ffparams = osprey.ForcefieldParams()
mol = osprey.readPdb('./' + file_name)
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

protein = osprey.Strand(
    mol,
    templateLib=templateLib,
    residues=[chain_a_residues[0], chain_a_residues[-1]]
)

for res in prot_flex:
    protein.flexibility[res] \
        .setLibraryRotamers(osprey.WILD_TYPE) \
        .addWildTypeRotamers() \
        .setContinuous()

ligand = osprey.Strand(
    mol,
    templateLib=templateLib,
    residues=[chain_b_residues[0], chain_b_residues[-1]]
)

for res in lig_flex:
    ligand.flexibility[res] \
        .setLibraryRotamers(osprey.WILD_TYPE) \
        .addWildTypeRotamers() \
        .setContinuous()

three_letter_amino_acid_codes = [
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL"
]

for res in mut_residues:
    ligand.flexibility[res] \
        .setLibraryRotamers(three_letter_amino_acid_codes) \
        .addWildTypeRotamers() \
        .setContinuous()

trflex = osprey.c.confspace.StrandFlex.TranslateRotate(10, 2.0)
ligandForConf = [ligand, trflex]

proteinConfSpace = osprey.ConfSpace([protein])
ligandConfSpace = osprey.ConfSpace([ligandForConf])
complexConfSpace = osprey.ConfSpace([protein, ligandForConf])

parallelism = osprey.Parallelism(cpuCores=4)
ecalc = osprey.EnergyCalculator(complexConfSpace, ffparams, parallelism=parallelism)

eps = 0.683
print("Using epsilon:", eps)

run_tag = f"{file_name}_{mutated_residue_label}_flex_{protein_flexible_residue or 'none'}"

kstar = osprey.KStar(
    proteinConfSpace,
    ligandConfSpace,
    complexConfSpace,
    epsilon=eps,
    writeSequencesToFile=f"{run_tag}.tsv",
    maxSimultaneousMutations=100
)

for info in kstar.confSpaceInfos():
    eref = osprey.ReferenceEnergies(info.confSpace, ecalc)
    info.confEcalc = osprey.ConfEnergyCalculator(
        info.confSpace,
        ecalc,
        referenceEnergies=eref
    )

    emat = osprey.EnergyMatrix(
        info.confEcalc,
        cacheFile=f"emat.{run_tag}.{info.id}.dat"
    )

    def makePfunc(rcs, confEcalc=info.confEcalc, emat=emat):
        return osprey.PartitionFunction(
            confEcalc,
            osprey.AStarTraditional(emat, rcs, showProgress=False),
            osprey.AStarTraditional(emat, rcs, showProgress=False),
            rcs
        )

    info.pfuncFactory = osprey.KStar.PfuncFactory(makePfunc)

scoredSequences = kstar.run(ecalc.tasks)

analyzer = osprey.SequenceAnalyzer(kstar)

print()
print("Using results now.")
print()

for scoredSequence in scoredSequences:
    print("result:")
    print("\tsequence: %s" % scoredSequence.sequence)
    print("\tK* score: %s" % scoredSequence.score)

    numConfs = 10
    analysis = analyzer.analyze(scoredSequence.sequence, numConfs)
    print(analysis)

    analysis.writePdb(
        f"seq.{run_tag}.{scoredSequence.sequence}.pdb",
        f"Top {numConfs} conformations for sequence {scoredSequence.sequence}"
    )

print("Done!")