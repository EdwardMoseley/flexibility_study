import sys
from pdbfixer import PDBFixer
import osprey

osprey.start()

if len(sys.argv) < 4:
    print("Usage: python IAS_noproteinflex.py <file_name> <protein_flexible_residue_tag> <target_res>")
    sys.exit(1)

file_name = sys.argv[1]
protein_flexible_residue = sys.argv[2]  # should be A0
target_res = sys.argv[3]                # B1 ... B14

if protein_flexible_residue != "A0":
    raise ValueError(
        f"This no-protein-flexibility script expects protein_flexible_residue='A0', "
        f"but received {protein_flexible_residue}"
    )

print("Using pdb file:", file_name)
print("Protein flexible residue tag:", protein_flexible_residue)
print("Ligand mutation position:", target_res)

fixer = PDBFixer(filename=file_name)

residue_count_chain_a = 0
for chain in fixer.topology.chains():
    if chain.id == "A":
        for residue in chain.residues():
            residue_count_chain_a += 1

residue_count_chain_b = 0
for chain in fixer.topology.chains():
    if chain.id == "B":
        for residue in chain.residues():
            residue_count_chain_b += 1

print(f"Number of residues in Chain A: {residue_count_chain_a}")
print(f"Number of residues in Chain B: {residue_count_chain_b}")

chain_a_residues = ["A" + str(i) for i in range(1, residue_count_chain_a + 1)]
chain_b_residues = ["B" + str(i) for i in range(1, residue_count_chain_b + 1)]

if target_res not in chain_b_residues:
    raise ValueError(
        f"Target residue {target_res} is not present in Chain B residues: "
        f"{chain_b_residues[0]} through {chain_b_residues[-1]}"
    )

prot_flex = []
lig_flex = [target_res]
mut_residues = [target_res]

print("Mutated residue:", mut_residues[0])
print("Protein Flexibilities:")
print(prot_flex)
print("Ligand Flexibilities:")
print(lig_flex)

ffparams = osprey.ForcefieldParams()
mol = osprey.readPdb("./" + file_name)
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

protein = osprey.Strand(
    mol,
    templateLib=templateLib,
    residues=[chain_a_residues[0], chain_a_residues[-1]]
)

print("Skipping protein flexibility because protein_flexible_residue is A0.")

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
    "LEU", "LYS", "MET", "PHE", #"PRO",
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

run_tag = f"{file_name}_{mut_residues[0]}_flex_{protein_flexible_residue}"

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