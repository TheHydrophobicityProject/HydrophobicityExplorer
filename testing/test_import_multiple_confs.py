import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, rdFreeSASA

# Calculating avgerage properties from multiple confs
#read them from sdf file.

suppl = Chem.SDMolSupplier("pol_h_confs.sdf") #itterator that has all mols in the sdf file.
mhps = []

for mol in suppl: #calculate desired property for each molecule in the file.
    radii = Chem.rdFreeSASA.classifyAtoms(mol)
    sasa = Chem.rdFreeSASA.CalcSASA(mol, radii)
    logP = Chem.Descriptors.MolLogP(mol)

    mhp = logP / sasa
    mhps.append(mhp)

avgMHP = sum(mhps) / len(mhps) #get average.

print(f"{avgMHP = }")
print(mhps)
#this is super fast.
