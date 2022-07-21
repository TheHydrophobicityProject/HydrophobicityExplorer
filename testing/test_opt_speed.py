import rdkit, time
from rdkit import Chem
from rdkit.Chem import AllChem

times = []
smiles = "CC(C(=O)OC)(C)CC(c1ccccc1)CC(C)CC(C)CC(C)CC(C(=O)OC)(C)CC(C(=O)OC)(C)CC(C)CC(C(=O)OC)(C)CC(C)CC(c1ccccc1)CC(C)CC(C)CC(C)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(C)CC(C)CC(C(=O)OC)(C)CC(C)CC(c1ccccc1)CC(C)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(C)CC(C(=O)OC)(C)CC(C)"

def setup(smiles):
    pol = Chem.MolFromSmiles(smiles)
    #check mol
    Chem.SanitizeMol(pol)
    #opt steps
    pol_h = Chem.AddHs(pol)
    return pol_h

def initialOpt(smiles):
    pol_h = setup(smiles)
    #random coords lead to better geometries than using the rules rdkit has. Excluding this kwarg leads to polymers that do not fold properly.
    AllChem.EmbedMolecule(pol_h, useRandomCoords=True) 
    AllChem.MMFFOptimizeMolecule(pol_h, maxIters=2500) 
    return pol_h
    # 24.1528112411499

def optMoreConfs(smiles):
    pol_h = setup(smiles)
    ids = AllChem.EmbedMultipleConfs(pol_h, numConfs=5, useRandomCoords=True, numThreads=0)
    touple_list = AllChem.MMFFOptimizeMoleculeConfs(pol_h, numThreads=0, maxIters=2500) #default 200 iterations.
    return touple_list, pol_h
    # 51.50327501296997 twice the time but 10X the data with 10 confs.
    # 25.218886041641234 when using 5 confs and 200 iters

NUMRUNS = 1
#optimize a fairly difficult molecule several times to see how good each optimization method is.
for i in range(NUMRUNS):
    start = time.time()
    touple_list, pol_h = optMoreConfs(smiles)
    end = time.time()
    duration = end - start
    print(i, duration)
    times.append(duration)

average = sum(times) / len(times)
print(f"{average = }")

print(touple_list)

#Writing all the conformations to a sdf file.
writer = Chem.SDWriter('pol_h_confs.sdf')
for cid in range(pol_h.GetNumConformers()):
    pol_h.SetProp('_Name', f'conformer_{cid}')
    # pol_h.SetProp('ID', f'conformer_{cid}') #Similar method can be used to print number of monomers for plot jobs.
    writer.write(pol_h, confId=cid)

# Calculating avgerage properties
