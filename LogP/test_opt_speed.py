from distutils.archive_util import make_archive
from MakePolymer import optPol
import time
times = []
smiles = "CC(C(=O)OC)(C)CC(c1ccccc1)CC(C)CC(C)CC(C)CC(C(=O)OC)(C)CC(C(=O)OC)(C)CC(C)CC(C(=O)OC)(C)CC(C)CC(c1ccccc1)CC(C)CC(C)CC(C)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(C)CC(C)CC(C(=O)OC)(C)CC(C)CC(c1ccccc1)CC(C)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(C)CC(C(=O)OC)(C)CC(C)"

#optimize a fairly difficult molecule several times to see how good each optimization method is.
for i in range(10):
    start = time.time()
    mol_h, mol = optPol(smiles)
    end = time.time()
    duration = end - start
    print(i, duration)
    times.append(duration)

average = sum(times) / len(times)
print(f"{average = }")

#single molecule opt with random coords max iters 2500 MMMFFOptimize molecule 24.1528112411499

