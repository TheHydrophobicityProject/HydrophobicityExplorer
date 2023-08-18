from hydrophobicity_explorer import nb, custom_input_to_mol_file
from rdkit import Chem
import pytest, os, json

def test_nb_write():
    nb.main()
    with open(nb.name, "r") as S:
        read_nb = json.load(S)
    
    os.remove(nb.name)
    assert read_nb == nb.notebook

def test_custom_polymer():
    mol_file = "styrene.mol"
    orig_smi = "CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)"
    custom_input_to_mol_file.main(smiles=orig_smi, file=mol_file)
    mol = Chem.MolFromMolFile(mol_file)
    read_smi = Chem.MolToSmiles(mol)

    os.remove(mol_file)

    orig_smi = Chem.CanonSmiles(orig_smi)
    read_smi = Chem.CanonSmiles(read_smi, useChiral=0)

    assert orig_smi == read_smi