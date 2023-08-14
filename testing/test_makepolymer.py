import hydrophobicity_explorer as hx
from hydrophobicity_explorer import MakePolymer
import pytest, os
from rdkit import Chem


# def test_of_Polymer_creation():
#     init = "Hydroxyl" # "O" 
#     term = "Benzyl_alcohol" # *Cc1ccccc1O
#     monomer = "Methylacrylate" # CC(C(=O)OC)C
#     n = 3

#     for boolean in [True, False]:
#         Returned_Stuff = MakePolymer.createPolymerObj(init, n, monomer, term, test = boolean)
#         if boolean:
#             assert len(Returned_Stuff) == 2
#             smiles, POL = Returned_Stuff
#             assert smiles == init + monomer + term
#         else:
#             POL = Returned_Stuff
        
#         assert POL.smiles == init + n * monomer + term

### Well Well Well, look at how the use of global variables came back to bite me!

def test_reading_from_files():
    
    n = 3
    smiles = "O" + n * "Cc1ccccc1O" + "CC(C(=O)OC)C"

    pol = Chem.MolFromSmiles(smiles)

    name="tmp.mol"
    Chem.MolToMolFile(pol, name)

    POL = MakePolymer.read_pol(name, n)
    os.remove(name)

    canon_read = Chem.CanonSmiles(POL.smiles)
    canon_original = Chem.CanonSmiles(smiles)

    
    assert canon_original == canon_read

def test_draw_pol():
    smiles = "O" + "Cc1ccccc1O" + "CC(C(=O)OC)C"
    pol = Chem.MolFromSmiles(smiles)
    name = "tmp.png"
    MakePolymer.drawPol(pol, drawName=name)
    pol2 = Chem.MolFromPNGFile(name)
    os.remove(name)

    inch = Chem.MolToInchi(pol)
    inch2 = Chem.MolToInchi(pol2)

    assert inch == inch2

