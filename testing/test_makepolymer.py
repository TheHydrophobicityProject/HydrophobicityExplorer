import hydrophobicity_explorer as hx
from hydrophobicity_explorer import MakePolymer
import pytest

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