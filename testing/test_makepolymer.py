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

class TestPolymerIO:
    n = 3
    smiles = "O" + n * "Cc1ccccc1O" + "CC(C(=O)OC)C"

    def test_reading_from_files(self):
        pol = Chem.MolFromSmiles(self.smiles)

        name="tmp.mol"
        Chem.MolToMolFile(pol, name)

        POL = MakePolymer.read_pol(name, self.n)
        os.remove(name)

        canon_read = Chem.CanonSmiles(POL.smiles)
        canon_original = Chem.CanonSmiles(self.smiles)

        
        assert canon_original == canon_read

    def test_draw_single_pol(self):
        pol = Chem.MolFromSmiles(self.smiles)
        name = "tmp.png"
        MakePolymer.drawPol(pol, drawName=name)
        pol2 = Chem.MolFromPNGFile(name)
        os.remove(name)

        inch = Chem.MolToInchi(pol)
        inch2 = Chem.MolToInchi(pol2)

        assert inch == inch2

    def test_draw_grid_of_pols(self):
        POL = MakePolymer.Polymer(self.n, self.smiles)
        pol = [POL]

        name = "tmp.png"
        MakePolymer.drawPol(pol, drawName=name)
        
        # We can't read from a grid image directly, so just check if file exists
        assert os.path.exists(name)
        os.remove(name)


    def test_write_pol_to_file(self):
        POL = MakePolymer.Polymer(self.n, self.smiles)

        nConfs = 5
        pol_list = []
        for _ in range(nConfs):
            pol_h = MakePolymer.optPol(POL.flat)
            pol_list.append(pol_h)

        write_name = "tmp.mol"
        MakePolymer.write_pol(write_name, pol_list)

        read_mol = Chem.MolFromMolFile(write_name)
        Chem.RemoveStereochemistry(read_mol) # FLAT won't have any

        os.remove(write_name)

        smi1 = Chem.MolToSmiles(POL.flat)
        smi2 = Chem.MolToSmiles(read_mol)

        smi1 = Chem.CanonSmiles(smi1)
        smi2 = Chem.CanonSmiles(smi2)

        assert smi1 == smi2

def test_make_simple_polymer():
    MakePolymer.main(n = 1, single_monomer = "Ethylene", initiator = "Hydroxy", terminator = "*CC(C(=O)OC)C", calculation = ["MHP", "RG", "MV"], quiet = True)

def test_plotting():
    MakePolymer.main(n = 2, single_monomer = "Ethylene", initiator = "Hydroxy", terminator = "*CC(C(=O)OC)C", calculation = ["MHP", "RG", "MV"], quiet = True, plot = True)
    plot_name = "Size-dependent-stats.png"
    assert os.path.exists(plot_name)
    os.remove(plot_name)

def test_pallendromic_end_piece_filter():
    with pytest.raises(ValueError): # terminator is non-pallendromic and needs to have "*"
        MakePolymer.main(n = 1, single_monomer = "Styrene", initiator = "Hydroxy", terminator = "CC(C(=O)OC)C", calculation = ["MHP", "RG", "MV"], quiet = True)

def test_make_block_polymer():
    MakePolymer.main(n = 2, comonomer_sequence = ["Vinylalcohol", "Vinylchloride"], calculation = ["MHP", "RG", "MV"])

def test_make_random_polymer():
    MakePolymer.main(n = 2, comonomer_sequence = ["Vinylalcohol", "Vinylchloride"], calculation = ["MHP", "RG", "MV"], random = True)