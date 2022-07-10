import rdkit, argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from random import choices
from smiles import monomer_dict
from MakePolymer import validate_end_group, inator_smi_lookup, add_inator_smiles

def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type=int, help="Total number of monomers desired in final product.")
    parser.add_argument("-i", "--initiator", type=str, default="Hydrogen", help="Initiator Key from initiator dict or SMILES. Defaults to Hydrogen.")
    parser.add_argument("-m", type=str, nargs='*',
            help="Space-separated list of monomer keys from smiles.py or a custom smiles string for a monomer unit. The head of the monomer must be to the left and the attachment points must be at the ends of the string. Integer coefficients can be used to add bias towards a given monomer. For example, 2 A B means A is twice as likely to be used as B when selecting the next monomer.")
    parser.add_argument("-t", "--terminator", type=str, default="Hydrogen", help="Terminator key taken from initiator dict or SMILES. Defaults to Hydrogen.")
    parser.add_argument("-f", type=str, help="Filename prefix you would like to use. The number of monomers will be added automatically.")
    args = parser.parse_args()
    return args

def prepFilename(filename, n):
    split = filename.split(".") #split off file extention in case provided.
    name = f"{split[0]}_{n}.mol"
    return name

def makePolymerBody(weighted_monomer_list, n):
    #first, let us expand the provided coefs to actually give a weighted list.
    #start with empty repeat unit and concatonate stuff we find in the list.
    expanded_list = []
    monomer_weights = []
    #ommission of a coeficient implies 1 copy
    repeat_coef = 1
    for element in weighted_monomer_list:
        try:
            repeat_coef = int(element) #is this a string of an integer?
        except:
            monomer_weights.append(repeat_coef)
            expanded_list.append(element)
            repeat_coef = 1 #reset coef.

    #now that we have an explicit list, we can pick at random from the list n times and make the polymer body.
    body_list = choices(expanded_list, weights=monomer_weights, k=n)

    #now merge list into smiles.
    smiles = "".join(body_list)
    print(f'body smiles is:\n{smiles}')

    return smiles

def main():
    args = getArgs()

    #get proper file name.
    file_name = prepFilename(args.f, args.n)

    #first get inators sorted out:
    #get their smiles
    init, term = inator_smi_lookup(args.initiator, args.terminator)
    #make sure they are formatted the right way.
    init = validate_end_group(init, Init=True)
    term = validate_end_group(term, Term=True)

    #replace any dict keys with corresponding smiles.
    deciphered_dict_keys = [monomer_dict[x] if x in monomer_dict else x for x in args.m]

    #now we need to generate polymer body
    polymer_body_smiles = makePolymerBody(deciphered_dict_keys, args.n)

    #now we need to attatch the end groups:
    total_smiles = add_inator_smiles(polymer_body_smiles, init, term)
    print("Finished adding end groups. Beginning optimization.")

    mol = Chem.MolFromSmiles(total_smiles)
    
    Chem.SanitizeMol(mol)
    #opt steps
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_h, useRandomCoords=True)
    AllChem.MMFFOptimizeMolecule(mol_h, maxIters=5000)
    #maybe this number of itterations should be specified with cli arguments (give option).
    
    Chem.MolToMolFile(mol_h, file_name)

    print(f"Done. Saved to {file_name}")

if __name__ == "__main__":
    main()
