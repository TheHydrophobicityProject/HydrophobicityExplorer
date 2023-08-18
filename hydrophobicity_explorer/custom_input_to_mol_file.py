import rdkit, argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from hydrophobicity_explorer.MakePolymer import getStaticSettings, write_pol, make_One_or_More_Polymers

def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--smiles", type=str, help="Complete smiles string to be converted.")
    parser.add_argument("-m", "--smarts", type=str, help="Complete smarts string to be converted.")
    parser.add_argument("-i", "--inchi", type=str, help="Complete inchi string to be converted.")
    parser.add_argument("-f", "--file", type=str, help="Filename you would like to save to. MOL, PDB and XYZ are acceptable.")
    args, _ = parser.parse_known_args() #second result is for unknown arguments
    return args

def main(**kwargs):
    args = getArgs()
    defaults = getStaticSettings()

    VARS = vars(args)
    given_args = {k: v for (k, v) in VARS.items() if v is not None}
    for key in kwargs: 
        given_args[key] = kwargs[key] #assign all keyword arguments to proper place in var dictionary

    print(given_args)

    if len(given_args) > 2 and args.file is not None:
        print("Please only provide one input molecule type.")
        quit()
    elif len(given_args) == 0:
        print("No arguments found. Please use `cutomPol -h` for options.")
        quit()

    if given_args["file"] is None:
        print("please provide a filename with the -f flag.")
        quit()

    if given_args["smiles"] is not None:
        smi = given_args["smiles"]
    elif given_args["smarts"] is not None:
        mol = Chem.MolFromSmarts(given_args["smarts"])
        smi = Chem.MolToSmiles(mol)
    elif given_args["inchi"] is not None:
        mol = Chem.MolFromInchi(given_args["inchi"])
        smi = Chem.MolToSmiles(mol)
    else:
        print("I am confused. Input format not recognized.")
        quit()

    POL_LIST = make_One_or_More_Polymers("", -1, smi, "", verbosity=True, custom=True)
    POL = POL_LIST[0]

    write_pol(given_args["file"], pol_list=POL.pol_list)

if __name__ == "__main__":
    main()
