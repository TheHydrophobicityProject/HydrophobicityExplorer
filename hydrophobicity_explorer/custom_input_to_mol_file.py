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
    args = parser.parse_args()
    return args

def main():
    args = getArgs()
    defaults = getStaticSettings()

    VARS = vars(args)
    given_args = {k: v for (k, v) in VARS.items() if v is not None}

    if len(given_args) > 2 and args.file is not None:
        print("Please only provide one input molecule type.")
        quit()

    if args.file is None:
        print("please provide a filename with the -f flag.")
        quit()

    if args.smiles is not None:
        smi = args.smiles
    elif args.smarts is not None:
        mol = Chem.MolFromSmarts(args.smarts)
        smi = Chem.MolToSmiles(mol)
    elif args.inchi is not None:
        mol = Chem.MolFromInchi(args.inchi)
        smi = Chem.MolToSmiles(mol)
    else:
        print("I am confused. Input format not recognized.")
        quit()

    POL_LIST = make_One_or_More_Polymers("", -1, smi, "", verbosity=True, custom=True)
    POL = POL_LIST[0]

    write_pol(args.file, pol_list=POL.pol_list)

if __name__ == "__main__":
    main()
