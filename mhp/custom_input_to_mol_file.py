import rdkit, argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from MakePolymer import optPol

def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--smiles", type=str, help="Complete smiles string to be converted.")
    parser.add_argument("-m", "--smarts", type=str, help="Complete smarts string to be converted.")
    parser.add_argument("-i", "--inchi", type=str, help="Complete inchi string to be converted.")
    parser.add_argument("-f", "--file", type=str, help="filename you would like to use. The format NAME_numberOfMonomers.mol is required.")
    args = parser.parse_args()
    return args

def checkFilename(filename):
    split = filename.split(".")

    if split[-1] != "mol":
        print("Bad file extention. Please use Name_n.mol")
        quit()

    n = split[0].split("_")[-1]
    try:
        int(n)
    except:
        print("Please add an N to your filename. Name_n.mol")
        quit()

def main():
    args = getArgs()
    VARS = vars(args)
    given_args = { k:v for (k, v) in VARS.items() if v is not None }

    if len(given_args) > 2 and args.file is not None:
        print("Please only provide one input molecule type.")
        quit()

    if args.file is None:
        print("please provide a filename with the -f flag.")
        quit()
    else:
        checkFilename(args.file)

    if args.smiles is not None:
        mol = Chem.MolFromSmiles(args.smiles)
    elif args.smarts is not None:
        mol = Chem.MolFromSmarts(args.smarts)
    elif args.inchi is not None:
        mol = Chem.MolFromInchi(args.inchi)
    else:
        print("I am confused")
        quit()

    #doing this extra step so opt is consistent with option used in primary script. Only need to change one set of parameters while finding best options.
    mol_h, mol = optPol(Chem.MolToSmiles(mol)) 
    
    Chem.MolToMolFile(mol_h, args.file)

if __name__ == "__main__":
    main()
