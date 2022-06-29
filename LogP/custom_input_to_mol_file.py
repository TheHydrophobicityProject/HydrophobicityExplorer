import rdkit, argparse
from rdkit import Chem

def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--smiles", type=str, help="Complete smiles string to be converted.")
    parser.add_argument("-m", "--smarts", type=str, help="Complete smarts string to be converted.")
    parser.add_argument("-i", "--inchi", type=str, help="Complete inchi string to be converted.")
    parser.add_argument("-f", "--file", type=str, help="filename you would like to use. The format NAME_numberOfMonomers is best.")
    args = parser.parse_args()
    return args

def main():
    args = getArgs()
    VARS=vars(args)
    given_args = { k:v for (k, v) in VARS.items() if v is not None }

    if len(given_args) > 2 and args.file is not None:
        print("Please only provide one argument.")
        quit()
    if args.file is None:
        print("please provide a filename with the -f flag.")
        quit()

    if args.smiles is not None:
        mol=Chem.MolFromSmiles(args.smiles)
    elif args.smarts is not None:
        mol=Chem.MolFromSmarts(args.smarts)
    elif args.inchi is not None:
        mol = Chem.MolFromInchi(args.inchi)
    else:
        print("I am confused")
        quit()
    
    Chem.MolToMolFile(mol, args.file)

if __name__ == "__main__":
    main()

#need to add checks for filename. Appropriate denotion of number of monomers is required.
#need to check file extention too
#add to README
#should optimize mol to get coordinates too.