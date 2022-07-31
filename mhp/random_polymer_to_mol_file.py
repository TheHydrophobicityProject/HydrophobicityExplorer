import rdkit, argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from random import choices, shuffle
from smiles import monomer_dict
from MakePolymer import validate_end_group, inator_smi_lookup, add_inator_smiles, optPol, getStaticSettings

def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type=int, help="Total number of monomers desired in final product.")
    parser.add_argument("-i", "--initiator", type=str, default="Hydrogen", help="Initiator Key from initiator dict or SMILES. Defaults to Hydrogen.")
    parser.add_argument("-m", type=str, nargs='*',
            help="Space-separated list of monomer keys from smiles.py or a custom smiles string for a monomer unit. The head of the monomer must be to the left and the attachment points must be at the ends of the string. Integer coefficients can be used to influence composition of the polymer. The exact method used depends on the use of the -p flag.")
    parser.add_argument("-t", "--terminator", type=str, default="Hydrogen", help="Terminator key taken from initiator dict or SMILES. Defaults to Hydrogen.")

    parser.add_argument("-f", type=str, help="Filename prefix you would like to use. The number of monomers will be added automatically. Input will be forced to use .sdf format.")
    parser.add_argument("-p", "--protocol", type=str, default="ratio",
            help="The method of polymer body generation. The valid options are \"weight\" or \"ratio\" [default]. The coefficients that lead each monomer will be converted to a percentage of the total monomers and a polymer with a random ordering of monomers of that specific composition will be generated. If you specify \"weight\" the polymer composition will be random with weights favoring the monomers with higher coefficients. This could result in polymers with a different emperical formula than the coefficient would indicate.")
    args = parser.parse_args()
    return args

def prepFilename(filename, n):
    split = filename.split(".") #split off file extention in case provided.
    name = f"{split[0]}_{n}.sdf"
    return name

def getCoeffs(Coefs_and_monomers):
    #first, let us expand the provided coefs to actually give a weighted list.
    coeffs_list = []
    monomer_list = []
    #ommission of a coeficient implies 1 copy
    repeat_coef = 1
    for element in Coefs_and_monomers:
        try:
            repeat_coef = int(element) #is this a string of an integer?
        except:
            #if not string of integer, it should be considered a smiles (there will be griping from rdkit if not.)
            coeffs_list.append(repeat_coef)
            monomer_list.append(element)
            repeat_coef = 1 #reset coef.

    return coeffs_list, monomer_list

def mergeList(lst):
    smiles = "".join(lst)
    print(f'body smiles is:\n{smiles}')
    return smiles


def makePolymerBody_weighted(weighted_monomer_list, n):
    monomer_weights, expanded_list = getCoeffs(weighted_monomer_list)

    #now that we have an explicit list, we can pick at random from the list n times and make the polymer body.
    body_list = choices(expanded_list, weights=monomer_weights, k=n)

    #now merge list into smiles.
    smiles = mergeList(body_list)

    return smiles

def makePolymerBody_ratio(formula_list, n):
    coeffs, monomers = getCoeffs(formula_list)
    print(monomers)
    sum_coeffs = sum(coeffs)
    body_list = []
    for i, coeff in enumerate(coeffs):
        real_coeff = round(coeff / sum_coeffs * n) #we can't have float coeffs. Even though we are rounding everything should even out unless two monomers' coeffs initially end in .5
        body_list += [item for item in [monomers[i]] for j in range(real_coeff)] #we repeat the monomer an appropriate number of times

    #if n happens to be != to length of list, maybe allow user to modify by popping a random monomer from list or chosing a type of monomer to remove (and we remove a random one of that type.)
    #This isn't super important right now since they have a smiles they can also modify if they want exactly the right number.
    if n != len(body_list):
        print(f"Due to rounding the length of the polymer is {len(body_list)} and the n specified was {n}.")

    #Now we have a list of monomers in the correct relative ammounts.
    #we just need to randomize the order.
    shuffle(body_list)
    print(body_list)
    smiles = mergeList(body_list)
    return smiles

def main():
    args = getArgs()
    defaults = getStaticSettings()

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
    if args.protocol == "weight":
        polymer_body_smiles = makePolymerBody_weighted(deciphered_dict_keys, args.n)
    elif args.protocol == "ratio":
        polymer_body_smiles = makePolymerBody_ratio(deciphered_dict_keys, args.n)
    else:
        argparse.ArgumentError("Unknown protocol selected. Use only \"ratio\" or \"weight\". \"ratio\" is the default if the -p flag is not used.")

    #now we need to attatch the end groups:
    total_smiles = add_inator_smiles(polymer_body_smiles, init, term)
    print("Finished adding end groups. Beginning optimization.")

    pol, suppl = optPol(total_smiles, name=file_name, 
        nConfs=defaults["opt_numConfs"], threads=defaults["opt_numThreads"], iters=defaults["opt_maxIters"]) #this function will also save the file.

    print(f"Done. Saved to {file_name}")

if __name__ == "__main__":
    main()
