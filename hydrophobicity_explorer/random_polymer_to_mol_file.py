from random import choices, shuffle
import rdkit, argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from functools import reduce
from math import gcd


def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type=int, help="Total number of monomers desired in final product.")
    parser.add_argument("-i", "--initiator", type=str, default="Hydrogen", help="Initiator Key from initiator dict or SMILES. Defaults to Hydrogen.")
    parser.add_argument("-m", type=str, nargs='*',
            help="Space-separated list of monomer keys from smiles.py or a custom smiles string for a monomer unit. The head of the monomer must be to the left and the attachment points must be at the ends of the string. Integer coefficients can be used to influence composition of the polymer. The exact method used depends on the use of the -p flag.")
    parser.add_argument("-t", "--terminator", type=str, default="Hydrogen", help="Terminator key taken from initiator dict or SMILES. Defaults to Hydrogen.")

    parser.add_argument("-f", type=str, help="Filename prefix you would like to use. The number of monomers will be added automatically. Formats of PDB, XYZ and Mol are acceptable.")
    parser.add_argument("-p", "--protocol", type=str, default="ratio",
            help="The method of polymer body generation. The valid options are \"weight\" or \"ratio\" [default]. The coefficients that lead each monomer will be converted to a percentage of the total monomers and a polymer with a random ordering of monomers of that specific composition will be generated. If you specify \"weight\" the polymer composition will be random with weights favoring the monomers with higher coefficients. This could result in polymers with a different emperical formula than the coefficient would indicate.")
    parser.add_argument("-a", "--array", default = False, action = "store_true", help="generate an array of polymers from n=1 to the n specified with the -n flag.")
    args = parser.parse_args()
    return args


def getCoeffs(Coefs_and_monomers):
    #first, let us expand the provided coefs to actually give a weighted list.
    coeffs_list = []
    monomer_list = []
    #ommission of a coeficient implies 1 copy
    repeat_coef = 1
    for element in Coefs_and_monomers:
        try:
            repeat_coef = int(element)
        except:
            #if not string of integer, it should be considered a smiles (there will be an error from rdkit if not.)
            coeffs_list.append(repeat_coef)
            monomer_list.append(element)
            repeat_coef = 1

    return coeffs_list, monomer_list


def makePolymerBody_weighted(weighted_monomer_list, n):
    monomer_weights, expanded_list = getCoeffs(weighted_monomer_list)
    #now that we have an explicit list, we can pick at random from the list n times and make the polymer body.
    body_list = choices(expanded_list, weights=monomer_weights, k=n)
    #now merge list into smiles.
    smiles = "".join(body_list)

    return smiles


def makePolymerBody_ratio(formula_list, n, verbo=False, mpn=1):
    coeffs, monomers = getCoeffs(formula_list)
    sum_coeffs = sum(coeffs)
    body_list = []
    real_coeffs = []
    roundup = True
    for i, coeff in enumerate(coeffs):
        unrounded_coeff = (coeff / sum_coeffs) * n * mpn
        #We need to make the case where 2 monomers have 0.5 split one to go up and the other down. However, "Integers" should not be changed.
        if unrounded_coeff % 0.5 == 0 and unrounded_coeff % 1 != 0:
            if roundup:  #the first one we see is rounded up
                unrounded_coeff += 0.1
                roundup = False
            else:  #the second one is rounded down and the flag is reset
                roundup = True
                unrounded_coeff -= 0.1
        #now we round and hopefully the length should always be correct
        real_coeff = round(unrounded_coeff) #we can't have float coeffs. Even though we are rounding everything should even out unless two monomers' coeffs initially end in .5
        body_list += [item for item in [monomers[i]] for _ in range(real_coeff)] #we repeat the monomer an appropriate number of times
        real_coeffs.append(real_coeff)

    #if n happens to be != to length of list, maybe allow user to modify by popping a random monomer from list or chosing a type of monomer to remove (and we remove a random one of that type.)
    #This isn't super important right now since they have a smiles they can also modify if they want exactly the right number.
    if n * mpn != len(body_list):
        print(f"WARNING: Due to rounding the length of the polymer is {len(body_list)} and the n specified was {n}.")
    if len(body_list) == 0:
        return None, "0:0"

    den = reduce(gcd, real_coeffs)
    ratio = ":".join(str(int(i / den)) for i in real_coeffs)

    #Now we have a list of monomers in the correct relative ammounts.
    #we just need to randomize the order.
    shuffle(body_list)
    smiles = "".join(body_list)
    if verbo:
        print(f"{n = }")
        print(f"Ratio of monomers used is {ratio}")
    return smiles, ratio


def main():
    #done in every instance.
    args = getArgs()
    defaults = getStaticSettings()
    #first get inators sorted out:
    #get their smiles
    init, term = inator_smi_lookup(args.initiator, args.terminator)
    #make sure they are formatted the right way.
    init = validate_end_group(init, Init=True)
    term = validate_end_group(term, Term=True)
    #replace any dict keys with corresponding smiles.
    deciphered_dict_keys = parse_smiles_dict_keys(args.m, monomer_dict)

    #do these steps multiple times if array of files is requested.
    if not args.array:
        n_iter = [args.n]
    else:
        n_iter = range(1, args.n + 1)

    for n in n_iter:
        #get proper file name.
        filename = args.f
        split = filename.split(".")
        file_name = f"{split[0]}_{n}.{split[1]}"
        #now we need to generate polymer body
        if args.protocol == "weight":
            polymer_body_smiles = makePolymerBody_weighted(
                deciphered_dict_keys, n)
        elif args.protocol == "ratio":
            polymer_body_smiles, ratio = makePolymerBody_ratio(
                deciphered_dict_keys, n)
        else:
            argparse.ArgumentError("Unknown protocol selected. Use only \"ratio\" or \"weight\". \"ratio\" is the default if the -p flag is not used.")

        #now we need to attatch the end groups:
        total_smiles = add_inator_smiles(polymer_body_smiles, init, term)
        print("Finished adding end groups. Beginning optimization.")
        if args.array:
            POL = Polymer(n = args.n, smiles=total_smiles)
        else:
            POL = Polymer(n = n, smiles=total_smiles)

        POL.suppl = optPol(total_smiles,
                            nConfs=defaults["opt_numConfs"],
                            threads=defaults["opt_numThreads"],
                            iters=defaults["opt_maxIters"]
                            )  #this function will also save the file.

        write_pol(name=file_name, verbosity=True, suppl=POL.suppl)
        print(f"Done. Saved to {file_name}")
        
        if args.protocol == "ratio":
            print(f"Ratio of monomers used is {ratio}")


if __name__ == "__main__":
    from hydrophobicity_explorer.smiles import monomer_dict
    from hydrophobicity_explorer.MakePolymer import write_pol, validate_end_group, inator_smi_lookup, add_inator_smiles, optPol, getStaticSettings, parse_smiles_dict_keys, Polymer
    main()
