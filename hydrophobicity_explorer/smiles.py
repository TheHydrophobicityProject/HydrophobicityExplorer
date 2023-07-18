import argparse, pandas, os

monomer_dict = {

    #amino acids (use OH for termination group)
    'L-ala': 'N[C@@H](C)C(=O)',
    'beta-ala': 'NCCC(=O)',
    'L-cystiene': 'N[C@@H](CS)C(=O)',
    'glycine': 'NCC(=O)',
    'L-isoleucine': 'N[C@@H](C(C)CC)C(=O)',
    'L-leucine': 'N[C@@H](CC(C)C)C(=O)',
    'L-phenylala': 'N[C@@H](CC1=CC=CC=C1)C(=O)',
    'L-phenylgly': 'N[C@@H](C1=CC=CC=C1)C(=O)',
    'L-valine': 'N[C@@H](C(C)C)C(=O)',

    #monomers for polyamides (use H for termination group)
    'Nylon6': 'CCCCCC(=O)N',

    #monomers and comonomers for polycarbonates
    'BPA_Carbonate': 'c1ccc(cc1)C(C)(C)c1ccc(cc1)OC(=O)O',

    #monomers and comonomers for polyesters
    'Butylene_Adip': 'C(=O)CCCCC(=O)OCCCCO',
    'Butylene_Succinate': 'C(=O)CCC(=O)OCCCCO',
    'Butylene_Terephthalate': 'CCCCOC(=O)c1ccc(cc1)C(=O)O',
    'Ethylene_Terephthalate': 'CCOC(=O)c1ccc(cc1)C(=O)O',
    '3HBV': 'C(CC)CC(=O)O',
    '3HB': 'C(C)CC(=O)O',
    'Lactic_acid': 'C(C)C(=O)O',

    #monomers for polyethers
    'Ethylene_oxide': 'CCO',
    'EO': 'CCO',
    'Propylene_oxide': 'CC(C)O',

    #Vinyl monomers written to depict primary addition of alkene (i.e. substituent is on second carbon of alkene)
    'Acrylamide': 'CC(C(=O)N)',
    'Butylacrylate': 'CC(C(=O)OCCCC)',
    'Dimethylacrylamide': 'CC(C(=O)N(C)(C))',
    'DMA': 'CC(C(=O)N(C)(C))',
    'Dimethylaminoethylacrylate': 'CC(C(=O)OCC(N(C)(C)))',
    'Ethylene': 'CC',
    'Ethylacrylate': 'CC(C(=O)OCC)',
    'Ethylhexylacrylate': 'CC(C(=O)OCC(CC)CCCC)',
    '1-Hexene': 'CC(CCCC)',
    'Hydroxybutylacrylate': 'CC(C(=O)OCCCCO)',
    'Hydroxyethylacrylate': 'CC(C(=O)OCCO)',
    'Hydroxyethylmethacrylate': 'CC(C(=O)OCCO)(C)',
    'Methoxyethylacrylate': 'CC(C(=O)OCCOC)',
    'Methylacrylate': 'CC(C(=O)OC)C',
    'Methylmethacrylate': 'CC(C(=O)OC)(C)',
    'MMA': 'CC(C(=O)OC)(C)',
    'NVP': 'CC(N1CCCC1=O)',
    'Propylene': 'CC(C)',
    'Styrene': 'CC(c1ccccc1)',
    'Vinylalcohol': 'CC(O)',
    'Vinylchloride': 'CC(Cl)'
}

#initiator dictionary
init_dict = {
    'Benzyl': '*Cc1ccccc1',
    'Benzyl_alcohol': '*Cc1ccccc1O',
    'Benzoyl': 'c1ccc(cc1)C(=O)O*',
    'Butyl': 'CCCC',
    'Hydroxyl': 'O',
    'Hydroxy': 'O',
    'Hydrogen': '',
    'n-Butyl': 'CCCC',
    'Methoxy': 'CO*',
    'Ethoxy': 'CCO*',
    'Methyl': 'C',
    'Vinyl': 'C=C'
}

template = {
    "end_groups": {
        "KEY0": "SYMMETRICSMILES",
        "KEY1": "*ASYMMETRIC SMILES"
    },
    "monomers": {
        "KEY": "HEAD_SMILES_TAIL"
    }
}


def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--end_group', default = False, action = "store_true", help = "Prints all end group key-value pairs.")
    parser.add_argument('-m', '--monomer', default = False, action = "store_true", help = "Prints all monomer key-value pairs.")
    parser.add_argument('-w', '--write', default = False, action = "store_true", help = "Writes a json file that can be filled with user-specified smiles.")
    args = parser.parse_args()
    return args


def _addUserSmiles(user_dict, endgroup_dict=None, mnmr_dict=None):
    if endgroup_dict is not None:
        user_eg_dict = user_dict["end_groups"]
        for key in user_eg_dict:
            if "KEY" not in key:
                endgroup_dict[key] = user_eg_dict[key]
        return endgroup_dict
    if mnmr_dict is not None:
        user_monomer_dict = user_dict["monomers"]
        for key in user_monomer_dict:
            if "KEY" not in key:
                mnmr_dict[key] = user_monomer_dict[key]
        return mnmr_dict


def showDict(raw_dict):
    dictionary = {"KEYS": [], "SMILES": []}
    for key in raw_dict:
        dictionary['KEYS'].append(key)
        dictionary['SMILES'].append(raw_dict[key])
    df = pandas.DataFrame(dictionary)
    print(df)
    # print(df.to_string(index=False))
    print("\n")


def checkAndMergeSMILESDicts(egs, mnmrs):
    if os.path.exists("smiles.json"):
        print("User-made key-value pairs will be shown as well.\n")
        from hydrophobicity_explorer.settings import readJson
        user_dict = readJson("smiles.json")
        init_dict = egs
        monomer_dict = mnmrs
        init_dict = _addUserSmiles(user_dict, endgroup_dict=egs)
        monomer_dict = _addUserSmiles(user_dict, mnmr_dict=mnmrs)
        return init_dict, monomer_dict
    else:
        return egs, mnmrs


def main():
    args = getArgs()
    if args.end_group or args.monomer:
        egs, mnmrs = checkAndMergeSMILESDicts(init_dict, monomer_dict)
    if args.end_group:
        showDict(egs)
    if args.monomer:
        showDict(mnmrs)
    if args.write:
        from hydrophobicity_explorer.settings import writeJson
        name = "smiles.json"
        if os.path.exists(name):
            inp = input(f"{name} exists. Should it be overwritten? [Y/n]: ")
            if inp.lower() == "y" or inp == "":
                writeJson(template, name)
                print("created smiles.json")
            else:
                print("Please rename the existing notebook so it is not overwritten.")
        else:
            writeJson(template, name)
            print("created smiles.json")


if __name__ == "__main__":
    main()
