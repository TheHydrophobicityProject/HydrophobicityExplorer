from rdkit import Chem
from rdkit.Chem import AllChem
import argparse
from threading import local

###Built-in dict. Can be divorced from this file later

monomer_dict={

    #monomers for polyamides
    'Nylon6': 'CCCCCC(=O)N',

    #monomers and comonomers for polycarbonates
    'BPA_Carbonate':  'c1ccc(cc1)C(C)(C)c1ccc(cc1)OC(=O)O',

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
    'Propylene_oxide': 'CC(C)O',

    #Vinyl monomers written to depict primary addition of alkene (i.e. substituent is on second carbon of alkene)
    'Butylacrylate': 'CC(C(=O)OCCCC)',
    'Dimethylacrylamide': 'CC(C(=O)N(C)(C))',
    'Ethylene': 'CC',
    'Hydroxyethylmethacrylate': 'CC(C(=O)OCCO)(C)',
    'Methylacrylate': 'CC(C(=O)OC)C',
    'Methylmethacrylate': 'CC(C(=O)OC)(C)',
    'Propylene': 'CC(C)',
    'Styrene': 'CC(c1ccccc1)',
    'Vinylalcohol': 'CC(O)'
}

#initiator dictionary (i.e. endgroup of polymer chain)
init_dict={'Benzyl': 'c1ccccc1CO', 
           'Butyl': 'CCCC',
           'Hydroxyl': 'O', 
           'Hydrogen': '',
           'Methoxy': 'CO', 
           'Methyl': 'C', 
           'Vinyl': 'C=C'
           }

def getArgs():
    parser=argparse.ArgumentParser()
    parser.add_argument("-n", type=int, help="number of monomer or super-monomer units.")
    parser.add_argument("-i", "--initiator", type=str, default="Hydrogen", help="initiator Key from initiator dict. Defaults to Hydrogen.")
    parser.add_argument("-t", "--terminator", type=str, default="Hydrogen", help="terminator key taken from initiator dict. Defaults to Hydrogen.")
    parser.add_argument("-m","--single_monomer", type=str, help="monomer key from any of the included monomer dicts. Use -s instead to specify a monomer that is not included.")
    #parser.add_argument("-M","--multiple_monomers", type=bool, help="Use multiple monomers? Any argument means True.", default=False)
    parser.add_argument("-s","--super_monomer", type=str, nargs='*',
                        help="a series of space-separated monomer SMILES arranged in their repeating sequence. You can add an int preceeding any monomer to represent multiple copies of that monomer. e.g. 2 A B means AAB is the repeating super-monomer. Use quotes surrounding SMILES with problematic characters like = or ()")
    args=parser.parse_args()

    return args

def createPolymerSMILES(i,n,m,t):
    print(i,n,m,t)
    
    given_inators = [i,t]
    #gets from dict if available. Otherwise assume SMILES and continue.
    smiles_inators=[init_dict[x] if x in init_dict else x for x in given_inators]

    if type(m)==list:
        repeat_unit=""
        repeat=1
        for element in m:
            try:
                repeat=int(element)
            except:
                repeat_unit+=repeat*element
                repeat=1
    else:
        repeat_unit=m

    polymer_SMILES=smiles_inators[0]+n*repeat_unit+smiles_inators[1]
    
    # print(repeat_unit)
    # print(smiles_inators)
    # print(polymer_SMILES)

    return polymer_SMILES

def main():
    args=getArgs()
    #two cases: one monomer or supermonomer.
    #if both are specified something is wrong.
    if args.single_monomer is not None and args.super_monomer is not None:
        raise argparse.ArgumentError("Cannot specify both single and super monomers")
    
    #This gives a list of components of a super-monomer or just the string used for single monomer in dict
    repeat_unit=list(filter(None,[args.single_monomer,args.super_monomer]))[0]

    polSMILES=createPolymerSMILES(args.initiator,args.n,repeat_unit,args.terminator)

    #make Mol object:
    pol=Chem.MolFromSmiles(polSMILES)
    #opt steps
    pol_h=Chem.AddHs(pol)
    AllChem.EmbedMolecule(pol_h)
    AllChem.MMFFOptimizeMolecule(pol_h)

    return pol_h, pol, polSMILES

main()
