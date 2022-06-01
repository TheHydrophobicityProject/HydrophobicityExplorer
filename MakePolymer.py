from rdkit import Chem
from rdkit.Chem import AllChem,Draw,Descriptors,rdFreeSASA
import argparse

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
    parser.add_argument("-s","--super_monomer", type=str, nargs='*',
                        help="a series of space-separated monomer SMILES arranged in their repeating sequence. You can add an int preceeding any monomer to represent multiple copies of that monomer. e.g. 2 A B means AAB is the repeating super-monomer. Use quotes surrounding SMILES with problematic characters like = or ()")
    parser.add_argument("-d","--draw", type=str, help="Filename for polymer image.")
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
        repeat_unit=monomer_dict[m]

    polymer_SMILES=smiles_inators[0]+n*repeat_unit+smiles_inators[1]
    
    #print(repeat_unit)
    #print(smiles_inators)
    print(polymer_SMILES)

    return polymer_SMILES

def bestConformer(pol_h,numConfs,seed,threads): #currently unused. See notes or question in optPol().
    ids = AllChem.EmbedMultipleConfs(pol_h, numConfs=numConfs, randomSeed=seed, useExpTorsionAnglePrefs=True, numThreads=threads)
    best=[]
    for id in ids:
        prop = AllChem.MMFFGetMoleculeProperties(pol_h)
        ff = AllChem.MMFFGetMoleculeForceField(pol_h, prop, confId=id)
        ff.Minimize()
        en = float(ff.CalcEnergy())
        econf = (en, id)
        best.append(econf)
    best.sort()
    if len(best==0):
        raise Exception("Error: No valid conformations found for the molecule. Try increasing the number of conformations.")
    # The best conformer is the first tuple and the ID of that conformer is the second value in the tuple
    best_id = int(best[0][1])
    # Return the best ID
    return best_id
    
def optPol(smiles):
    #make Mol object:
    pol=Chem.MolFromSmiles(smiles)
    #check mol
    Chem.SanitizeMol(pol)
    #opt steps
    pol_h=Chem.AddHs(pol)
    #pol_h=bestConformer(pol_h,100,42,2) #add support for changing these with cli arguments
    AllChem.EmbedMolecule(pol_h,useRandomCoords=True)
    AllChem.MMFFOptimizeMolecule(pol_h, maxIters=100) #does this repeated optimization obviate the need for the bestConformer function (copied from polymer LogP v4_4_4_alla*) 
    #maybe this number of itterations should be specified with cli arguments (give option).
    return pol_h, pol

def drawPol(pol,drawName):
    #will allow specified names later
    Chem.Draw.MolToFile(pol,drawName)

def LogP_Sasa(pol_h):#,best_conf_id):
    # Now calculate LogP and SASA
    # Calculate SASA based on the best conformer
    # classifyAtoms CRASHED when I tried it with , confIdx=best_conf_id
    # but someone needs to go back and make sure it's actually OK to use it without
    # and that that won't cause problems!
    radii = Chem.rdFreeSASA.classifyAtoms(pol_h)
    #sasa = Chem.rdFreeSASA.CalcSASA(pol_h, radii, confIdx=best_conf_id)
    sasa = Chem.rdFreeSASA.CalcSASA(pol_h, radii)
    # LogP does NOT have an option to feed in a conformer so just calculate it for the overall molecule
    logP = Chem.Descriptors.MolLogP(pol_h)
    # Now return LogP and SASA
    return logP, sasa

def main():
    args=getArgs()
    #two cases: one monomer or supermonomer.
    #if both are specified something is wrong.
    if args.single_monomer is not None and args.super_monomer is not None:
        raise argparse.ArgumentError("Cannot specify both single and super monomers")
    
    #This gives a list of components of a super-monomer or just the string used for single monomer in dict
    repeat_unit=list(filter(None,[args.single_monomer,args.super_monomer]))[0]

    polSMILES=createPolymerSMILES(args.initiator,args.n,repeat_unit,args.terminator)

    pol_h,pol=optPol(polSMILES)

    if args.draw is not None:
        drawName=args.draw.split(".")[0]+".png"
        drawPol(pol,drawName)

    logP,sasa=LogP_Sasa(pol_h)

    print(logP,sasa)
    
    return pol_h, pol, polSMILES, logP, sasa

main()
