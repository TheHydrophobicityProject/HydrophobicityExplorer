from rdkit import Chem
from rdkit.Chem import AllChem,Draw,Descriptors,rdFreeSASA
import argparse,os

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
    parser.add_argument("-n", type=int, help="The number of monomer or super-monomer units.")
    parser.add_argument("-i", "--initiator", type=str, default="Hydrogen", help="Initiator Key from initiator dict or SMILES. Defaults to Hydrogen.")
    parser.add_argument("-t", "--terminator", type=str, default="Hydrogen", help="Terminator key taken from initiator dict or SMILES. Defaults to Hydrogen.")
    parser.add_argument("-m","--single_monomer", type=str, help="Monomer key from the included monomer dict. See the -s flag for specifying a monomer that is not included.")
    parser.add_argument("-s", "--super_monomer", type=str, nargs='*',
                        help="A series of space-separated monomer SMILES arranged in their repeating sequence. You can add an int preceeding any monomer to represent multiple copies of that monomer. e.g. 2 A B means AAB is the repeating super-monomer. Use quotes surrounding SMILES with problematic characters like = or ()")
    parser.add_argument("-d", "--draw", type=str, help="Filename for polymer image.")
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="Set increased verbosity. Will draw polymer to polymer.png unless alternate name set by -d option.")
    parser.add_argument("-c","--calculation", type=str, nargs='*', help="Type of calculation(s) to be performed input as a space-separated list. Options are LogP, SA (surface area) and RG (radius of gyration).")
    parser.add_argument("-f","--file", type=str, help="The name/path of the file you wish to save the mol to. Supported formats are .pdb, .xyz and .mol")
    parser.add_argument("-r", "--read", type=str, help="The name/path to file you wish to import. Supported formats are .pdb and .mol")
    args=parser.parse_args()
    return args

def createPolymerSMILES(i,n,m,t):
    given_inators = [i,t]
    #gets from dict if available. Otherwise assume SMILES and continue. There will eventually be an error if this isn't the case.
    smiles_inators=[init_dict[x] if x in init_dict else x for x in given_inators]

    if type(m)==list:
        repeat_unit=""
        repeat=1 #ommission of a coeficient implies 1 copy
        for element in m:
            try:
                repeat=int(element) #is this a string of an integer?
            except:
                repeat_unit+=repeat*element #if not, repeat the SMILES as many times as specified (or once if no coef. provided).
                repeat=1 #reset coef.
    else:
        repeat_unit=monomer_dict[m] #if not a list, look for the corresponding smiles in the dictionary, will throw error if not included.
    
    polymer_SMILES=smiles_inators[0]+n*repeat_unit+smiles_inators[1] #concatonate all parts as a big SMILES.
    return polymer_SMILES

def bestConformer(pol_h,numConfs,seed): #currently unused. See notes or question in optPol().
    ids = AllChem.EmbedMultipleConfs(pol_h, numConfs=numConfs, randomSeed=seed, useExpTorsionAnglePrefs=True, numThreads=0)
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
    Chem.Draw.MolToFile(pol,drawName)

def write_or_read_pol(pol_h_or_str,name):
    ext=name.split(".")[1]
    if type(pol_h_or_str)==str: #"Read" if in read mode.
        if os.path.exists(name):
            if ext == "pdb":
                pol_h = Chem.MolFromPDBFile(name)
            elif ext == "mol":
                pol_h = Chem.MolFromMolFile(name)
            else:
                print("unsuported extention:",ext,"Please use .pdb, or .mol") #.xyz cannot be read by rdkit.
                quit()
            return pol_h
        else:
            raise FileNotFoundError(name)
    else:
        exit_status=0
        if ext == "xyz":
            Chem.MolToXYZFile(pol_h_or_str,name)
        elif ext == "pdb":
            Chem.MolToPDBFile(pol_h_or_str,name)
        elif ext == "mol":
            Chem.MolToMolFile(pol_h_or_str,name)
        else:
            print("unsuported extention:",ext,"Please use .pdb, .xyz or .mol")
            exit_status=1
        return exit_status

def Sasa(pol_h):#,best_conf_id):
    # Now calculate LogP and SASA
    # Calculate SASA based on the best conformer
    # classifyAtoms CRASHED when I tried it with, confIdx=best_conf_id
    # but someone needs to go back and make sure it's actually OK to use it without
    # and that that won't cause problems!
    
    radii = Chem.rdFreeSASA.classifyAtoms(pol_h)
    #sasa = Chem.rdFreeSASA.CalcSASA(pol_h, radii, confIdx=best_conf_id)
    sasa = Chem.rdFreeSASA.CalcSASA(pol_h, radii)    
    return sasa

def LogP(pol_h):
    # LogP does NOT have an option to feed in a conformer so just calculate it for the overall molecule
    logP = Chem.Descriptors.MolLogP(pol_h)
    #I've seen that there is a way to include or exclude hydrogens from the LogP calculation.
    #Descriptors.MolLogP(z_no_Hs, includeHs=*bool*)). Do we want to include this?
    return logP

def RadGyration(pol_h):
    RG=Chem.rdMolDescriptors.CalcRadiusOfGyration(pol_h)
    #Chem.Descriptors3D.RadiusOfGyration(pol_h)
    #both seem to give identical results based on "SMILES to Rg.ipynb"
    return RG

def doCalcs(pol_h,calcs):
    #Calcs are only done if requested.
    #Not a fan of nested if statements. Open to suggestions on improvments.
    data={}
    for calc in calcs:
        if calc == "SA":
            sasa=Sasa(pol_h)
            data["SA"]=sasa
        elif calc == "LogP":
            logP=LogP(pol_h)
            data["LogP"]=logP
        elif calc == "RG":
            rg=RadGyration(pol_h)
            data["RG"]=rg
        else:
            print("unrecognized calculation:",calc+". Use SA, LogP or RG")
    return data

def main():
    args=getArgs()
    if args.read is None: #then get polymer parameters from CLI arguments.
        #Two cases: one monomer or supermonomer.
        #If both are specified something is wrong.
        if args.single_monomer is not None and args.super_monomer is not None:
            raise argparse.ArgumentError("Cannot specify both single and super monomers")
        
        #This gives a list of components of a super-monomer or just the string used for single monomer in dict
        repeat_unit=list(filter(None,[args.single_monomer,args.super_monomer]))[0]

        polSMILES=createPolymerSMILES(args.initiator,args.n,repeat_unit,args.terminator)

        if args.verbose: #extra info if requested
            print("Polymer interpreted as:",args.initiator,str(args.n),"*",str(repeat_unit),args.terminator)
            print("This gives the following SMILES:",polSMILES)

        pol_h,pol=optPol(polSMILES) #both are mol objects
    else: #get mol from file
        pol_h=write_or_read_pol("Read",args.read)

    if args.verbose: #more extra information
        print("requested calculations are",args.calculation)

    #saving the polymer to a file.
    if args.file is not None: #technically nothing wrong with using this as a roundabout way of converting between filetypes
        if args.verbose:
            print("attempting to save molecule to",args.file)
        stat=write_or_read_pol(pol_h,args.file)
        if args.verbose and stat == 0:
            print("success.")

    #drawing an picture of the polymer.
    #drawings with optimized geoms are ugly in 2D. Fix so it "flattens out" so we can get drawing of unknown file
    if args.draw is not None and args.read is None:
        drawName=args.draw.split(".")[0]+".png"
        drawPol(pol,drawName)
    else:
        if args.verbose and args.read is None:
            #if we can flatten out the mol change this trigger too.
            #produce image if increased verbocity is requested even if no name is set.
            drawPol(pol,"polymer.png")

    #doing only the specified calculations.
    if args.calculation is not None:
        data=doCalcs(pol_h,args.calculation)
        print(data)

main()
