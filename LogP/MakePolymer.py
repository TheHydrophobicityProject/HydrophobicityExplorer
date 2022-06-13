from functools import cache
from rdkit import Chem
from rdkit.Chem import AllChem,Draw,Descriptors,rdFreeSASA
import argparse,os,csv,json
import matplotlib.pyplot as plt
from smiles import monomer_dict, init_dict

def getJsonArgs(jsonFile, dict):
    with open(jsonFile, 'r') as J:
        runs_dict = json.load(J)
        for run in runs_dict["runs"]: #so few items the nested for loops shouldn't be a big deal
            run_keys = run.keys()
            for dict_key in dict.keys():
                if dict_key not in run_keys:
                    run[dict_key] = dict[dict_key]
    run_list = runs_dict["runs"]
    return run_list

def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type = int, help = "The number of monomer or super-monomer units.")
    parser.add_argument("-i", "--initiator", type = str, default = "Hydrogen", help = "Initiator Key from initiator dict or SMILES. Defaults to Hydrogen.")
    parser.add_argument("-t", "--terminator", type = str, default = "Hydrogen", help = "Terminator key taken from initiator dict or SMILES. Defaults to Hydrogen.")
    parser.add_argument("-m","--single_monomer", type = str, help = "Monomer key from the included monomer dict. See the -s flag for specifying a monomer that is not included.")
    parser.add_argument("-s", "--super_monomer", type = str, nargs = '*',
                        help = "A series of space-separated monomer SMILES arranged in their repeating sequence. You can add an int preceeding any monomer to represent multiple copies of that monomer. e.g. 2 A B means AAB is the repeating super-monomer. Use quotes surrounding SMILES with problematic characters like = or ()")
    parser.add_argument("-d", "--draw", type = str, help = "Filename for polymer image.")
    parser.add_argument("-v", "--verbose", default = False, action = "store_true", help = "Set increased verbosity. Will draw polymer to polymer.png unless alternate name set by -d option.")
    parser.add_argument("-c","--calculation", type = str, nargs = '*', 
                        help = "Type of calculation(s) to be performed input as a space-separated list. Options are LogP, SA (surface area), MV (Molecular Volume), MHP (Mathers Hydrophobicity Parameter (LogP/SA; each of which will also be reported. Use XMHP to exclude those plots)) and RG (radius of gyration).")
    parser.add_argument("-f","--file", type = str, help = "The name/path of the file you wish to save the mol to. Supported formats are .pdb, .xyz and .mol")
    parser.add_argument("-r", "--read", type = str, help = "The name/path to file you wish to import. Supported formats are .pdb and .mol")
    parser.add_argument("-p", "--plot", default = False, action = "store_true", 
                        help = "Include this option to generate a plot of whatever calculations are specified with -c on polymers from 1 to the n specified with the -n flag. This means the molecule cannot be read from a file with the -r flag. If used with the -f flag multiple files will be saved with names based off the one provided.")
    parser.add_argument("-e", "--export", type = str, help = "Include this option to export a .csv file of all data calculations. Specify the name here.")
    parser.add_argument("-j", "--json", type = str, help = "The path to a compatible .json file with any of the above arguments.")
    args = parser.parse_args()
    #get additional arguments from json file if provided or by default if no args provided.
    vardict = vars(args)
    if args.json is not None:
        run_list = getJsonArgs(args.json, vardict)
    else:
        run_list = [vardict]
    return run_list

def getRepeatUnit(single, super):
    #Two cases: one monomer or supermonomer.
    #If both are specified something is wrong.
    if single is not None and super is not None:
        raise TypeError("Cannot specify both single and super monomers")
    #This gives a list of components of a super-monomer or just the string used for single monomer in dict
    repeat_unit = list(filter(None, [single, super]))[0]
    return repeat_unit

@cache #avoid multiple lookups if multiple runs with same inputs
def monomer_smi_lookup(m):
    repeat_unit = monomer_dict[m]
    return repeat_unit

@cache #avoid multiple lookups if multiple runs with same inputs
def inator_smi_lookup(i,t):
    given_inators = [i,t]
    #gets from dict if available. Otherwise assume SMILES and continue. There will eventually be an error if this isn't the case.
    smiles_inators = [init_dict[x] if x in init_dict else x for x in given_inators]
    return smiles_inators

def get_building_blocks(i,t,m):
    smiles_inators = inator_smi_lookup(i,t)
    if type(m) == list:
        repeat_unit = ""
        repeat = 1 #ommission of a coeficient implies 1 copy
        for element in m:
            try:
                repeat = int(element) #is this a string of an integer?
            except:
                repeat_unit += repeat * element #if not, repeat the SMILES as many times as specified (or once if no coef. provided).
                repeat = 1 #reset coef.
    else:
        repeat_unit = monomer_smi_lookup(m) #if not a list, look for the corresponding smiles in the dictionary, will throw error if not included.
    
    init = smiles_inators[0]
    term = smiles_inators[1]
    
    return init, term, repeat_unit

def createPolymerSMILES(i,n,m,t):
    init, term, repeat_unit = get_building_blocks(i,t,m)
    polymer_SMILES = init + n * repeat_unit + term #concatonate all parts as a big SMILES.
    return polymer_SMILES

def bestConformer(pol_h,numConfs,seed): #currently unused. See notes or question in optPol().
    ids = AllChem.EmbedMultipleConfs(pol_h, numConfs=numConfs, randomSeed=seed, useExpTorsionAnglePrefs=True, numThreads=0)
    best = []
    for id in ids:
        prop = AllChem.MMFFGetMoleculeProperties(pol_h)
        ff = AllChem.MMFFGetMoleculeForceField(pol_h, prop, confId=id)
        ff.Minimize()
        en = float(ff.CalcEnergy())
        econf = (en, id)
        best.append(econf)
    best.sort()
    if len(best == 0):
        raise Exception("Error: No valid conformations found for the molecule. Try increasing the number of conformations.")
    # The best conformer is the first tuple and the ID of that conformer is the second value in the tuple
    best_id = int(best[0][1])
    # Return the best ID
    return best_id
    
def optPol(smiles):
    #make Mol object:
    pol = Chem.MolFromSmiles(smiles)
    #check mol
    Chem.SanitizeMol(pol)
    #opt steps
    pol_h = Chem.AddHs(pol)
    #pol_h = bestConformer(pol_h,100,42,2) #add support for changing these with cli arguments
    AllChem.EmbedMolecule(pol_h, useRandomCoords = True)
    AllChem.MMFFOptimizeMolecule(pol_h, maxIters = 100) #does this repeated optimization obviate the need for the bestConformer function (copied from polymer LogP v4_4_4_alla*) 
    #maybe this number of itterations should be specified with cli arguments (give option).
    return pol_h, pol

def make_One_or_More_Polymers(i,n,r,t, *, verbosity=False, plot=False):
    POL_LIST = []
    SMI_LIST = []
    Unopt_pols = []
    if plot:
        N_array = range(1, n+1)
        for j in N_array:
            smi = createPolymerSMILES(i,j,r,t)
            if verbosity:
                print(f"Done generating SMILES with n = {j} now: {smi}")
                print("Converting to mol now.")
            pol_h,pol = optPol(smi)
            POL_LIST.append(pol_h)
            SMI_LIST.append(smi)
            Unopt_pols.append(pol)
        return POL_LIST, SMI_LIST, Unopt_pols
    else:
        smi = createPolymerSMILES(i,n,r,t)
        if verbosity:
            print(f'Polymer interpreted as: {i} {n} * {r} {t}')
            print(f"This gives the following SMILES: {smi}")
        pol_h,pol = optPol(smi) #both are mol objects
        return pol_h, smi, pol

def drawPol(pol, drawName):
    if type(pol) == list: #save a grid image instead
        img=Chem.Draw.MolsToGridImage(pol, legends = [f"n = {i+1}" for i, mol in enumerate(pol)], subImgSize = (250, 250))
        img.save(drawName)
    else:
        Chem.Draw.MolToFile(pol, drawName)

def write_or_read_pol(name, *, verbosity=False, read=False, mol=None):
    ext=name.split(".")[1]
    if read:
        if os.path.exists(name):
            if ext == "pdb":
                pol_h = Chem.MolFromPDBFile(name)
            elif ext == "mol":
                pol_h = Chem.MolFromMolFile(name)
            else:
                print(f"unsuported extention: {ext} Please use .pdb, or .mol") #.xyz cannot be read by rdkit.
                quit()
            polSMILES = Chem.MolToSmiles(pol_h)
            pol = Chem.MolFromSmiles(polSMILES)
            return pol_h, polSMILES, pol
        else:
            raise FileNotFoundError(name)
    else:
        if verbosity:
            print(f'attempting to save molecule to {name}')
        if ext == "xyz":
            Chem.MolToXYZFile(mol, name)
        elif ext == "pdb":
            Chem.MolToPDBFile(mol, name)
        elif ext == "mol":
            Chem.MolToMolFile(mol, name)
        else:
            print(f"Unsuported extention: {ext} Please use .pdb, .xyz or .mol")
            quit()
        if verbosity:
            print(f'Success')

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
    RG = Chem.rdMolDescriptors.CalcRadiusOfGyration(pol_h)
    #Chem.Descriptors3D.RadiusOfGyration(pol_h)
    #both seem to give identical results based on "SMILES to Rg.ipynb"
    return RG

def MolVolume(pol_h):
    MV = Chem.AllChem.ComputeMolVolume(pol_h, confId=-1, gridSpacing=0.2, boxMargin=2.0)
    return MV

def doCalcs(pol_h, calcs):
    #the variable calcs is a set
    #Calcs are only done if requested.
    data = {}
    if "SA" in calcs or "MHP" in calcs or "XMHP" in calcs:
        sasa = Sasa(pol_h)
        if not "XMHP" in calcs:
            data["SA"] = sasa
        calcs.discard("SA")
    if "LogP" in calcs or "MHP" in calcs or "XMHP" in calcs:
        logP = LogP(pol_h)
        if not "XMHP" in calcs:
            data["LogP"] = logP
        calcs.discard("LogP")
    if "RG" in calcs:
        rg = RadGyration(pol_h)
        data["RG"] = rg
        calcs.discard("RG")
    if "MV" in calcs:
        mv  =  MolVolume(pol_h)
        data["MV"] = mv
        calcs.discard("MV")
    if "MHP" in calcs or "XMHP" in calcs:
        mhp = logP / sasa
        data["MHP"] = mhp
        calcs.discard("MHP")
        calcs.discard("XMHP")
    if len(calcs) > 0:
        print(f"Unrecognized calculation(s): {calcs}. Use SA, LogP, MV, MHP, XMHP or RG")
    return data

def makePlot(pol_list, calculations, smiles_list, *, verbosity=False):
    dicts = []
    for i, pol in enumerate(pol_list):
        calcs = set(calculations)
        pol_data = doCalcs(pol,calcs)
        pol_data["N"] = i + 1
        pol_data["smi"] = smiles_list[i]
        dicts.append(pol_data)
    data = {k: [d[k] for d in dicts] for k in dicts[0]}
    
    ncols = len(data) - 2

    if ncols == 1: #matplotlib got angry at me for trying to make a plot with only one subplot. Use plt.plot to avoid this.
        calc_key = [k if k != "XMHP" else "MHP" for k in calculations][0] #use given calc as key unless XMHP, then use MHP.
        plt.plot(data["N"], data[calc_key],'o')
        plt.title(f'{calc_key} vs n')
        plt.xlabel('n') 
        plt.ylabel(f'{calc_key}')
    else:
        figure, axis = plt.subplots(ncols = ncols)
        series = 0
        for key in data:
            if key != "N" and key != "smi":
                axis[series].scatter(data["N"], data[key])
                axis[series].set_title(f"{key} vs n")
                series += 1
    figname = "Size-dependent-stats.png"
    plt.savefig(figname, bbox_inches = 'tight')
    print(f'Saved plot to {figname}')
    if verbosity:
        plt.show()
    return data, dicts

def exportToCSV(exptName, data, dicts_list, verbosity=False):
    with open(exptName, "w", newline = "") as c:
        cols = list(data.keys())
        writer = csv.DictWriter(c, fieldnames = cols)
        writer.writeheader()
        writer.writerows(dicts_list)
    print("Done exporting data to .csv file.")
    if verbosity:
        print(data)

def main():
    run_list = getArgs()

    for vardict in run_list:
        if vardict["read"] is None: #then get polymer parameters from CLI arguments.
            repeat_unit = getRepeatUnit(vardict["single_monomer"], vardict["super_monomer"])
                        
            if vardict["plot"]:
                POL_LIST, SMI_LIST, UNOPT_POL_LIST = make_One_or_More_Polymers(vardict["initiator"], vardict["n"],
                                                        repeat_unit, vardict["terminator"], verbosity=vardict["verbose"], plot=vardict["plot"])
            else:
                pol_h, polSMILES, pol = make_One_or_More_Polymers(vardict["initiator"], vardict["n"],
                                                        repeat_unit, vardict["terminator"], verbosity=vardict["verbose"], plot=vardict["plot"])
        else: #get mol from file
            if vardict["plot"]:
                raise TypeError("You may not plot data read from a file.") #we should be able to check for other files with name convention "{name}_{n}.{ext}"
            pol_h, polSMILES, pol = write_or_read_pol(vardict["read"], read=True)
            #pol_h is the as-is (probably 3-D) structure of the molecule. pol is the 2D structure.

        #saving the polymer to a file.
        if vardict["file"] is not None: #technically nothing wrong with using this as a roundabout way of converting between filetypes                
            if vardict["plot"]:
                base = vardict["file"].split(".")[0]
                ext = vardict["file"].split(".")[1]
                for i, mol in enumerate(POL_LIST):
                    name = f"{base}_{i+1}.{ext}"
                    write_or_read_pol(name, mol=mol)
            else:
                write_or_read_pol(vardict["file"], mol=pol_h)

        #drawing a picture of the polymer.
        if vardict["plot"]:
            pol = UNOPT_POL_LIST #submit this list of mols for use in grid image.
        if vardict["draw"] is not None:
            drawName = f'{vardict["draw"].split(".")[0]}.png'
            drawPol(pol, drawName)
        else:
            if vardict["verbose"]:
                #produce image if increased verbosity is requested even if no name is set.
                print("Saving image to polymer.png by default.")
                drawPol(pol, "polymer.png")

        #CALCULATIONS
        if vardict["verbose"]:
            print(f'requested calculations are {vardict["calculation"]}')
        if vardict["calculation"] is not None:
            if not vardict["plot"]:
                calcs = set(vardict["calculation"])
                data = doCalcs(pol_h, calcs) #use set to remove duplicates
                print(data)
                data["N"] = vardict["n"]
                data["smi"] = polSMILES
                dicts = [data]
                print(data)
            else:
                data, dicts = makePlot(POL_LIST, vardict["calculation"], SMI_LIST, verbosity=vardict["verbose"])
                
            if vardict["export"] is not None:
                exportToCSV(vardict["export"], data, dicts, verbosity=vardict["verbose"])

        print("\n") #separating runs visually if more than one.

if __name__ == "__main__":
    main()
