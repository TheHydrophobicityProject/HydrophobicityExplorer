from functools import cache
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, rdFreeSASA
import argparse, os, csv, json
import matplotlib.pyplot as plt
from smiles import monomer_dict, init_dict

def getStaticSettings():
    with open("settings.json", "r") as S: #this is where many defaults are set so they can easily be changed.
        settings_dict = json.load(S)
        return settings_dict

def getJsonArgs(jsonFile, dict):
    with open(jsonFile, 'r') as J: #open json file
        runs_dict = json.load(J) #read it
        for run in runs_dict["runs"]: #so few items the nested for loops shouldn't be a big deal
            run_keys = run.keys()
            for dict_key in dict.keys(): #the keys submitted to the func (derrived from CLI arguments)
                if dict_key not in run_keys: #if there is a key provided by user not in the dict derrived from the json file
                    run[dict_key] = dict[dict_key] #then add it
    run_list = runs_dict["runs"] #now we have a list of runs with all arguments from file and command line.
    return run_list

def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type = int, help = "The number of monomer or comonomer blocks.")
    parser.add_argument("-i", "--initiator", type = str, default = "Hydrogen", help = "Initiator Key from initiator dict or SMILES. Defaults to Hydrogen.")
    parser.add_argument("-t", "--terminator", type = str, default = "Hydrogen", help = "Terminator key taken from initiator dict or SMILES. Defaults to Hydrogen.")
    parser.add_argument("-m","--single_monomer", type = str, help = "Monomer key from the included monomer dict. See the -s flag for specifying a monomer that is not included.")
    parser.add_argument("-b", "--comonomer_block", type = str, nargs = '*',
                        help = "A series of space-separated monomer SMILES arranged in their repeating sequence. You can add an int preceeding any monomer to represent multiple copies of that monomer. e.g. 2 A B means AAB is the repeating super-monomer. Use quotes surrounding SMILES with problematic characters like = or ()")
    parser.add_argument("-d", "--draw", type = str, help = "Filename for polymer image.")
    parser.add_argument("-v", "--verbose", default = False, action = "store_true", help = "Set increased verbosity. Will draw polymer to polymer.png unless alternate name set by -d option.")
    parser.add_argument("-c","--calculation", type = str, nargs = '*', 
                        help = "Type of calculation(s) to be performed input as a space-separated list. Options are LogP, SA (surface area), MV (Molecular Volume), MHP (Mathers Hydrophobicity Parameter (LogP/SA; each of which will also be reported. Use XMHP to exclude those plots)) and RG (radius of gyration).")
    parser.add_argument("-f","--file", type = str, help = "The name/path of the file you wish to save the mol to. Supported formats are .pdb, .xyz and .mol")
    parser.add_argument("-r", "--read", type = str, help = "The name/path to file you wish to import. Supported formats are .pdb, .mol") #and .sdf")
    parser.add_argument("-p", "--plot", default = False, action = "store_true", 
                        help = "Include this option to generate a plot of whatever calculations are specified with -c on polymers from 1 to the n specified with the -n flag. This means the molecule cannot be read from a file with the -r flag. If used with the -f flag multiple files will be saved with names based off the one provided.")
    parser.add_argument("-e", "--export", type = str, help = "Include this option to export a .csv file of all data calculations. Specify the name here.")
    parser.add_argument("-j", "--json", type = str, help = "The path to a compatible .json file with any of the above arguments.")
    parser.add_argument("-q", "--quiet", default = False, action = "store_true", help = "Add this option to suppress the confirmation step which by default prevents calculations from running until the structure of the polymer is approved.")
    args = parser.parse_args()
    #get additional arguments from json file if provided or by default if no args provided.
    vardict = vars(args)
    if args.json is not None:
        run_list = getJsonArgs(args.json, vardict)
    else:
        run_list = [vardict] #args assumed to be in list later due to above.
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
    init = smiles_inators [0]
    term = smiles_inators[1]
    return init, term

def validate_end_group(inator, *, Init=False, Term=False, verbosity=False):
    if not Init and not Term:
        raise ValueError("Need to specify wether end group is terminator or initiator.")
    
    if Init:
        idx = 0 #look at first char of initiator
    else:
        idx = -1 #look at last character of terminator

    if inator != "" and inator[idx] == "*": #the attatchment point does not face the rest of polymer
        if verbosity:
            print("initiator smiles in wrong direction. Converting to mol object.")
        inator = Chem.MolFromSmiles(inator)
    elif inator != inator[::-1] and "*" not in inator:
        raise ValueError("end group smiles is not palendromic and has no attatchment point specified.")
    else:
        inator = inator.replace("*", "") #remove asterisk if not using rdkit method

    return inator

def get_building_blocks(i,t,m,*, verbosity = False):
    init, term = inator_smi_lookup(i, t)
    
    if type(m) == list:
        #replace any dict keys with corresponding smiles.
        deciphered_dict_keys = [monomer_dict[x] if x in monomer_dict else x for x in m]
        
        #we need an accurate count for number of monomers since a grouping specified by -s can be AABBB.
        #In this example, 1 unit of n is really 5 monomers. We want proper notation in figures.
        explicit_coefs = [int(x) if len(str(x)) == 1 else 0 for x in deciphered_dict_keys] #all monomers are (hopefully) 2 atoms or more. Assume others are coefs.
        sum_implicit_coefs = len(deciphered_dict_keys) - len(explicit_coefs) #number of implicit coefs of 1.
        monomers_per_n = sum(explicit_coefs) + sum_implicit_coefs #total monomer unit count per super-monomer.

        #start with empty repeat unit and concatonate stuff we find in the list.
        repeat_unit = ""
        #ommission of a coeficient implies 1 copy
        repeat_coef = 1
        for element in deciphered_dict_keys:
            try:
                repeat_coef = int(element) #is this a string of an integer?
            except:
                repeat_unit += repeat_coef * element #if not, repeat the SMILES as many times as specified (or once if no coef. provided).
                repeat_coef = 1 #reset coef.
    else:
        repeat_unit = monomer_smi_lookup(m) #if not a list, look for the corresponding smiles in the dictionary, will throw error if not included.
        monomers_per_n = 1

    init = validate_end_group(init, Init=True, verbosity=verbosity)
    term = validate_end_group(term, Term=True, verbosity=verbosity)
    
    return init, term, repeat_unit, monomers_per_n

def attatch_frags(polymer_smiles, *, add_initiator = (False, None), add_terminator = (False, None)): #the initiator and terminator are the kwargs
    pol = Chem.MolFromSmiles(polymer_smiles)
    #get indicies of fake atoms ("*")
    fake_atoms = [a.GetIdx() for a in pol.GetAtoms() if a.GetAtomicNum() == 0]
    #and their neighbors (to which we will actually be attatching.)
    conn_atoms = [pol.GetAtomWithIdx(x).GetNeighbors()[0].GetIdx() for x in fake_atoms]

    #lable the head and tail, accounting for possible absense of one or both inators.
    inators = []
    if add_initiator[0]:
        head = pol.GetAtomWithIdx(conn_atoms[0])
        head.SetProp("atomNote", "head")
        inators.append(add_initiator[1])
        if add_terminator[0]:
            tail = pol.GetAtomWithIdx(conn_atoms[1])
            tail.SetProp("atomNote", "tail")
            inators.append(add_terminator[1])
    elif add_terminator[0]:
        tail = pol.GetAtomWithIdx(conn_atoms[0])
        tail.SetProp("atomNote", "tail")
        inators.append(add_terminator[1])
    else:
        raise Exception(f"Unknown combination of inators {add_initiator = }, {add_terminator = }.")

    #set name to what will be used after one loop completes.
    mergedrw = pol
    for inator in inators:
        #see above.
        fake_atoms = [a.GetIdx() for a in inator.GetAtoms() if a.GetAtomicNum() == 0]
        #this time we just isolate atom object instead of index.
        attatch = [inator.GetAtomWithIdx(x).GetNeighbors()[0] for x in fake_atoms][0]
        #lable.
        attatch.SetProp("atomNote", "attatch")
        #put the two mols into the same object (still no bond between them.)
        merged = Chem.CombineMols(inator, mergedrw)
        #change to rwmol object which can be changed.
        mergedrw = Chem.RWMol(merged)

        #indicies of atoms with notes
        attachments = [a.GetIdx() for a in mergedrw.GetAtoms() if a.HasProp('atomNote')]
        #isolating proper index from list to use in bond formation.
        inator_attatchment = [i for i in attachments if mergedrw.GetAtomWithIdx(i).GetProp('atomNote') == "attatch"][0]

        if inator == add_initiator[1]: #inator changes with each pass through the loop
            bond_here = [i for i in attachments if mergedrw.GetAtomWithIdx(i).GetProp('atomNote') == "head"][0]

        if inator == add_terminator[1]:
            bond_here = [i for i in attachments if mergedrw.GetAtomWithIdx(i).GetProp('atomNote') == "tail"][0]
            
        #make bond
        mergedrw.AddBond(bond_here, inator_attatchment, Chem.rdchem.BondType.SINGLE)
        #change label so that atom is not targeted a second time for bond formation.
        mergedrw.GetAtomWithIdx(inator_attatchment).ClearProp('atomNote')

    #count up number of dummy atoms ("*")
    dummies = [a for a in mergedrw.GetAtoms() if a.GetAtomicNum() == 0]
    numDummies = len(dummies)

    #remove the dummy atoms (need to do one at a time because index changes each time).
    for i in range(numDummies):
        mergedrw.RemoveAtom([a.GetIdx() for a in mergedrw.GetAtoms() if a.GetAtomicNum() == 0][0])

    smi = Chem.MolToSmiles(mergedrw)
    return smi

def add_inator_smiles(smi, init, term, *, verbosity=False):
    if verbosity:
        print(f"polymer smiles is {smi} before any end groups")

    if type(init) != str: #i.e. a mol object instead
        smi = "*" + smi
        add_initiator = True # we will attatch with mol-based methods
    else:
        add_initiator = False
        smi = init + smi #just attatch as string instead
        if verbosity and init != "":
            print(f"polymer smiles is {smi} after adding initiator smiles")
            
    if type(term) != str: #i.e. a mol object instead
        smi = smi + "*" #same as above but for terminator. Attachment point is at end this time.
        add_terminator = True
    else:
        add_terminator = False
        smi = smi + term #same as above but for terminator
        if verbosity and init != "": #no need to print this extra info if adding "" (Hydrogen)
            print(f"polymer smiles is {smi} after adding terminator smiles")

    if add_terminator or add_initiator:
        if verbosity:
            print(f"converting polymer body {smi} to mol object to add frags")
        smi = attatch_frags(smi, add_initiator=(add_initiator, init), add_terminator=(add_terminator, term))

    return smi

def createPolymerSMILES(i,n,r,t,*, verbosity = False, test = False):
    init, term, repeat_unit, m_per_n = get_building_blocks(i,t,r, verbosity=verbosity)

    polymer_SMILES = n * repeat_unit
    
    if test: # a parameter used to generate an n=1 image where it is easy to see where end groups attatch
        #if you don't do this and have n=15, the image is very hard to parse visually and some parts of pol will overlap.
        test_smi = repeat_unit
        test_smi = add_inator_smiles(test_smi, init, term, verbosity=verbosity)
        verbosity = False #turn off verbosity for the next generation because we already display info about endgroup connections the first time.
    
    full_smiles = add_inator_smiles(polymer_SMILES, init, term, verbosity=verbosity)

    if test:
        return test_smi, full_smiles, m_per_n
        #return test smiles too so it can be previewed. It is fast to make both before confirmation
        #but we do the confirmation before optimizing geometry.
    else:
        return full_smiles, m_per_n
   
def optPol(smiles, *, name=None, nConfs=5, threads=0, iters=1500): #name is provided my supplemental scripts.
    #make Mol object:
    pol = Chem.MolFromSmiles(smiles)
    #check mol
    Chem.SanitizeMol(pol)
    #opt steps
    pol_h = Chem.AddHs(pol)
    #random coords lead to better geometries than using the rules rdkit has. Excluding this kwarg leads to polymers that do not fold properly.
    ids = AllChem.EmbedMultipleConfs(pol_h, numConfs=nConfs, useRandomCoords=True, numThreads=threads) #get multiple conformers for better stats
    # AllChem.EmbedMolecule(pol_h, useRandomCoords=True) 
    # AllChem.MMFFOptimizeMolecule(pol_h, maxIters=2500) 
    touple_list = AllChem.MMFFOptimizeMoleculeConfs(pol_h, numThreads=threads, maxIters=iters) #default 200 iterations.
    for i, tup in enumerate(touple_list):
        if tup[0] == 1:
            pol_h.RemoveConformer(i)
            # print(f"removing conf {i}")
    if pol_h.GetNumConformers() == 0: #tell the user to change something if none of the confs are good.
        raise Exception("Optimization failed to converge. Rereun with higher maxIters.")
    
    #calculations are inconsistent if using conf ids instead of just single-conf mols. Translate to sdf mol supplier to make it easy to integrate with reading files.
    if name is None:
        sdfFilename = "tmp.sdf"
    else:
        ext = name.split(".")[1]
        if ext != "sdf":
            raise Exception("Filename must use .sdf format.")
        sdfFilename = name
    writer = Chem.SDWriter(sdfFilename)
    for conf in pol_h.GetConformers(): #loop through all conformers that still exist. We only write the conformations that converged.
        cid = conf.GetId() #The numbers may not be sequential.
        pol_h.SetProp('_Name', f'conformer_{cid}') #when sdf is read each conf is separate mol object.
        # pol_h.SetProp('ID', f'conformer_{cid}') #Similar method can be used to print number of monomers for plot jobs.
        writer.write(pol_h, confId=cid)
    
    suppl = Chem.SDMolSupplier(sdfFilename) #itterator that has all mols in the sdf file.
    
    if name is None:
        os.remove(sdfFilename)

    return pol, suppl #suppl now has each conformation as a separate mol obj when we iterate thru it.

def confirmStructure(smi, *, proceed=None):
    #save image to temporary file
    drawPol(Chem.MolFromSmiles(smi), "tmp_confirm.png")
    img = Image.open("tmp_confirm.png")
    #show it to user
    img.show()
    inp = input("Does this look right? [Y/n]")
    
    if os.path.exists("tmp_confirm.png"):
        os.remove("tmp_confirm.png")
        #delete the file

    #affirmation is y, Y or just hitting enter
    if inp.lower() == "y" or inp == "":
        inp = True
        print("Great! If you wish to bypass this confirmation step, use the -q flag when running this script.")
    else:
        # inp = False #not actually used since program quits.
        print("Please try adjusting input and try again.")
        quit()
        #aborts so user can retry

    if proceed is not None:
        return inp #used to stop plotting jobs from asking for confirmation for each pol those jobs generate.

def make_One_or_More_Polymers(i, n, r, t, *, verbosity=False, plot=False, confirm=False):
    POL_LIST = []
    SMI_LIST = []
    Unopt_pols = []
    if plot: #make molecules from n=1 to n specified by user.
        N_array = range(1, n+1)
        #this allows us to confirm only once for plotting jobs
        if confirm == True:
            proceed = False
        else:
            proceed = True

        for j in N_array:
            if j == 1 and confirm and not proceed:
                test_smi, smi, m_per_n = createPolymerSMILES(i,j,r,t, verbosity=verbosity, test=True)
                verbosity = False
                proceed = confirmStructure(test_smi, proceed=proceed)
            
            if j > 1 or not confirm: #do not test if j is large or if we ask not to test at all.
                smi, m_per_n = createPolymerSMILES(i, j, r, t, verbosity=verbosity)
            
            if verbosity:
                print(f"Done generating SMILES with n = {j} now: {smi}")
                print("Converting to RDkit mol now.")
            #get opt and unopt molecules.
            pol, pol_h = optPol(smi)
            POL_LIST.append(pol_h)
            Unopt_pols.append(pol)
            #save smiles
            SMI_LIST.append(smi)
        return POL_LIST, SMI_LIST, Unopt_pols, m_per_n
    else: #just one polymer.
        test_smi, full_smi, m_per_n = createPolymerSMILES(i, n, r, t, verbosity=verbosity, test=True)
        if verbosity:
            print(f'Polymer interpreted as: {i} {n} * {r} {t}')
            print(f"This gives the following SMILES: {full_smi}")

        if confirm:
            print("Showing structure with n=1 to confirm correct end groups")
            confirmStructure(test_smi)
        
        pol, pol_h = optPol(full_smi) #both are mol objects
        return pol_h, full_smi, pol, m_per_n

def drawPol(pol, drawName, *, mpn=1, image_size=250):
    if type(pol) == list: #save a grid image instead
        img = Chem.Draw.MolsToGridImage(pol, legends = [f"n = {(i + 1) * mpn}" for i in range(len(pol))], subImgSize=(image_size, image_size))
        #mpn is the number of monomers per "n". This is > 1 when -s is used and multiple monomers or copies of the same monomer are specified.
        img.save(drawName)
    else:
        Chem.Draw.MolToFile(pol, drawName)

def write_or_read_pol(name, *, verbosity=False, read=False, mol=None):
    ext = name.split(".")[1]
    if read:
        pol_h = None
        suppl = None
        if os.path.exists(name):
            #is the file type valid?
            if ext == "pdb":
                pol_h = Chem.MolFromPDBFile(name)
            elif ext == "mol":
                pol_h = Chem.MolFromMolFile(name)
            elif ext == "sdf":
                suppl = Chem.SDMolSupplier(name)
                for pol in suppl: #grab one conf so we can visualize
                    pol_h = pol
                    break #we break the loop since only one conf is needed to visualize
            else:
                print(f"unsuported extention: {ext} Please use .pdb, .mol or .sdf") #.xyz cannot be read by rdkit.
                quit()

            if suppl is None:
                sdf_name="tmp.sdf"
                writer = Chem.SDWriter(sdf_name)
                writer.write(pol_h)
                suppl = Chem.SDMolSupplier(sdf_name) #iterator that has all mols in the sdf file.
                os.remove(sdf_name)

            #convert to smiles so it can be visualized
            polSMILES = Chem.MolToSmiles(pol_h)
            pol = Chem.MolFromSmiles(polSMILES)
            #but visualization needs to come from unoptimized polymer. (could also do RemoveAllConformers() but we still need smiles anyway.)
            return suppl, polSMILES, pol #These are used for calcs, smiles part of csv and 2D visualization.
        else:
            raise FileNotFoundError(name)
    else: #i.e. write file.
        if verbosity:
            print(f'attempting to save molecule to {name}')

        #what are we dealing with?
        if type(mol) == type(Chem.SDMolSupplier()):
            suppl = True
            for pol in mol:
                first_conf = pol
                cid = -1
                break #we will only be writing the first conf in the iterable for non-sdf files.
        else:
            suppl = False
            for conf in mol.GetConformers():
                first_conf = mol #keep name convention consistent.
                cid = conf.GetId()
                break #we will only be writing the first conf in the iterable for non-sdf files.

        #is the file type valid?
        if ext == "sdf":
            writer = Chem.SDWriter(name)
            if suppl:
                for pol in mol:
                  writer.write(pol)
            else: #is a mol opject        
                for conf in mol.GetConformers(): #loop through all conformers that still exist. We only write the conformations that converged.
                    cid = conf.GetId() #The numbers may not be sequential.
                    mol.SetProp('_Name', f'conformer_{cid}') #when sdf is read each conf is separate mol object.
                    # mol.SetProp('ID', f'conformer_{cid}') #Similar method can be used to print number of monomers for plot jobs.
                    writer.write(mol, confId=cid)
        
        elif ext == "xyz":
            Chem.MolToXYZFile(first_conf, name, confId = cid)
        elif ext == "pdb":
            Chem.MolToPDBFile(first_conf, name, confId = cid)
        elif ext == "mol":
            Chem.MolToMolFile(first_conf, name, confId = cid)
        
        else:
            print(f"Unsuported extention: {ext} Please use .pdb, .xyz or .mol")
            quit()

        if verbosity:
            print(f'Success')

def avg_stat(list_of_stats):
    return sum(list_of_stats) / len(list_of_stats)

def Sasa(pol_h):
    # Calculate SASA
    sasa_lst = []

    for mol in pol_h: #calculate desired property for each molecule in the file.
        radii = Chem.rdFreeSASA.classifyAtoms(mol)
        sasa = Chem.rdFreeSASA.CalcSASA(mol, radii)
        sasa_lst.append(sasa)

    return avg_stat(sasa_lst)

def LogP(pol_h):
    logP_lst = []
    for mol in pol_h:
        logP_lst.append(Chem.Descriptors.MolLogP(mol))

    return avg_stat(logP_lst)

def RadGyration(pol_h):
    rg_list = []

    for mol in pol_h:
        rg_list.append(Chem.rdMolDescriptors.CalcRadiusOfGyration(mol))
    #Chem.Descriptors3D.RadiusOfGyration(pol_h)
    #both seem to give identical results based on "SMILES to Rg.ipynb"
    return avg_stat(rg_list)

def MolVolume(pol_h, *, grid_spacing=0.2, box_margin=2.0):
    mv_list = []
    for mol in pol_h:
        mv_list.append(Chem.AllChem.ComputeMolVolume(mol, gridSpacing=grid_spacing, boxMargin=box_margin))

    return avg_stat(mv_list)

def doCalcs(pol_h, calcs):
    #The type of variable /calcs/ is set
    #Calcs are only done if requested.
    #remove entries from set after each calculation and print the unrecognized ones at the end.
    data = {}
    if "SA" in calcs or "MHP" in calcs or "XMHP" in calcs:
        sasa = Sasa(pol_h)
        if not "XMHP" in calcs: #if XMHP is included user eXcluisively wants MHP, so we don't return this data.
            data["SA"] = sasa
        calcs.discard("SA")
    if "LogP" in calcs or "MHP" in calcs or "XMHP" in calcs:
        logP = LogP(pol_h)
        if not "XMHP" in calcs: #if XMHP is included user eXcluisively wants MHP, so we don't return this data.
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
        data["LogP/SA"] = mhp
        calcs.discard("MHP")
        calcs.discard("XMHP")
    if len(calcs) > 0:
        print(f"Unrecognized calculation(s): {calcs}. Use SA, LogP, MV, MHP, XMHP or RG")
    return data

def makePlot(pol_list, calculations, smiles_list, *, verbosity=False, mpn=1, data_marker='o', fig_filename="Size-dependent-stats.png"):
    units = { "LogP/SA":"Angstroms^-2", "LogP":"", "RG":"Angstroms", "SA":"Angstroms^2", "MV":"Molar Volume" }
    dicts = []
    for i, pol in enumerate(pol_list):
        calcs = set(calculations)
        pol_data = doCalcs(pol, calcs)
        pol_data["N"] = (i + 1) * mpn
        pol_data["smi"] = smiles_list[i]
        dicts.append(pol_data)
    data = {k: [d[k] for d in dicts] for k in dicts[0]}
    
    ncols = len(data) - 2

    if ncols == 1: #matplotlib got angry at me for trying to make a plot with only one subplot. Use plt.plot to avoid this.
        calc_key = [k if k != "XMHP" or k != "MHP" else "LogP/SA" for k in data.keys()][0] #use given calc as key unless XMHP, then use MHP.
        plt.plot(data["N"], data[calc_key], data_marker)
        plt.title(f'{calc_key} vs n')
        plt.xlabel('n') 
        plt.ylabel(f'{calc_key} ({units[calc_key]})')
    else:
        #need to make multiple subplots if multiple calcs requested.
        figure, axis = plt.subplots(ncols = ncols)
        series = 0
        for key in data:
            #we can't plot N vs N or anything to do with smiles
            if key != "N" and key != "smi":
                axis[series].scatter(data["N"], data[key])
                axis[series].set_title(f"{key} vs n")
                series += 1
    figname = fig_filename
    plt.savefig(figname, bbox_inches = 'tight')
    print(f'Saved plot to {figname}')
    if verbosity:
        print(data)
        plt.show()
    return data, dicts

def exportToCSV(exptName, data, dicts_list, verbosity=False):
    with open(exptName, "w", newline = "") as c:
        #set column names as dict keys
        cols = list(data.keys())
        writer = csv.DictWriter(c, fieldnames = cols)
        writer.writeheader()
        #write the data.
        writer.writerows(dicts_list)
    print(f"Done exporting data to {exptName}.")
    if verbosity: #this is turned off by main() if plotting is also turned on since both functions can print data and that is only needed once.
        print(data)

def main():
    default_dict = getStaticSettings()
    run_list = getArgs()

    for vardict in run_list:
        if vardict["read"] is None: #then get polymer parameters from CLI arguments.
            repeat_unit = getRepeatUnit(vardict["single_monomer"], vardict["comonomer_block"])
                        
            if vardict["plot"]:
                POL_LIST, SMI_LIST, UNOPT_POL_LIST, M_PER_N = make_One_or_More_Polymers(vardict["initiator"], vardict["n"],
                    repeat_unit, vardict["terminator"], verbosity=vardict["verbose"], plot=vardict["plot"], confirm = not vardict["quiet"])
            else:
                pol_h, polSMILES, pol, M_PER_N = make_One_or_More_Polymers(vardict["initiator"], vardict["n"],
                    repeat_unit, vardict["terminator"], verbosity=vardict["verbose"], plot=vardict["plot"], confirm = not vardict["quiet"])
        else: #get mol from file
            if vardict["plot"]:
                raise TypeError("You may not plot data read from a file.") #we should be able to check for other files with name convention "{name}_{n}.{ext}"
            pol_h, polSMILES, pol = write_or_read_pol(vardict["read"], read=True)
            M_PER_N = 1
            #pol_h is the as-is (probably 3-D) structure of the molecule. pol is the 2D structure.

        #saving the polymer to a file.
        if vardict["file"] is not None: #technically nothing wrong with using this as a roundabout way of converting between filetypes                
            if vardict["plot"]:
                base = vardict["file"].split(".")[0]
                ext = vardict["file"].split(".")[1]
                for i, mol in enumerate(POL_LIST):
                    name = f"{base}_{i + 1}.{ext}"
                    write_or_read_pol(name, mol=mol)
            else:
                write_or_read_pol(vardict["file"], mol=pol_h)

        #drawing a picture of the polymer.
        if vardict["plot"]:
            pol = UNOPT_POL_LIST #submit this list of mols for use in grid image.
        if vardict["draw"] is not None:
            drawName = f'{vardict["draw"].split(".")[0]}.png'
            drawPol(pol, drawName, mpn=M_PER_N)
        else:
            if vardict["verbose"]:
                defaultName = default_dict["drawing_default"]
                #produce image if increased verbosity is requested even if no name is set.
                print(f"Saving image to {defaultName} by default.")
                drawPol(pol, defaultName, mpn=M_PER_N)

        #CALCULATIONS
        if vardict["verbose"]:
            print(f'requested calculations are {vardict["calculation"]}')
        if vardict["calculation"] is not None:
            if not vardict["plot"]:
                calcs = set(vardict["calculation"])
                data = doCalcs(pol_h, calcs) #use set to remove duplicates
                data["N"] = vardict["n"] * M_PER_N
                data["smi"] = polSMILES
                dicts = [data]
                print(data)
            else:
                data, dicts = makePlot(POL_LIST, vardict["calculation"], SMI_LIST, verbosity=vardict["verbose"], mpn=M_PER_N)
                
            if vardict["export"] is not None:
                if vardict["plot"]: #we don't need to print data twice if both -p and -e use verbosity=True
                    verbo = False
                else:
                    verbo = vardict["verbose"]        
                exportToCSV(vardict["export"], data, dicts, verbosity=verbo)

        print("\n") #separating runs visually if more than one.

if __name__ == "__main__":
    main()
