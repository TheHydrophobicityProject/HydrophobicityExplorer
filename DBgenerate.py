# DO NOT UPLOAD THIS PROJECT
import pymongo
from pymongo import MongoClient
import concurrent.futures
import time
import random
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdFreeSASA
import rdkit
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import xlsxwriter


database_IP = "192.168.41.8"
database_port = 27017

# Get the ID of the best conformer given a molecule with conformers and a set of conformer IDs
# Note that mol_h MUST be pre-populated with the desired conformers
def GetBestConformerID(mol_h, ids, returnMinusOne=False):
    best = []
    for id in ids:
        prop = AllChem.MMFFGetMoleculeProperties(mol_h)
        ff = AllChem.MMFFGetMoleculeForceField(mol_h, prop, confId=id)
        ff.Minimize()
        en = float(ff.CalcEnergy())
        econf = (en, id)
        best.append(econf)
    best.sort()
    print("Found "+str(len(best))+" conformations.")
    if len(best) == 0 and returnMinusOne is True:
        return -1
    elif len(best) == 0:
        raise Exception("Error: No valid conformations found for the molecule. Try increasing the number of conformations.")
    # The best conformer is the first tuple and the ID of that conformer is the second value in the tuple
    best_id = int(best[0][1])
    # Return the best ID
    return best_id

# Utility method to save a conformer to a file
# Accepts a molecule, conformer ID, and filename
# Note that this function has a HARDCODED PREFIX FOLDER that it saves files inside of.
def SaveConformationToFile(molecule, confId=-1, filename="scratch"):
    # Hardcoded prefix where the files are saved - change if necessary
    prefix = "molecules/"
    if confId == -1:
        raise Exception("You must provide a conformation ID to use this function.")
    # Get the conformation with the given ID
    best_conf = Chem.MolToMolBlock(molecule,confId=confId)
    # Open the file, write, then close the file.
    f = open(prefix + filename, 'w+')
    # MAKE SURE this doesn't just append to the file over and over again, which would be bad
    f.write(best_conf)
    f.close()
    print("Successfully wrote conformation to file.")
    return prefix + filename


# Given a SMILES string and other parameters, this function generates the required number of conformers and gets the LogP
# and SASA for the best one
# Calculates logP for the OVERALL MOLECULE since I couldn't work out how to do it for a conformer
# DANGER - the output of this function is NOT 100% consistent; if you get your LogP and SASA from two different runs
# you WILL run into problems. Get ALL your values at the SAME TIME.


# Initiate a new collection to the conformers DB - should be thread safe
def NewConformersDBConnection():
    client = MongoClient(database_IP,port=database_port)
    collection = client.cocalc.conformers
    return collection

# Copy from here

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# Utility functions for database version
# These should DEFINITELY get moved to another file; they are only here for debugging


# Thread safe function to add conformers to database
def AsyncAddConformerToDB(conformer):
    # Set up the print prefix
    print_prefix = conformer["initiatorName"]+" "+conformer["repeatUnitName"]+" "+str(conformer["numRepeatUnits"])+"-mer: "
    #print_prefix = "" #print_prefix = "("+conformer["initiatorName"]+\
#    "-"+conformer["monomerName"]+"-"+str(conformer["numRepeatUnits"])+"mer): "
    #print(print_prefix+"Began database update")
    # Set the update check to false
    updateDB = False
    # open a database connection
    conformersConnector = MongoClient(database_IP, port=database_port)
    cocalcDB = conformersConnector.cocalc
    conformerCollection = cocalcDB.conformers
    # check to make sure the conformer isn't already in the DB
    #print("Connected to database")
    numResults = conformerCollection.count_documents({'smiles': conformer['smiles'],
                               'precision': conformer['precision']})
    if numResults > 0:
        print(print_prefix+"Found existing result")
        existingResult = True
        # We found a result - get the DB result to check the highest version
        result = conformerCollection.find_one({'smiles': conformer['smiles'],
                         'precision': conformer['precision']},
                                              sort=[('version', pymongo.DESCENDING)])
        if result["version"] is not None:
            # result found with non-null version, check if it needs updating
            if result['version'] < conformer['version']:
                # DB result is lower version, so we should update it
                print("DB value is older version, updating")
                updateDB = True
            else:
                # DB result is equal or higher version, don't update
                print(print_prefix+"DB value is already up to date.")
        else:
            # If the version in the DB is null, always update
            updateDB = True
    else:
        print(print_prefix+"No existing result")
        # No results, update the DB
        updateDB = True
        existingResult = False

    # Update database if needed
    #print("Reached updateDB evaluation")
    if updateDB == True:
        # Unfortunately I can't seem to get upsert to work with update_one...
        if existingResult == True:
            # Result already exists, update it
            #print("Reached database update_one")
            conformerCollection.update_one(
                {"smiles": conformer["smiles"],
                        "precision": conformer["precision"],
                        "version":conformer["version"]},
                {"$set": conformer})
            #print(print_prefix+"Updated database value")
        else:
            # Result does not already exist, insert it
            #print("Reached database insert_one")
            #print(conformer)
            conformerCollection.insert_one(conformer)
            print(print_prefix+"Inserted new database value")
    else:
        # No update necessary; do nothing
        print(print_prefix+"Database update not needed; skipping.")
    return 0

def SaveLogPAndSASAForAllConformersToDB(initiators, #initiator is a JSON dictionary
                               monomers, # molecule parameters
                               initiatorName=None,
                               start=1, stop=10, step=1, # range parameters for monomer generation
                               precision=100, # number of conformations to examine
                               version=1, randomSeed=1, numThreads = 1,
                               useExpTorsionAnglePrefs=True, useRandomCoords=True,
                                     returnMinusOneIfError=True, updateDB=True):
    all_polymers = GenerateAllPolymersInRange(initiators=initiators,monomers=monomers,start=start,stop=stop,step=step,
                                              precision=precision,randomSeed=randomSeed,numThreads=numThreads,
                                              useExpTorsionAnglePrefs=useExpTorsionAnglePrefs,
                                              useRandomCoords=useRandomCoords)
    polymers = DeduplicatePolymersWithDB(all_polymers)
    masterExecutionStartTime = time.time()
    # Set up the execution values
    polymersForExecution = []
    for polymer in polymers:
        # Set up the parameters needed by the asynchronous worker function
        parameters = {}
        parameters["numConfs"] = precision
        parameters["version"] = version
        parameters["randomSeed"] = randomSeed
        parameters["numThreads"] = 1  # DO NOT change this
        parameters["useRandomCoords"] = useRandomCoords
        parameters["useExpTorsionAnglePrefs"] = useExpTorsionAnglePrefs
        parameters["returnMinusOne"] = returnMinusOneIfError
        parameters["updateDB"] = updateDB
        # Now lump the parameters together to pass down to executor.submit
        # If you know how you're SUPPOSED to pass multiple arguments to executor.submit, please fix this.
        allparams = {}
        allparams["polymer"] = polymer
        allparams["parameters"] = parameters
        polymersForExecution.append(allparams)
    # Multithreading!
    with concurrent.futures.ThreadPoolExecutor(max_workers=numThreads) as executor:
    #    # polymers block used to go here
        for polymer in polymersForExecution:
            executor.submit(AsyncGetLogPAndSASAForBestConformer, polymer)
            #AsyncGetLogPAndSASAForBestConformer(polymer)
        print("Queued "+str(len(polymers))+" polymers for conformation searching.")
    print("Completed execution in "+str(time.time()-masterExecutionStartTime)+" seconds.")

# This function is NOT SAFE FOR CONCURRENT ACCESS
# If you use this WHILE YOU ARE RUNNING SOMETHING ELSE THAT UPDATES THE DATABASE
# you may encounter errors, so don't do that
def RemoveAlreadyExstingPolymers(polymers):
    # Set counter for number of polymers removed to zero
    numRemoved = 0
    numOriginalPolymers = len(polymers)
    # Open a new database connection
    conformersConnector = MongoClient(database_IP,port=database_port)
    collection = conformersConnector.cocalc.conformers
    # iterate over the polymers we are given
    numProcessed = 0
    #print("Length is "+str(len(polymers)))
    for polymer in polymers:
        numProcessed += 1
        # Look up whether there is an existing entry with same SMILES, precision, and version
        numResults = collection.count_documents({'smiles':polymer['smiles'],
                                  'precision':polymer['precision'],
                                  'version':{'$gte':int(polymer["version"])}})
        #print(numResults)
        if numResults > 0:
            #print("Found a DB result")
            # We found a result - get the highest database value
            result = collection.find_one({'smiles':polymer['smiles'],
                                  'precision':polymer['precision'],
                                  'version':{'$gte':int(polymer['version'])}},
                                         sort=[('version', pymongo.DESCENDING)])
            #print("DB check completed")
            # check the version
            if result['version'] >= polymer['version']:
                # DB result is >= what we've been asked for, so drop this one
                polymers.remove(polymer)
                # Update removal counter
                #print("DB version is not lower")
                numRemoved += 1
            else:
                # DB version is not high enough; keep this one
                print("DB version for polymer with SMILES:"+
                      polymer['smiles']+"and precision "+str(polymer['precision']+" is lower; will update"))
        elif numResults == 0:
            # No results, keep this one
            print("Found 0 results, keeping")
            #pass
        else:
            # Something went wrong.
            print("Error retrieving database value for polymer "+polymer['smiles'])
            raise Exception("Error retrieving database value - returned invalid result count.")
        #print("Reached end of inside of for loop")
    # Now return the pruned LIST
    print("Pruning complete.")
    print("Removed "+str(numRemoved)+" pre-existing polymers of "+str(numOriginalPolymers)+" original polymers.")
    print("There are "+str(len(polymers))+" remaining polymers")
    #print(polymers)
    return polymers

# Given the described parameters, generates all valid polymers as a JSON dictionary
def GenerateAllPolymersInRange(initiators, #initiator is a JSON dictionary
                               monomers, # molecule parameters
                               initiatorName=None,
                               start=1, stop=10, step=1, # range parameters for monomer generation
                               precision=100, # number of conformations to examine
                               version=1, randomSeed=1, numThreads = 1,
                               useExpTorsionAnglePrefs=True, useRandomCoords=True # parameters to pass down
                               ):
    # Check which kind of value we were given for the initiators
    if isinstance(initiators, dict):
        # Dictionary - we will iterate over this
        initiatorIsDictionary = True
    elif isinstance(initiators, str):
        # String - need to convert to dictionary
        initiatorIsDictionary = False
        if initiatorName is not None:
            # If we are given an initiator name, pass it
            initiators = {initiatorName:initiators}
        else:
            # Otherwise set a backup value just in case
            initiators = {'NotProvided':initiators}
    else:
        # Error!
        raise Exception("You must provide either a dict or a string for the initiator.")
    polymers = []
    # Iterate over the list of possible initiators
    for initiatorName, initiatorSMILES in initiators.items():
        # For each initiator, iterate over the list of possible monomer types
        for monomerName, monomerSMILES in monomers.items():
            # For each monomer type, iterate over each possible length within the given parameters
            for monomerLength in range(start, stop+1, step):
                # Construct the JSON object for each individal polymer
                polymer = {}
                polymer["smiles"] = (initiatorSMILES + (monomerLength * monomerSMILES))
                polymer["precision"] = precision
                polymer["version"] = version
                polymer["initiator"] = initiatorSMILES
                if initiatorIsDictionary == True and initiatorName is not None:
                    polymer["initiatorName"] = initiatorName
                polymer["repeatUnit"] = monomerSMILES
                polymer["repeatUnitName"] = monomerName
                polymer["numRepeatUnits"] = monomerLength
                polymer["randomSeed"] = randomSeed
                # Now add it to the dictionary
                polymers.append(polymer)
    # Return the completed dictionary
    print("Generated "+str(len(polymers))+" polymers")
    #print(polymers)
    return polymers


# Experimental multithreaded version
# DO NOT pass numThreads greater than 1 in the parameters
def AsyncGetLogPAndSASAForBestConformer(allparams):
    # Start the timer
    executionStartTime = time.time()
    # Set error tracking value to False
    noValidConformers = False
    # Unpack the parameters since I couldn't work out how to pass more than one parameter to the ThreadPoolExecutor
    # change this (and the bit in the parent function) if you know how to do that
    polymer = allparams["polymer"]
    params = allparams["parameters"]
    # Set up the print prefix
    print_prefix = polymer["initiatorName"] + " " + polymer["repeatUnitName"] + " " + str(
        polymer["numRepeatUnits"]) + "-mer: "
    # Print that we have begun execution - for multithread troubleshooting
    print(print_prefix+"Began execution")
    # First process the SMILES string
    mol = Chem.MolFromSmiles(polymer["smiles"])
    mol_h = Chem.AddHs(mol)
    # Start conformer generation timer
    embedStartTime = time.time()
    # Now generate the conformations
    ids = AllChem.EmbedMultipleConfs(mol_h, numConfs=params["numConfs"], randomSeed=params["randomSeed"],
                                     useExpTorsionAnglePrefs=params["useExpTorsionAnglePrefs"],
                                     useRandomCoords=params["useRandomCoords"],
                                     numThreads=params["numThreads"])
    # Find the best one
    #print("Embedded the conformations")
    try:
        best_conf_id = GetBestConformerID(mol_h, ids)
    except Exception as e:
        if params["returnMinusOne"] is True:
            print(print_prefix+"Error: Found no valid conformations for " + polymer["smiles"] + ".")
            noValidConformers = True
        else:
            raise e
    print(print_prefix+"Best conformer ID is " + str(best_conf_id) + ".")
    # print the time taken for profiling
    print(print_prefix+"Embedded conformations in "+str(time.time()-embedStartTime)+" seconds.")
    if best_conf_id == -1 and params["returnMinusOne"] is True:
        return 0, 1
    
    radii = rdFreeSASA.classifyAtoms(mol_h)
    if noValidConformers is False:
        # Now calculate LogP and SASA
        # Calculate SASA based on the best conformer
        # classifyAtoms CRASHED when I tried it with , confIdx=best_conf_id
        # but someone needs to go back and make sure it's actually OK to use it without
        # and that that won't cause problems!
        #print(print_prefix+"Found at least one conformer")
        polymer["sasa"] = rdFreeSASA.CalcSASA(mol_h, radii, confIdx=best_conf_id)
        # LogP does NOT have an option to feed in a conformer so just calculate it for the overall molecule
        polymer["LogP"] = Descriptors.MolLogP(mol_h)
        polymer["LogPoverSA"] = polymer["LogP"]/polymer["sasa"]
        polymer["bestConformerID"] = best_conf_id
        polymer["bestConformerMolBlock"] = Chem.MolToMolBlock(mol_h,confId=best_conf_id)
        # Calculate radius of gyration
        polymer["Rg"] = Chem.Descriptors3D.RadiusOfGyration(mol_h, confId=best_conf_id)
        # Calculate the number of rotatable bonds
        polymer["numRotatableBonds"] = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol_h)
        # Track calculation time for analysis
        polymer["calcTime"] = (time.time() - executionStartTime)
    else:
        print(print_prefix+"Found 0 conformers")
        # No valid LogP/SASA, so set error state
        polymer["noConformersFound"] = True
        # And set sasa/LogP to placeholder values
        polymer["sasa"] = 1
        polymer["LogP"] = 0
    # Now for ALL molecules we want to save ALL the RDKit data, just in case
    # DO NOT DO THIS BELOW
    # polymer["allConfsRDKitMoleculeObjectWithH"] = mol_h
    # This is a pointer and the database HATES it
    # if you try to pass down a memory pointer like this the insert will silently fail
    # If database update was requested, do that
    if params["updateDB"] == True:
        AsyncAddConformerToDB(polymer)
    # All operations completed, so print the finish time
    print("Completed conformer search for " + polymer["initiatorName"] + "-" +
          str(polymer["numRepeatUnits"]) + "-" + polymer["repeatUnitName"] +
          " in " + str(time.time() - executionStartTime) + " seconds.")
    # Now return LogP and SASA
    return polymer["LogP"], polymer["sasa"]

# Get the ID of the best conformer given a molecule with conformers and a set of conformer IDs
# Note that mol_h MUST be pre-populated with the desired conformers
def GetBestConformerID(mol_h, ids, returnMinusOneIfError=False):
    best = []
    for id in ids:
        prop = AllChem.MMFFGetMoleculeProperties(mol_h)
        ff = AllChem.MMFFGetMoleculeForceField(mol_h, prop, confId=id)
        ff.Minimize()
        en = float(ff.CalcEnergy())
        econf = (en, id)
        best.append(econf)
    best.sort()
    print("Found "+str(len(best))+" conformations.")
    if len(best) == 0 and returnMinusOneIfError is True:
        return -1
    elif len(best) == 0:
        raise Exception("Error: No valid conformations found for the molecule. Try increasing the number of conformations.")
    # The best conformer is the first tuple and the ID of that conformer is the second value in the tuple
    best_id = int(best[0][1])
    # Return the best ID
    return best_id


# Takes savePath as optional argument for the directory to save in which you MUST NOT include a trailing slash
def SaveDBConformationsToFile(query, savePath="."):
    cocalcDB = MongoClient(database_IP, port=database_port)
    conformersCollection = cocalcDB.cocalc.conformers
    numResults = conformersCollection.count_documents(query)
    print("Found "+str(numResults)+" results.")
    conformationsWritten = 0
    results = conformersCollection.find(query)
    for result in results:
        filePath = savePath + "/" +str(result["precision"]) + "_conformations" + \
            "/" + result["repeatUnitName"] + "_" + result["initiatorName"] + "_" + \
            str(result["numRepeatUnits"]) + "-mer_(" + \
            str(result["precision"]) + "_conformations).mol"
        f = open(filePath, 'w+')
        f.write(result["bestConformerMolBlock"])
        f.close()
        conformationsWritten += 1
    print("Wrote "+str(conformationsWritten)+" conformations to files.")
    return 0


def DeduplicatePolymersWithDB(polymers):
    cocalcDB = MongoClient(database_IP, port=database_port)
    conformersCollection = cocalcDB.cocalc.conformers
    nonduplicated = []
    removedResults = 0
    retainedResults = 0
    totalResults = 0
    for polymer in polymers:
        totalResults += 1
        comparePolymer = {}
        comparePolymer["smiles"] = polymer["smiles"]
        comparePolymer["version"] = polymer["version"]
        comparePolymer["precision"] = polymer["precision"]
        numResults = conformersCollection.count_documents(comparePolymer)
        if numResults > 0:
            removedResults += 1
        else:
            retainedResults += 1
            nonduplicated.append(polymer)
    print("Processed "+str(totalResults)+" polymers with "+str(removedResults)+" removed and "+str(retainedResults)+" retained.")
    return nonduplicated


# HOW TO USE
# specify the precision and version of the data you want to get
# Xval and Yval can be any field in the database (check mongo.falcoperegrin.us for field names)
# given as a STRING. If you set showGraphs=True the function will display the graphs in your
# notebook (not recommended for large datasets). If you set a "filename" the function will
# save your graphs as a PDF.
def GetAllGraphs(precision=100, version=1, Xval="numRepeatUnits", Yval="LogPoverSA",
                 showGraphs=False, filename="test_graphs.pdf"):
    cocalcDB = MongoClient(database_IP, port=database_port)
    conformersCollection = cocalcDB.cocalc.conformers
    # Set up the figure counter
    figuresProcessed = 0
    # Get all the valid initiator names from the database
    all_initiators = conformersCollection.distinct('initiatorName',
                                                   filter={'precision':precision, 'version':version})
    # Using PDF functions...
    with PdfPages(filename) as pdf:
        for initiator in all_initiators:
            # For each valid initiator name in the database, get all the valid monomer names
            all_monomers = conformersCollection.distinct('repeatUnitName', filter={'precision':precision,
                                                                                'version':version,
                                                                                'initiatorName':initiator})
            for monomer in all_monomers:
                # For each valid monomer, pull the data and graph it.
                # Reset the arrays
                Xvalues = []
                Yvalues = []
                # First get the relevant data points from the database
                data_points = conformersCollection.find({'precision':precision, 'version':version,
                                                         'initiatorName':initiator, 'repeatUnitName':monomer},
                                                        sort=[('numRepeatUnits',pymongo.ASCENDING)])
                for data_point in data_points:
                    # For each data point, append the requested X and Y values to the array
                    Xvalues.append(data_point[Xval])
                    Yvalues.append(data_point[Yval])
                # Now graph the data
                fig_main = plt.figure()
                fig = fig_main.add_subplot()
                # Set X axis label to the user's request
                fig.set_xlabel(Xval+" ("+initiator+" "+monomer+")")
                # Set Y axis label to the user's request
                fig.set_ylabel(Yval)
                # Scatter the data
                fig.scatter(Xvalues, Yvalues, color='tab:red')
                fig.tick_params(axis='y')
                fig_main.tight_layout()
                # Save the figure to the PDF
                pdf.savefig(fig_main)
                # Close the figure so it doesn't contaminate the next one
                plt.close(fig_main)
                figuresProcessed += 1
                print("Output figure "+str(figuresProcessed))
    print("Saved "+str(figuresProcessed)+" figures to "+filename)
    return False


# Function to export spreadsheets from the database
# Accepts:
# filename - path and filename where you want to save the data
# query - a MongoDB query using the PyMongo syntax describing the data you want to retrieve
# fields - an array of the database field names you want to export
# in the same order you want them written to the spreadsheet
# fieldnames (optional) - the names that you want displayed for the fields in the spreadsheet
# (must by an array in the same order as fields)
# sort (optional) - an array of tuples with the names of fields you want sorted and the sort value 
# (1 for ascending, =1 for descending)
def MakeSpreadsheet(filename="test.xlsx", query={}, fields=[], fieldnames=[], sort={}):
    # Connect to the database
    DBconnect = MongoClient(database_IP, port=database_port)
    cocalcDB = DBconnect.cocalc
    conformersCollection = cocalcDB.conformers

    # start the timer
    executionStart = time.time()

    # Fetch the query, getting only the fields specified - if a sort was specified, also use that
    if sort:
        results = conformersCollection.find(query, fields).sort(sort)
    else:
        results = conformersCollection.find(query, fields)

    # Set up the file
    workbook = xlsxwriter.Workbook(filename)
    worksheet = workbook.add_worksheet()

    # output the header - write out all the names of the fields in the top row of the spreadsheet
    # if fieldnames is not empty, use that instead
    if len(fieldnames) != 0:
        for index, fieldname in enumerate(fieldnames):
            worksheet.write(0, index, fieldname)
    else:
        for index, field in enumerate(fields):
            worksheet.write(0, index, field)

    # start the row counter in the second row
    rows = 1

    # loop through the results
    for result in results:
        # for each result, loop through the fields
        for index, field in enumerate(fields):
            # for each field, write the contents to the corresponding row of the spreadsheet - use try/except in case some data points don't have all values
            try:
                worksheet.write(rows, index, result[field])
            except:
                print("Value not found for row "+str(rows)+" column "+str(index))
        # after all the fields have been written, increment the row counter
        rows += 1
    # we are finished so close the workbook
    workbook.close()
    print("Wrote out "+str(rows)+" rows with "+str(rows*len(fields))+" total cells")
    print("Completed output in "+str(time.time()-executionStart)+" seconds.")
    return None

# Included for reference
def MakePDFPrecisionComparison(precision=[10, 100, 200], version=1, Xval="numRepeatUnits", Yval="LogPoverSA",
                 showGraphs=False, filename="test_graphs_mv.pdf"):
    cocalcDB = MongoClient(database_IP, port=database_port)
    conformersCollection = cocalcDB.cocalc.conformers
    # Set up the figure counter
    figuresProcessed = 0
    # Get all valid initiators from the DB and concatenate them.
    # Get the initial values
    all_initiators = conformersCollection.distinct('initiatorName',
                                              filter={'precision': precision[0], 'version': version})
    for i in range(len(precision)-1):
        # for each i, set all initiators equal to the intersection of itself and the values for the next precision step
        all_initiators = intersection(all_initiators, conformersCollection.distinct('initiatorName',
                                                  filter={'precision': precision[i+1], 'version': version}))
    # Using PDF functions...
    with PdfPages(filename) as pdf:
        for initiator in all_initiators:
            # Get the first valid set of values
            all_monomers = conformersCollection.distinct('repeatUnitName',
                                                      filter={'precision': precision[0], 'version': version})
            # For each precision after the first one, fetch the data and get the intersection with the current data
            for i in range(len(precision)-1):
                all_monomers = intersection(all_monomers, conformersCollection.distinct('repeatUnitName', filter={'precision': precision[i+1],
                                                                                     'version': version,
                                                                                      'initiatorName': initiator}))
            ## We are going to construct a 3D array
            # Outer level - list of all monomers
            # Middle level - all precisions for a given monomer
            # inner level - data points for a given precision
            # we need to graph the inner two
            # Set up the master arrays

            all_colors = ['tab:red', 'tab:orange', 'y', 'tab:green', 'tab:blue', 'tab:purple', 'xkcd:light purple']
            for monomer in all_monomers:
                all_X_values = []
                all_Y_values = []
                # For each valid monomer, pull the data and graph it.
                # for each precision value given, get the data points
                for index, thisPrecision in enumerate(precision):
                    thisX = []
                    thisY = []
                    data_points = conformersCollection.find({'precision':thisPrecision, 'version':version,
                                                         'initiatorName':initiator, 'repeatUnitName':monomer},
                                                        sort=[('numRepeatUnits',pymongo.ASCENDING)])
                    # now extract the desired values and add them to the arrays
                    for data_point in data_points:
                        thisX.append(data_point[Xval])
                        thisY.append(data_point[Yval])
                    # Add the arrays to the master array
                    all_X_values.append(thisX)
                    all_Y_values.append(thisY)
                # Now graph the data
                fig_main = plt.figure()
                fig = fig_main.add_subplot()
                # Set X axis label to the user's request
                fig.set_xlabel(Xval+" ("+initiator+" "+monomer+")")
                # Set Y axis label to the user's request
                fig.set_ylabel(Yval)
                # Scatter the data
                # Loop over all the precision values for each point...
                for j in range(len(precision)):
                    # For each one, graph the X/Y using the next value from the color array
                    # Use modulo so if we run out of values we loop over and start again
                    fig.scatter(all_X_values[j], all_Y_values[j], color=all_colors[j % (len(all_colors))])
                    fig.plot(all_X_values[j], all_Y_values[j], color=all_colors[j % (len(all_colors))], label=(str(precision[j])))
                #fig.scatter(Xvalues1, Yvalues1, color='tab:red')
                #fig.plot(Xvalues1, Yvalues1, color='tab:red')
                #fig.scatter(Xvalues2, Yvalues2, color='tab:purple')
                #fig.plot(Xvalues2, Yvalues2, color='tab:purple')
                fig.tick_params(axis='y')
                fig_main.legend(loc="lower left")
                fig_main.tight_layout()
                # Save the figure to the PDF
                pdf.savefig(fig_main)
                # Close the figure so it doesn't contaminate the next one
                plt.close(fig_main)
                figuresProcessed += 1
                print("Output figure " + str(figuresProcessed))
                #print(len(Xvalues1)-len(Xvalues2))
    print("Saved "+str(figuresProcessed)+" figures to "+filename)
    return False

#upsert syntax

# db = client.acme
# c = db.pyshoot

# #testdata
# datatoinsert=[]
# #random data
# for i in range(0,20):
#     l=random.choice(string.ascii_letters[:5])
#     a=random.choice(string.ascii_letters[:5])

#     data={'monomer':l,
#         'triple':{
#             'version':i,
#             'smiles':random.choice(string.ascii_letters[:5]),
#             'perscion':100
#         },
#         'special':a+str(i),
#         'uniqueness':{
#                     'monomer':l,
#                     'smiles':a,
#                     'perscion':100},
#         'upsert': 'true'}
    
#     #enforces uniqueness across 3 different fields
#     c.update({'uniqueness':{
#             'monomer':l,
#             'smiles':a,
#             'perscion':100
#         }},data,upsert=True)

# #mongo "mongodb+srv://cluster0.kkvny.mongodb.net/acme" --username root

