from hydrophobicity_explorer.settings import writeJson
import os
notebook = {
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {},
   "outputs": [],
   "source": [
    "from hydrophobicity_explorer.MakePolymer import main as makePol\n",
    "from hydrophobicity_explorer.custom_input_to_mol_file import main as customPol"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Available Keyword Arguments. More details on the github page.\\\n",
    "{'n': 0, 'initiator': 'Hydrogen', 'terminator': 'Hydrogen', 'single_monomer': 'None', 'comonomer_sequence': None, 'draw': None, 'verbose': False, 'calculation': None, 'save': None, 'read': None, 'plot': False, 'export': None, 'json': None, 'quiet': False, 'random': False}\\\n",
    "\n",
    "Only one of single_monomer or comonomer_sequence can be specified at a time\\\n",
    "initiator, terminator and single_monomer can be dict keys from smiles.py or smiles that are formatted correctly (see github page)\\\n",
    "draw is the filename to which a picture of the polymer can be saved\\\n",
    "verbose increases text output\\\n",
    "calculation is a list of any of the supported calculation types\\\n",
    "save is a path to a file to which to save the polymer as a one of a few supported molecule file formats\\\n",
    "read is the path to a file from which the program will read a polymer\\\n",
    "export is the name/path to a .csv file to which to export the data\\\n",
    "json is the name/path to a .json file which has parameters for several runs\\\n",
    "quiet suppesses the prompts to check proper connectivity of end groups\\\n",
    "random allows comonomer_sequence to be interpreted as a desired ratio of comonomers randomly ordered rather than a repeating block.\\"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Polymer interpreted as: Hydrogen 3 * Styrene Hydrogen\n",
      "This gives the following SMILES: CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)\n",
      "failed to remove tmp file.\n",
      "Saving image to polymer.png by default.\n",
      "requested calculations are ['XMHP', 'RG']\n",
      "         RG   LogP/SA  N                                   smi\n",
      "0  3.539339  0.023622  3  CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)\n"
     ]
    }
   ],
   "source": [
    "makePol(n=3, single_monomer=\"Styrene\", calculation=[\"XMHP\", \"RG\"], verbose=True)"
   ]
  }
 ],
    "metadata": {
        "kernelspec": {
            "display_name": "Python 3.10.6 ('HX-env')",
            "language": "python",
            "name": "python3"
        },
        "language_info": {
            "name": "python",
            "version": "3.10.6"
        },
        "orig_nbformat": 4
    },
    "nbformat": 4,
    "nbformat_minor": 2
}

name = "hydrophobicity_explorer.ipynb"

def main():    
    if os.path.exists(name):
        inp = input(f"{name} exists. Should it be overwritten? [Y/n]: ")
        if inp.lower() != "y" and inp != "":
            print("Please rename the existing notebook so it is not overwritten.")
            quit()
    
    writeJson(notebook, name)
    print(f"Done creating {name}")

if __name__ == "__main__":
    main()