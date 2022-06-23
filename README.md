# Polymer Calculation Suite

This project enables users to perform several calulations on a limitless scope of polymers as long as the primary structure is known.

The script `MakePolymer.py` has a wide range of command-line options that allow fine control over polymer specification and output format. These can be obtained by running `python3 MakePolymer.py -h`. A large, but non-comprehensive list of examples will be covered in the next section.

## Usage and Examples

Note that the following examples all use the `-q` flag. This suppresses the default behavior, which asks the user for confirmation that the polymer they have specified has been interpreted correctly. This may be more important when using smiles inputs that do not come stock with this program. Addtional details about adding your own smiles to the included dictionaries can be found [here](README.md#modifying-the-smiles-dictionary). This is the default behavior to protect users from optimizing the geometry for incorrect polymers, which can be a time-consuming process for large values of `n`, but this can be dissabled for use in scripting or batch jobs.

### Specifying Polymer Components

There are dictionaries of monomers and terminal units in `smiles.py` The composition of a polymer containing units in these dictionaries can be spelled out in the following manner. The `-v` flag increases verbosity.

```bash
$ python3 MakePolymer.py -n 3 -m Styrene -v -q   # -n is number of monomers
polymer smiles is CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1) before any end groups
Polymer interpreted as: Hydrogen 3 * Styrene Hydrogen
This gives the following SMILES: CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)
Saving image to polymer.png by default.
requested calculations are None
```
Because verbosity was enabled an image of the polymer was saved with a default name. The name of the image can be specified with the `-d` flag.
The initiator and terminal groups default to Hydrogen if none are specified.

Here is another example with a more complex set of arguments:
```bash
$ python3 MakePolymer.py -n 2 -s 2 "CC(C(=O)OCCCC)" "CC(C)" -i Methoxy -t Benzyl -v -q
polymer smiles is CC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C)CC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C) before any end groups
polymer smiles is COCC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C)CC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C) after adding initiator smiles
polymer smiles is COCC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C)CC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C)Cc1ccccc1 after adding terminator smiles
Polymer interpreted as: Methoxy 2 * ['2', 'CC(C(=O)OCCCC)', 'CC(C)'] Benzyl
This gives the following SMILES: COCC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C)CC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C)Cc1ccccc1
Saving image to polymer.png by default.
requested calculations are None
```
The `-s` flag here allows us to define comonomers which are repeating sequences of smaller monomers in a specific order. The `-i` and `-t` flags are used to define initiators and terminators from either the dictionary or from SMILES. The `-s` flag can also be used to define monomers not in the dictionary with SMILES.

#### Modifying The SMILES Dictionary

Monomers should be added to the `monomer_dict` in `smiles.py` with the tail of the monomer at the left of the SMILES string and the head at the right. For example, propylene would be written 'CC(C)'. This allows easy construction of the polymer body by simply repeating this string `n` times.

Initiator and terminator groups should be added to `init_dict` and the atom to which the rest of the polymer should attatch must be denoted with `*`. Additionally, the SMILES must be written such that the `*` is the first or last character in the string. If the end-group is palindromic no asterisk is necessary. The existing dictionary has examples of each of these conditions.

### Reading a Polymer From A File

You will notice with the second example the run time is noticable since there are several conformations being compared to make the final mol object in rdkit. Additionally, this process is not perfectly reproducible. If desired, one can load a premade .mol or .pdb file instead of spelling out the polymer with the `-m` or `-s` flag. Polymers spelled out with the previously demonstrated methods can be converted to files as well with the `-f` flag. See the following section for details.

```bash
$ python3 MakePolymer.py -r pol.mol -c SA RG LogP -q
{'SA': 911.5262248851872, 'LogP': 14.510599999999974, 'RG': 7.430526236202889, 'N': None, 'smi': 'CCCCOC(=O)C(COC)CC(C)CC(C)CC(CC(C)CC(C)CC(CC(C)CC(C)CC(CC(C)CC(C)c1ccccc1CO)C(=O)OCCCC)C(=O)OCCCC)C(=O)OCCCC'}
```
The above example also shows how calculations are specified. Each calculation has a short string associated with it that can be use with the `-c` flag so only the desired calculations are performed. These can be found by using the `-h` flag. The data dictionary shows `'N' : None` because the smiles is not analyzed in any way in this configuration. However, this dictionary entry can be filled if the `-n` flag is used.

### Saving Polymer to File

The name or path of the file can be specified with the `-f` flag. Valid extentions are `.mol`, `.pdb` and `.xyz`. Be aware that `.xyz` files cannot be read back into this program.

```bash
$ python3 MakePolymer.py -n 4 -m Styrene -c MHP -v -f pol.mol -q
polymer smiles is CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1) before any end groups
Polymer interpreted as: Hydrogen 4 * Styrene Hydrogen
This gives the following SMILES: CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)
Saving image to polymer.png by default.
requested calculations are ['MHP']
{'SA': 605.1670849486483, 'LogP': 8.770700000000003, 'MHP': 0.014493022205172056, 'N': 4, 'smi': 'CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)'}
```
In the next plotting section it is shown that many polymers of different lengths can be generated with the `-p` flag. When the `-f` flag is specified as well, each of those molecules will be saved to its own file with a name based off the one specified as a CLI argument. The number of repeat units will be used as a suffix.

### Plotting Size-Dependent Calculations

The size-dependent plots of any calculations performed can be generated with the `-p` flag. The sizes plotted will range from 1 repeat unit to the number specified by the `-n` flag. Because the repeat unit needs to be well-defined, this plotting option is unavailable if the polymer is being read from a file. Since the `-v` flag is used, a grid image of all the generated molecules will be created as well.

The data can be exported with the `-e` flag. 

```bash
$ python3 MakePolymer.py -n 2 -m Styrene -c XMHP -p -e data.csv -v -q #XMHP requests that the MHP data be eXclusively returned.
polymer smiles is CC(c1ccccc1) before any end groups
Done generating SMILES with n = 1 now: CC(c1ccccc1)
Converting to mol now.
polymer smiles is CC(c1ccccc1)CC(c1ccccc1) before any end groups
Done generating SMILES with n = 2 now: CC(c1ccccc1)CC(c1ccccc1)
Converting to mol now.
Saving image to polymer.png by default.
requested calculations are ['XMHP']
Saved plot to Size-dependent-stats.png
{'MHP': [0.012376251236427096, 0.013573507027943717], 'N': [1, 2], 'smi': ['CC(c1ccccc1)', 'CC(c1ccccc1)CC(c1ccccc1)']}
#popup of plot appears.
Done exporting data to .csv file.
```
Nearly all of this output, including the plot popup and the polymer grid image, are excluded if the `-v` flag is excluded, but the image of the plot is still saved to a file. Because the calculation of MHP requires surface area and LogP values, they will always be included when MHP is specified in the list of desired calculations. When in this mode, if a name for an polymer file is specified, each molecule will be saved to a file based off the provided name. However, the number of mers in each molecule will be appended to the filename.

### Running Jobs With Config Files

Json files can be used instead of or in conjunction with any of the above command-line arguments.
Examples of valid json files are provided, but the most important aspect is the `runs` array:

```json
{
    "runs":
    [
        {
            "n": 4,
            "read": "pol.mol",
            "calculation": [
                "SA",
                "RG",
                "LogP"
            ]
        }
    ]
}
```

The list of dictionaries created by importing the paramaters from the json file is filled in with the default values for parameters like `"plot"` or `"export"`, then each dictonary of inputs is run in succession. The provided command-line arguments are overwritten by the corresponding argument in the json file.

## Dependencies

This project has been tested with `Python 3.10.4` and the following dependencies:

```python3
>>> rdkit.__version__
'2020.03.2'
>>> matplotlib.__version__
'3.5.2'
>>> argparse.__version__
'1.1'
>>> csv.__version__
'1.0'
>>> json.__version__
'2.0.9'
>>> PIL.__version__
'9.1.1'
```
