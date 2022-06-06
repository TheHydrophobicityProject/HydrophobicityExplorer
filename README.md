# Polymer Calculation Suite

This project enables users to perform several calulations on a limitless scope of polymers as long as the primary structure is known.

The script `MakePolymer.py` has a wide range of command-line options that allow fine control over polymer specification and output format. These can be obtained by running `python3 MakePolymer.py -h`. A large, but non-comprehensive list of examples will be covered in the next section.

## Usage and Examples

There is a dictionary of monomers and terminal units builtin `MakePolymer.py` The composition of a polymer containing units in these dictionaries can be spelled out in the following manner. The -v flag increases verbocity.

```python
python3 MakePolymer.py -n 3 -m Styrene -v
Polymer interpreted as: Hydrogen 3 * Styrene Hydrogen
This gives the following SMILES: CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)
Saving image to polymer.png by default.

requested calculations are None
```
Because verbocity was enabled an image of the polymer was saved with a default name. The name of the image can be specified with the -d flag.
The initiator and terminal groups default to Hydrogen if none are specified.

Here is another example with a more complex set of arguments:
```python
python3 MakePolymer.py -n 4 -s 2 "CC(C(=O)OCCCC)" "CC(C)" -i Methoxy -t Benzyl -v
Polymer interpreted as: Methoxy 4 * ['2', 'CC(C(=O)OCCCC)', 'CC(C)'] Benzyl
This gives the following SMILES: COCC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C)CC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C)CC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C)CC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C)c1ccccc1CO
Saving image to polymer.png by default.
requested calculations are None
```
The `-s` flag here allows us to define a "super-monomer" which is a repeating sequence of smaller monomers in a specific order. The `-i` and `-t` flags are used to define initiators and terminators from either the dictionary or from SMILES. The `-s` flag can also be used to define monomers not in the dictionary with SMILES.

You will notice with the second example the run time is noticable since there are several conformations being compared to make the final mol object in rdkit. Additionally, this process is not perfectly reproducible. If desired, one can load a premade .mol or .pdb file instead of spelling out the polymer with the `-m` or `-s` flag. Polymers spelled out with the previously demonstrated methods can be converted to files as well with the `-f` flag. See the builtin help for details on that.

```python
python3 MakePolymer.py -r pol.mol -c SA RG LogP
{'SA': 911.5262248851872, 'LogP': 14.510599999999974, 'RG': 7.430526236202889, 'N': None, 'smi': 'CCCCOC(=O)C(COC)CC(C)CC(C)CC(CC(C)CC(C)CC(CC(C)CC(C)CC(CC(C)CC(C)c1ccccc1CO)C(=O)OCCCC)C(=O)OCCCC)C(=O)OCCCC'}
```
The above example also shows how calculations are specified. Each calculation has a short string associated with it that can be use with the `-c` flag so only the desired calculations are performed. These can be found by using the `-h` flag. The data dictionary shows `'N' : None` because the smiles is not analyzed in any way in this configuration. However, this dictionary entry can be filled if the `-n` flag is used.

The size-dependent plots of any calculations performed can be generated with the `-p` flag and the data can be exported with the `-e` flag. 

```python
python3 MakePolymer.py -n 4 -m Styrene -c MHP -p -e data.csv -v
Done generating SMILES with n = 1 now: CC(c1ccccc1)
Converting to mol now.
Done generating SMILES with n = 2 now: CC(c1ccccc1)CC(c1ccccc1)
Converting to mol now.
Done generating SMILES with n = 3 now: CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)
Converting to mol now.
Done generating SMILES with n = 4 now: CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)
Converting to mol now.
Saving image to polymer.png by default.
requested calculations are ['MHP']
Saved plot to Size-dependent-stats.png
<popup of plot appears>
done exporting data to .csv file.
{'SA': [181.71900012667044, 325.84799130354446, 466.9956239244017, 607.5159276234424], 'LogP': [2.2490000000000006, 4.422900000000004, 6.596800000000006, 8.770700000000003], 'MHP': [0.012376251236427096, 0.013573507027943717, 0.014126042433896363, 0.01443698774172117], 'N': [1, 2, 3, 4], 'smi': ['CC(c1ccccc1)', 'CC(c1ccccc1)CC(c1ccccc1)', 'CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)', 'CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)']}
```
Nearly all of this output, including the plot popup and the polymer grid image, are excluded if the `-v` flag is excluded, but the image of the plot is still saved to a file. Because the calculation of MHP requires surface area and LogP values, they will always be included when MHP is specified in the list of desired calculations. When in this mode, if a name for an polymer file is specified, each molecule will be saved to a file based off the provided name. However, the number of mers in each molecule will be appended to the filename.

## Dependencies

This project has been tested with the following dependencies:

```python
>>> rdkit.__version__
'2020.03.2'
>>> matplotlib.__version__
'3.5.2'
>>> argparse.__version__
'1.1'
>>> csv.__version__
'1.0'
```
