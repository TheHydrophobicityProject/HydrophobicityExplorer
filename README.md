# Polymer Calculation Suite

This project enables users to perform several calulations on a limitless scope of polymers as long as the primary structure is known.

The script `MakePolymer.py` has a wide range of command-line options that allow fine control over polymer specification and output format. These can be obtained by running `python3 MakePolymer.py -h`. A large, but non-comprehensive list of examples will be covered in the next section.

See [here](#installation) for installation instructions.

## Usage and Examples

Note that the following examples all use the `-q` flag. This suppresses the default behavior, which asks the user for confirmation that the polymer they have specified has been interpreted correctly. For example:

<img src="images/preview_example.png">

Is the confirmation image that would appear when [this example](#specifying-multiple-comonomers) is run without the `-q` flag.

This confirmation process may be more important when using smiles inputs that do not come stock with this program for the first time. Addtional details about adding your own smiles to the included dictionaries can be found [here](#modifying-the-smiles-dictionary). This is the default behavior to protect users from optimizing the geometry for incorrect polymers, which can be a time-consuming process for large values of `n`, but this can be dissabled for use in scripting or batch jobs.

Regardless of the provided value of `n`, the preview will only show one monomer or [block of comonomers](#specifying-multiple-comonomers) for simplicity.

This prompt will not be shown if both end groups are Hydrogen.

### Specifying Polymer Components

There are dictionaries of monomers and terminal units in `smiles.py` The composition of a polymer containing units in these dictionaries can be spelled out in the following manner. The `-v` flag increases verbosity.

```bash
$ python3 MakePolymer.py -n 3 -m Styrene -v -q   # -n specifies the number of monomers
polymer smiles is CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1) before any end groups
Polymer interpreted as: Hydrogen 3 * Styrene Hydrogen
This gives the following SMILES: CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)
Saving image to polymer.png by default.
requested calculations are None
```
Because verbosity was enabled an image of the polymer was saved with a default name. The name of the image can be specified with the `-d` flag.
The initiator and terminal groups default to Hydrogen if none are specified.

## Specifying Multiple Comonomers
Here is another example with a more complex set of arguments:
```bash
$ python3 MakePolymer.py -n 2 -b 2 Butylacrylate "CC(C)" -i Benzoyl -t Benzyl -v -q
polymer smiles is CC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C)CC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C) before any end groups
polymer smiles is c1ccc(cc1)C(=O)OCC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C)CC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C) after adding smiles string for initiator
polymer smiles is c1ccc(cc1)C(=O)OCC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C)CC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C)Cc1ccccc1 after adding terminator smiles
Polymer interpreted as: Benzoyl 2 * ['2', 'CC(C(=O)OCCCC)', 'CC(C)'] Benzyl
This gives the following SMILES: c1ccc(cc1)C(=O)OCC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C)CC(C(=O)OCCCC)CC(C(=O)OCCCC)CC(C)Cc1ccccc1
Saving image to polymer.png by default.
requested calculations are None
```
The `-b` flag defines the block of comonomers in a specific repeating pattern. The `-i` and `-t` flags are used to define initiators and terminators from either the dictionary or from SMILES. The `-b` flag can also be used to define monomers not in the dictionary with SMILES, but it also accepts dictonary keys. When a coefficient is provided in the list of arguments defined by the `-b` flag, this changes the number of monomers per unit defined by the `-n` flag. In the above example, each unit of n refers to 3 monomers. The number of monomers per n will be used for plots and image labels.

## Performing Calculations

The `-c` flag is followed by abbreviations for calulations that are desired. Availible options are:

LogP, SA (surface area), MV (Molecular Volume), MHP (LogP/SA; each of which will also be reported. Use XMHP to exclude those constituent calculations) and RG (radius of gyration).

When RG is selected, an exponential regression is performed. Polymer RG scales by n^(1/3) with the MMFF95 force field. This gives a sense of how reasonable the optimization steps were.

A run with Styrene to n=10 had the following regression:

<img src="images/RG_regression.png">

#### Modifying The SMILES Dictionary

Monomers should be added to the `monomer_dict` in `smiles.py` with the tail of the monomer at the left of the SMILES string and the head at the right. For example, propylene would be written `CC(C)`. This allows easy construction of the polymer body by simply repeating this string `n` times.

Initiator and terminator groups should be added to `init_dict` and the atom to which the rest of the polymer should attatch must be denoted with `*`. The use of `*` is inspired by polymergenome.org. Additionally, the SMILES must be written such that the `*` is the first or last character in the string. If the end-group is palindromic no asterisk is necessary. The existing dictionary has examples of each of these conditions.

### Reading a Polymer From A File

You will notice with the second example the run time is noticable since there are several conformations being compared to make the final mol object in rdkit. Additionally, this process is not perfectly reproducible. If desired, one can load a premade `.sdf`, `.mol` or `.pdb` file instead of spelling out the polymer with the `-m` or `-b` flag. Polymers spelled out with the previously demonstrated methods can be converted to files as well with the `-f` flag. See the following section for details.

```bash
$ python3 MakePolymer.py -r pol.mol -c SA RG LogP -q
{'SA': 911.5262248851872, 'LogP': 14.510599999999974, 'RG': 7.430526236202889, 'N': None, 'smi': 'CCCCOC(=O)C(COC)CC(C)CC(C)CC(CC(C)CC(C)CC(CC(C)CC(C)CC(CC(C)CC(C)c1ccccc1CO)C(=O)OCCCC)C(=O)OCCCC)C(=O)OCCCC'}
```
The data dictionary shows `'N' : None` because the smiles is not analyzed in any way in this configuration. However, this dictionary entry can be filled if the `-n` flag is used.

## Alternative Input Methods

### Custom Input
If the methods contained within this program are inadequet for the type of molecule desired, the accessory script `custom_input_to_mol_file.py` may be useful. It can read Smiles, Smarts or Inchi strings and produce a .mol file that can be read for calculations with the master script.

Use the following to show instructions for this script.
```bash
$ python3 custom_input_to_mol_file.py -h
```
### Random Composition
The accessory script `random_polymer_to_mol_file.py` can be used to interpret a ratio of monomers and develop a polymer that satisfy the user's desired monomer ratio. The monomers will be in a random order.

For example, the command `python3 random_polymer_to_mol_file.py -n 20 -m 2 Styrene Vinylalcohol -f rand.sdf` will generate a randomly ordered 20 unit-long polymer with a 2:1 ratio of Styrene to Vinylalcohol and save it to `rand_20.sdf` (The number of monomers is added to the filename automatically), which can be read by `MakePolymer.py`.

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

## Plotting Size-Dependent Calculations

The size-dependent plots of any calculations performed can be generated with the `-p` flag. The sizes plotted will range from 1 repeat unit to the number specified by the `-n` flag. Because the repeat unit needs to be well-defined, this plotting option is unavailable if the polymer is being read from a file. Since the `-v` flag is used, a grid image of all the generated molecules will be created as well.

The data can be exported to a `.csv` file with the `-e` flag. 

```bash
#XMHP requests that the MHP data be eXclusively returned instead of including the LogP and SA values as well.
$ python3 MakePolymer.py -n 2 -m Styrene -c XMHP -p -e data.csv -v -q
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

An example plot of Styrene LogP/SA with n <= 10:
<img src="images/plot_example.png">

This plot has a lot of noise at higher n. It is possible to [mitigate this](#reducing-noise-in-plots).

### Improving Data Quality in Plots

A plot generated with this tool may not always show a smooth curve, especially at higher n, where it is more difficult to optimize the polymer geometry. Calculations that depend on geometry like surface area or LogP/SA will be affected by poor geometry optimizations. This can be partially fixed by [changing the default settings](#changing-default-settings) to increase the number of conformations used for calculations or the maximum number of iterations a conformer is allowed to use before it is either accepted or discarded for not converging (i.e. the change in energy between optimization steps is still not yet small enough to be considered 'done'). Both of these in effect increase the population of conformers used to calculate an average of a requested property, which reduces noise on the graph.

The default settings are meant to find a balance between run time and reasonable results. These parameters are better suited to some monomers than others, so tweaking the parameters may be necessary not only to reduce noise on the plot, but to get the run to complete at all. The following is an example of how changing the default settings may benefit the plots produced by this program:

Styrene | MaxIters 1500 | numConfs 5
<img src="images/Styrene_default.png">

Styrene | MaxIters 2000 | numConfs 5
<img src="images/Styrene_2000i_5c.png">

Styrene | MaxIters 1500 | numConfs 10
<img src="images/Styrene_1500i_10c.png">

Styrene | MaxIters 2000 | numConfs 10
<img src="images/Styrene_2000i_10c.png">

In this case, increasing both parameters independently or together improved the quality of the graph. However, each increase in either parameter increased the run time.

## Changing Default Settings

Some settings are not accessible with command-line arguments. They can be changed in the file `settings.json`. Comments are not allowed in json files so each of these options are explained here.

```python
{
    "opt_numConfs":5, #The number of conformations you would like to generate. Increasing this greatly increases run time.
    "opt_numThreads":0, #The number of threads you would like to use for the optimization and conf generation. 0 means maximum possible.
    "opt_maxIters":1500, #The number of iterations used in optimization. Increase this if a job fails to optimize.
    "drawing_subImgSize_edge":250, #The side length of a subimage when saved. 
    "drawing_default":"polymer.png", #The name of an image saved by default (verbosity turned on with no image name specified.)
    "MV_gridSpacing":0.2, #Used for molar volume calculation.
    "MV_boxMargin" :2.0, #Used for molar volume calculation.
    "plot_dataPoint":"o", #Matplotlib argument to change plot point appearence.
    "plot_Filename":"Size-dependent-stats.png" #name of plot.
}
```
Some of the values are used by external functions and others by functions in this program. If it is desired to change the appearence of the plot points from blue circles to red pluses, change `"plot_dataPoint":"o"` to `"plot_dataPoint":"r+"`, as per the [matplotlib documentation](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

If the settings file cannot be found hardcoded defaults will be used instead.

## Running Jobs With Config Files

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
>>> PIL.__version__
'9.1.1'
```

# Installation

## General Steps
This project uses conda to manage dependencies. 

First:\
Install conda using the [official guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Then, if you are using Windows, follow the additional steps for that operating system. For other Linux or [windows subsystem for linux](https://docs.microsoft.com/en-us/windows/wsl/about), just open a terminal and skip to [the next section](#steps-for-all-users).

### Extra Steps For Windows

1. Install [git for Windows](https://git-scm.com/download/win)

2. Press the Windows Key and search for and open "Anaconda Prompt"

### Steps For All Users

1. Clone the repository\
`git clone https://github.com/scohenjanes5/MHP.git`

2. Set up the conda environment\
`conda create -c conda-forge -n mhp rdkit scipy matplotlib`

    Allow conda to install the dependencies

3. Activate the environment\
`conda activate mhp`

You can now run any of the scripts shown above with\
`python3 PATH/TO/SCRIPT -arg1 -arg2 ...`

