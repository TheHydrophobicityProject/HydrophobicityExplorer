import setuptools

setuptools.setup(
    name="hydrophobicity_explorer",
    version="0.1.4.0",
    url="https://github.com/TheHydrophobicityProject/HydrophobicityExplorer",
    author="Sander Cohen-Janes",
    author_email="sander.cohen-janes@yale.edu",
    description="Facilitates solubility calculations on a wide range of polymers.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=["rdkit", "scipy", "pandas", "matplotlib", "rich"],
    packages=setuptools.find_packages(),
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6'
    ],
    entry_points={
        'console_scripts': [ #this section allows us to use functions from the command line.
            # 'command = package.module:function',
            'makePol = hydrophobicity_explorer.MakePolymer:main',
            'customPol = hydrophobicity_explorer.custom_input_to_mol_file:main',
            'HXSettings = hydrophobicity_explorer.settings:main',
            'HXNB = hydrophobicity_explorer.nb:main',
            'HXLib = hydrophobicity_explorer.smiles:main'
        ],
    },
)
