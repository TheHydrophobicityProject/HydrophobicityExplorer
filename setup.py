import setuptools

setuptools.setup(
    name="mhp",
    version="0.1.3.4",
    url="https://github.com/TheHydrophobicityProject/HydrophobicityExplorer",
    author="Sander Cohen-Janes",
    author_email="sander.cohen-janes@yale.edu",
    description="Facilitates solubility calculations on a wide range of polymers.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=["rdkit", "scipy", "pandas", "matplotlib"],
    packages=setuptools.find_packages(),
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6'
    ],
    entry_points={
        'console_scripts': [ #this section allows us to use functions from the command line.
            # 'command = package.module:function',
            'makePol = mhp.MakePolymer:main',
            'customPol = mhp.custom_input_to_mol_file:main',
            'mhpSettings = mhp.settings:main',
            'mhpNB = mhp.nb:main',
            'mhpLib = mhp.smiles:main'
        ],
    },
)
