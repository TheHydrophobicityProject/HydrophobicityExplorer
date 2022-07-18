import setuptools

setuptools.setup(
    name="mhp",
    version="0.1.0",
    url="https://github.com/scohenjanes5/MHP",
    author="Sander Cohen-Janes",
    author_email="scohenjanes@brandeis.edu",
    description="Facilitates solubility calculations on a wide range of polymers.",
    long_description=open('README.md').read(),
    packages=setuptools.find_packages(),
    install_requires=["rdkit", "argparse", "random", "csv"],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
    ],
)