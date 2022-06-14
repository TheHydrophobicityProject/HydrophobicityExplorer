monomer_dict={

    #amino acids
    'ala': 'C[C@@H](C(=O)O)N',
    'ala-prot': 'C[C@@H](C(=O)[O-])[NH3+]',
    'beta-ala': 'NCCC(O)=O',
    'beta-ala-prot': '[NH3+]CCC([O-])=O',
    'L-cystiene': 'C([C@@H](C(=O)O)N)S',
    'L-cystiene-prot': 'C([C@@H](C(=O)[O-])[NH3+])S',
    'L-DOPA': 'N[C@H](C(O)=O)CC1=CC=C(O)C(O)=C1',
    'L-DOPA-prot': '[NH3+][C@H](C([O-])=O)CC1=CC=C(O)C(O)=C1',
    'glycine': 'C(C(=O)O)N',
    'glycine-prot': 'C(C(=O)[O-])[NH3+]',
    'L-isoleucine': 'N[C@@H](C(C)CC)C(O)=O',
    'L-isoleucine-prot': '[NH3+][C@@H](C(C)CC)C([O-])=O',
    'L-leucine': 'N[C@@H](CC(C)C)C(O)=O',
    'L-leucine-prot': '[NH3+][C@@H](CC(C)C)C([O-])=O',
    'L-phenylala': 'N[C@@H](CC1=CC=CC=C1)C(O)=O',
    'L-phenylala-prot': '[NH3+][C@@H](CC1=CC=CC=C1)C([O-])=O',
    'L-phenylala_F': 'N[C@@H](CC1=CC=C(F)C=C1)C(O)=O',
    'L-phenylala_F-prot': '[NH3+][C@@H](CC1=CC=C(F)C=C1)C([O-])=O',
    'L-phenylala_Cl': 'N[C@@H](CC1=CC=C(Cl)C=C1)C(O)=O',
    'L-phenylala_Cl-prot': '[NH3+][C@@H](CC1=CC=C(Cl)C=C1)C([O-])=O',
    'L-phenylala_N': 'N[C@@H](CC1=CC=C(N)C=C1)C(O)=O',
    'L-phenylala_N-prot': '[NH3+][C@@H](CC1=CC=C(N)C=C1)C([O-])=O',
    'L-phenylala_CN': 'N[C@@H](CC1=CC=C(C#N)C=C1)C(O)=O',
    'L-phenylala_CN-prot': '[NH3+][C@@H](CC1=CC=C(C#N)C=C1)C([O-])=O',
    'L-phenylgly': 'N[C@@H](C1=CC=CC=C1)C(O)=O',
    'L-phenylgly-prot': '[NH3+][C@@H](C1=CC=CC=C1)C([O-])=O',
    'L-proline': 'C1C[C@H](NC1)C(=O)O',
    'L-proline-prot': '[O-]C(=O)[C@H](CCC2)[NH2+]2',
    'L-threonine': 'C[C@H]([C@@H](C(=O)O)N)O',
    'L-threonine-prot': 'C[C@H]([C@@H](C(=O)[O-])[NH3+])O',
    'L-valine': 'CC(C)[C@@H](C(=O)O)N',
    'L-valine-prot': 'CC(C)[C@@H](C(=O)[O-])[NH3+]',

    #monomers for polyamides
    'Nylon6': 'CCCCCC(=O)N',

    #monomers and comonomers for polycarbonates
    'BPA_Carbonate':  'c1ccc(cc1)C(C)(C)c1ccc(cc1)OC(=O)O',

    #monomers and comonomers for polyesters
    'Butylene_Adip': 'C(=O)CCCCC(=O)OCCCCO',
    'Butylene_Succinate': 'C(=O)CCC(=O)OCCCCO', 
    'Butylene_Terephthalate': 'CCCCOC(=O)c1ccc(cc1)C(=O)O',
    'Ethylene_Terephthalate': 'CCOC(=O)c1ccc(cc1)C(=O)O', 
    '3HBV': 'C(CC)CC(=O)O', 
    '3HB': 'C(C)CC(=O)O',
    'Lactic_acid': 'C(C)C(=O)O',

    #monomers for polyethers
    'Ethylene_oxide': 'CCO',
    'Propylene_oxide': 'CC(C)O',

    #Vinyl monomers written to depict primary addition of alkene (i.e. substituent is on second carbon of alkene)
    'Butylacrylate': 'CC(C(=O)OCCCC)',
    'Dimethylacrylamide': 'CC(C(=O)N(C)(C))',
    'Ethylene': 'CC',
    'Hydroxyethylmethacrylate': 'CC(C(=O)OCCO)(C)',
    'Methylacrylate': 'CC(C(=O)OC)C',
    'Methylmethacrylate': 'CC(C(=O)OC)(C)',
    'Propylene': 'CC(C)',
    'Styrene': 'CC(c1ccccc1)',
    'Vinylalcohol': 'CC(O)'
}

#initiator and terminator dictionary (i.e. endgroups of polymer chain)
init_dict={
    'Benzyl': 'Cc1ccccc1',
    'Benzyl_alcohol': 'cc1ccccc1CO', 
    'Butyl': 'CCCC',
    'Hydroxyl': 'O', 
    'Hydrogen': '',
    'Methoxy': 'CO', 
    'Methyl': 'C', 
    'Vinyl': 'C=C'
}
