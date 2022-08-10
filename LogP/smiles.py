monomer_dict = {

    #amino acids (use OH for termination group) 
    'L-ala': 'N[C@@H](C)C(=O)',
    'beta-ala': 'NCCC(=O)',
    'L-cystiene': 'N[C@@H](CS)C(=O)',
    'glycine': 'NCC(=O)',
    'L-isoleucine': 'N[C@@H](C(C)CC)C(=O)',
    'L-leucine': 'N[C@@H](CC(C)C)C(=O),
    'L-phenylala': 'N[C@@H](CC1=CC=CC=C1)C(=O),
    'L-phenylgly': 'N[C@@H](C1=CC=CC=C1)C(=O)',
    'L-valine': 'N[C@@H](C(C)C)C(=O)',
    
    #monomers for polyamides (use H for termination group)
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

#initiator dictionary
init_dict = {
    'Benzyl': '*Cc1ccccc1',
    'Benzyl_alcohol': '*Cc1ccccc1O', 
    'Benzoyl': 'c1ccc(cc1)C(=O)O',
    'Butyl': 'CCCC',
    'Hydroxyl': 'O', 
    'Hydrogen': '',
    'Methoxy': 'CO*',
    'Ethoxy': 'CCO*', 
    'Methyl': 'C', 
    'Vinyl': 'C=C'
}
