from functools import reduce
from math import gcd


def getCoeffs(Coefs_and_monomers):
    #first, let us expand the provided coefs to actually give a weighted list.
    coeffs_list = []
    monomer_list = []
    #ommission of a coeficient implies 1 copy
    repeat_coef = 1
    for element in Coefs_and_monomers:
        try:
            repeat_coef = int(element)
        except:
            #if not string of integer, it should be considered a smiles (there will be an error from rdkit if not.)
            coeffs_list.append(repeat_coef)
            monomer_list.append(element)
            repeat_coef = 1

    return coeffs_list, monomer_list

def makePolymerBody_ratio(formula_list, n, verbo=False, mpn=1):
    coeffs, monomers = getCoeffs(formula_list)
    sum_coeffs = sum(coeffs)
    body_list = []
    real_coeffs = []
    roundup = True
    for i, coeff in enumerate(coeffs):
        unrounded_coeff = (coeff / sum_coeffs) * n * mpn
        #We need to make the case where 2 monomers have 0.5 split one to go up and the other down. However, "Integers" should not be changed.
        if unrounded_coeff % 0.5 == 0 and unrounded_coeff % 1 != 0:
            if roundup:  #the first one we see is rounded up
                unrounded_coeff += 0.1
                roundup = False
            else:  #the second one is rounded down and the flag is reset
                roundup = True
                unrounded_coeff -= 0.1
        #now we round and hopefully the length should always be correct
        real_coeff = round(unrounded_coeff) #we can't have float coeffs. Even though we are rounding everything should even out unless two monomers' coeffs initially end in .5
        body_list += [item for item in [monomers[i]] for _ in range(real_coeff)] #we repeat the monomer an appropriate number of times
        real_coeffs.append(real_coeff)

    #if n happens to be != to length of list, maybe allow user to modify by popping a random monomer from list or chosing a type of monomer to remove (and we remove a random one of that type.)
    #This isn't super important right now since they have a smiles they can also modify if they want exactly the right number.
    if n * mpn != len(body_list):
        print(f"WARNING: Due to rounding the length of the polymer is {len(body_list)} and the n specified was {n}.")
    if len(body_list) == 0:
        return None, "0:0"

    den = reduce(gcd, real_coeffs)
    ratio = ":".join(str(int(i / den)) for i in real_coeffs)

    #Now we have a list of monomers in the correct relative ammounts.
    #we just need to randomize the order.
    shuffle(body_list)
    smiles = "".join(body_list)
    if verbo:
        print(f"{n = }")
        print(f"Ratio of monomers used is {ratio}")
    return smiles, ratio
