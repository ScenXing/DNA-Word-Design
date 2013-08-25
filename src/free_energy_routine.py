# Author: Nguyen Truong Duy
# Email: truongduy134@gmail.com

# Self-written modules
import m_poly

# Global variables
TYPE_ARR = ['A', 'C', 'G', 'T']
NUM_TYPE = 4

#############################################################################
# compute_free_energy(str1, pairwiseenergy):
#
# The way we compute the approximated free energy of a DNA string X is specified
#    Section 4 of the paper.
# freeenergy(X) = sum pairwiseenergy(X(i), X(i + 1)) for i = 0, ..., L - 2 
def compute_free_energy(str1, pairwiseenergy):
    """
    Returns the approximated free energy of the input DNA string.

    Inputs:
    + str1: a Python string that contains 'A', 'C', 'G', 'T' only
    + pairwiseenergy: a function that takes in two inputs X and Y where X, Y in {'A', 'C', 'G', 'T'} and returns a non-negative integer indicating the free energy between the ordered pair (X, Y).

    Output:
     + freeenergy: an integer indicating the approximated free energy of the given DNA word, computed by the following formula:
    freeenergy = sum of pairwiseenergy(str1[i], str1[i + 1]) for i = 0, ..., len(str1) - 2
    """

    freeenergy = 0
    len1 = len(str1)
    for ind in xrange(len1 - 1):
        freeenergy += pairwiseenergy(str1[ind], str1[ind + 1])
    return freeenergy

##############################################################################
def find_extreme_pairwise_energy(pairwiseenergy, extremefunc):
    """
    Finds an extreme value (min / max) among all 16 possible entries returned from the given pairwise free energy function.

    Inputs:
    + pairwiseenergy: a function that takes in two inputs X and Y where X, Y in {'A', 'C', 'G', 'T'} and returns a non-negative integer indicating the free energy between the ordered pair (X, Y).
    + extremefunc: a binary function over the set of integers (e.g. min or max)
 
    Output:
    + extremeval: an integer indicating the extreme value of the pairwise free energy function.    
    """

    extremeval = pairwiseenergy('A', 'A')
    for typeone in xrange(NUM_TYPE):
        for typetwo in xrange(NUM_TYPE):
            extremeval = extremefunc(extremeval, pairwiseenergy(TYPE_ARR[typeone], TYPE_ARR[typetwo]))

    return extremeval

################################################################################
def create_pairwise_energy_func(pairwiseenergymat, indarr = None):
    """
    Create a pairwise free energy function from the pairwise free energy matrix.

    Inputs:
    + indarr: a Python list which can be viewed as a permutation of 4 characters 'A', 'C', 'G', 'T'. The default value is ['A', 'C', 'G', 'T']
    + pairwiseenergymat: a 4 x 4 matrix of non-negative integers where rows and columns are indexed by indarr in that order, presenting the pairwise free energy values.

    Output:
    + pairwisefunc: a function that takes in two inputs X and Y where X, Y in {'A', 'C', 'G', 'T'} and returns a non-negative integer indicating the free energy between the ordered pair (X, Y).
    """

    indarr = indarr or ['A', 'C', 'G', 'T']
    maptypetoind = {}
    for ind in xrange(NUM_TYPE):
        maptypetoind[indarr[ind]] = ind
    
    def pairwisefunc(typeone, typetwo):
        return pairwiseenergymat[maptypetoind[typeone]][maptypetoind[typetwo]]
    return pairwisefunc

#################################################################################
# compute_min_free_energy(l, pairwiseenergy):
#
# Methodology:
#    + We use dynamic programming to compute the value.
#     + Define:
#        f(l, c) = the minimum free energy among all possible values of all
#        DNA strings of length l and the first character is c
#    + Then the value we want to compute is:
#        min{f(l, 'A'), f(l, 'C'), f(l, 'G'), f(l, 'T')}
#
#    + Recurrence relation:
#        f(l, c) = min{pairwiseenergy(c, d) + f(l - 1, d)} where d in {'A', 'C', 'G', 'T'}
#    + Base cases:
#        f(1, c) = 0 for all c
#        f(0, c) = 0 for all c
#
#    + We implment compute_min_free_energy_dp with recursion and memoization to compute f
#
# Time complexity: O(l)
#
# Note on the implementation:
#    + We use integers to represent the type of nucleotide, i.e. i represents TYPE_ARR[i]
#      (TYPE_ARR is a global variable) for i = 0, 1, 2, 3
def compute_min_free_energy(l, pairwiseenergy):
    """
    Compute the minimum free energy among all possible values of all DNA strings of a given length with respect to the given pairwise energy function.

    Inputs:
    + l: an integer indicating the length of DNA words we consider
    + pairwiseenergy: a function that takes in two inputs X and Y where X, Y in {'A', 'C', 'G', 'T'} and returns a non-negative integer indicating the free energy between the ordered pair (X, Y).

    Output:
    + minenergy: an integer indicating the minimum free energy of all DNA words of length l
    """

    memotable = [[-1] * l for i in xrange(NUM_TYPE)]

    minenergy = compute_min_free_energy_dp(l, 0, pairwiseenergy, memotable)
    for typeid in xrange(1, NUM_TYPE):
        minenergy = min(minenergy, compute_min_free_energy_dp(l, typeid, pairwiseenergy, memotable))

    return minenergy

#
# compute_min_free_energy_dp(l, typeid, pairwiseenergy, memotable):
#
#    Refer to the documentation of compute_min_free_energy
#
def compute_min_free_energy_dp(l, typeid, pairwiseenergy, memotable):
    """

    Inputs:
    + l: an integer indicating the length of a DNA word
    +
    + pairwiseenergy: a function that takes in two inputs X and Y where X, Y in {'A', 'C', 'G', 'T'} and returns a non-negative integer indicating the free energy between the ordered pair (X, Y).
    + memotable: 
    Output:
    """
    if l <= 1:
        return 0
    if memotable[typeid][l - 1] >= 0:
        return memotable[typeid][l - 1]

    
    minenergy = -1
    for nexttypeid in xrange(NUM_TYPE):
        energy = pairwiseenergy(TYPE_ARR[typeid], TYPE_ARR[nexttypeid]) + compute_min_free_energy_dp(l - 1, nexttypeid, pairwiseenergy, memotable)
        if minenergy < 0 or energy < minenergy:
            minenergy = energy

    memotable[typeid][l - 1] = minenergy
    return minenergy

###############################################################################################
# compute_exist_free_energy_dna_mat(maxlen, pairwiseenergy, maxenergy = -1):
#    + Computes and returns a binary matrix M of dimension NUM_TYPE x maxlen x maxenergy 
#      (where NUM_TYPE = 4) such that
#        M[c][L][E] = 1 means there exists a DNA string of length (L + 1) that 
#            starts with the character c and has free energy E
#            with respect to the given pairwiseenergy function.
#        M[c][L][E] = 0 otherwise 
#
# Input specification:
#    + maxlen: an integer indicating the maximum length of a DNA string
#    + pairwiseenergy: a function that takes in two inputs X and Y where
#        X, Y in {'A', 'C', 'G', 'T'} and returns the free energy
#        between the ordered pair (X, Y), which is an integer value.
#    + maxenergy: an integer indicating the maximum energy we are considering.
#        maxenergy = -1 if it is not specified.
#
# This is how we determine maxenergy (Note that our resulted matrix is of size maxlen x maxenergy):
#    + Let W be maximum value among 16 values of pairwiseenergy
#    + Let maxpossibleenergy = W * (maxlen - 1). maxpossibleenergy is the maximum free energy possible
#        for all DNA strings of length maxlen
#    + If maxenergy > 0 and maxenergy < maxpossibleenergy, then maxenergy = maxenergy
#    + Otherwise, maxenergy = maxpossibleenergy
#
# Methodology to compute the matrix:
#    + We use dynamic programming to solve this problem.
#    + Define the function
#        f(c, l, E) = 1 if there exists a DNA string of length l that 
#            starts with the character c and has free energy E
#        f(c, l, E) = 0 otherwise
#    + We see that M[c][L][E] = f(c, L + 1, E)
#    + Recurrence relation:
#        f(c, l, E) = OR{f(d, l - 1, E - pairwiseenergy(c, d)} for d in {'A', 'C', 'G', 'T'}
#    + Base cases:
#        f(c, l, 0) = 1 if l <= 1 for all c
#        f(c, l, E) = 0 if l <= 1 and E > 0
#        f(c, l, E) = 0 if E < 0
#
#    + We implment exist_free_energy_dna_dp with recursion and memoization to compute f
#
# Time complexity: O(L^2) (Note that maxenergy = O(L)
#
# Note on the implementation:
#    + We use integers to represent the type of nucleotide, i.e. i represents TYPE_ARR[i]
#      (TYPE_ARR is a global variable) for i = 0, 1, 2, 3
#
# Important note: M[c][L][E] = f(c, L + 1, E)
#
def compute_exist_free_energy_dna_mat(maxlen, pairwiseenergy, maxenergy = -1):
    """
    Computes and returns a binary matrix M of dimension NUM_TYPE x maxlen x maxenergy (where NUM_TYPE = 4) such that
    """

    # Find maximum possible energy which is equal to
    #    (maxlen - 1) * maxpairwiseenergy
    # where maxpairwiseenergy = max(pairwiseenergy(X, Y)) (X, Y can be 'A', 'C', 'G', 'T')
    maxpairwiseenergy = find_extreme_pairwise_energy(pairwiseenergy, max)
    maxpossibleenergy = (maxlen - 1) * maxpairwiseenergy
    if maxenergy < 0 or maxenergy > maxpossibleenergy:
        maxenergy = maxpossibleenergy

    # Declare and initialize existence matrix
    existmat = [[[-1] * (maxenergy + 1) for l in xrange(maxlen)] for typeid in xrange(NUM_TYPE)]

    # Repeatedly call the recursive routine (with memoization) to compute all entries of 
    # the existence matrix
    for typeid in xrange(NUM_TYPE):
        for l in xrange(maxlen):
            for energy in xrange(maxenergy + 1):
                if existmat[typeid][l][energy] < 0:
                    existmat[typeid][l][energy] = exist_free_energy_dna_dp(l + 1, typeid, energy, pairwiseenergy, existmat)

    return existmat

#
# exist_free_energy_dna_dp(l, typeid, freeenergy, pairwiseenergy, memotable):
#
#   Refer to the documentation of compute_exist_free_energy_dna_mat
#
def exist_free_energy_dna_dp(l, typeid, freeenergy, pairwiseenergy, memotable):
    """
    """
    if freeenergy < 0:
        return 0
    if l <= 1:
        if freeenergy == 0:
            return 1
        return 0

    if memotable[typeid][l - 1][freeenergy] >= 0:
        return memotable[typeid][l - 1][freeenergy]

    existflg = 0
    for nexttypeid in xrange(NUM_TYPE):
        existflg = existflg or exist_free_energy_dna_dp(l - 1, nexttypeid, freeenergy - pairwiseenergy(TYPE_ARR[typeid], TYPE_ARR[nexttypeid]), pairwiseenergy, memotable)
        if existflg:
            break

    memotable[typeid][l - 1][freeenergy] = existflg
    return existflg

##############################################################################################
# construct_bounded_energy_dna_list(triplelist, pairwiseenergy):
#     + Constructs and returns the list of DNA strings with length and bounded energy 
#      satisfied the conditions specified in triplelist. 
#
# Input specification:
#    + triplelist: a list of triples (l, A, B). Each triple (l, A, B) contains information
#        about the DNA we want to construct:
#            - l: the length of the DNA string
#            - A: the lower bound of the free energy of the DNA string
#            - B: the upper bound of the free energy of the DNA string
#    + pairwiseenergy: a function that takes in two inputs X and Y where
#        X, Y in {'A', 'C', 'G', 'T'} and returns the free energy
#        between the ordered pair (X, Y), which is an integer value.
#
# Detailed Flow:
#    + Let L be the maximum length of all DNA strings we want to construct
#      Let E be the maximum possible free energy 
#    + Construct the existence matrix M of dimension NUM_TYPE x L x E (using the function
#        compute_exist_free_energy_dna_mat) (NUM_TYPE = 4 typically) such that:
#            M[c][l][e] = 1 means there exists a DNA string of length (l + 1) that 
#                starts with the character c and has free energy e
#                with respect to the given pairwiseenergy function.
#            M[c][l][e] = 0 otherwise
#    + For each triple (l, A, B), we search if there exists some energy e BETWEEN A and B
#        such that M[c][l][e] = 1 for some type c
#      If such e exists, we call the function construct_free_energy_dna to construct
#        the DNA string of length l, starting with character c and has free energy e
#        which is bounded between A and B
#
# Time complexity: O(L^2 + n * L) since (E = O(L)) where
#    - O(L^2) is due to constructing the existence matrix
#    - O(n * L) is due to constructing n bounded free energy DNA strings of length at most L
def construct_bounded_energy_dna_list(triplelist, pairwiseenergy):
    """
    """
    numdna = len(triplelist)
    if numdna == 0:
        return []

    maxlen = max([triplelist[i][0] for i in xrange(numdna)])
    maxenergy = max([triplelist[i][2] for i in xrange(numdna)])
    maxenergy = max(maxenergy, max([triplelist[i][1] for i in xrange(numdna)]))

    existmat = compute_exist_free_energy_dna_mat(maxlen, pairwiseenergy, maxenergy)
    maxpossibleenergy = len(existmat[0][0]) - 1

    resultstrlist = []

    # Assumption: l >= 1, 0 <= energy <= maxpossibleenergy
    def exist_energy_dna(l, energy):
        for typeid in xrange(NUM_TYPE):
            if existmat[typeid][l - 1][energy]:
                return True
        return False

    for ind in xrange(numdna):
        # Preprocess / Standardize input data
        minenergy = triplelist[ind][1]
        if minenergy < 0:
            minenergy = 0
        maxenergy = triplelist[ind][2]
        if maxenergy > maxpossibleenergy:
            maxenergy = maxpossibleenergy
        if minenergy > maxenergy:
            temp = minenergy
            minenergy = maxenergy
            maxenergy = temp
        l = triplelist[ind][0]

        if l < 1:
            resultstrlist.append('')
            continue

        found = False
        for energy in xrange(minenergy, maxenergy + 1):
            if exist_energy_dna(l, energy):
                found = True
                resultstrlist.append(construct_free_energy_dna(l, energy, pairwiseenergy, existmat))
                break

        if not found:
            resultstrlist.append('')

    return resultstrlist        
            
    
##############################################################################################
# construct_free_energy_dna(l, energy, pairwiseenergy, existmat):
#    + Constructs and returns a DNA string of length l and having free energy as specified
#        in the second parameter with respect to the input pairwiseenergy function
#
# Input specfication:
#    + l: an integer indicating the length of the DNA string to be constructed
#    + energy: the free energy of the DNA string to be constructed
#    + pairwiseenergy: a function that takes in two inputs X and Y where
#        X, Y in {'A', 'C', 'G', 'T'} and returns the free energy
#        between the ordered pair (X, Y), which is an integer value.
#    + existmat: a 3D array of dimension NUM_TYPE x L x E (using the function
#        compute_exist_free_energy_dna_mat) (NUM_TYPE = 4 typically) such that:
#            - L >= l
#            - E > energy
#            - M[c][l][e] = 1 means there exists a DNA string of length (l + 1) that 
#                starts with the character c and has free energy e
#                with respect to the given pairwiseenergy function.
#            - M[c][l][e] = 0 otherwise
#
# Detail:
#    + Clearly, if l < 1, then we return the empty string
#    + For the first character c1 of the DNA, we choose any c1 among {'A', 'C', 'G', 'T'} such that
#        existmat[c1][l - 1][energy] = 1
#    + For the second character c2 of the DNA, our problem now becomes finding a DNA string of length
#        l - 1, starting with c2 and have energy equals to energy - pairwiseenergy(c1, c2). This 
#        means we choose any c2 among {'A', 'C', 'G', 'T'} such that
#            existmat[c2][l - 2][energy - pairwiseenergy(c1, c2)] = 1
#    + Continue in this manner until we find all n characters
#
# Time complexity: O(l)
def construct_free_energy_dna(l, energy, pairwiseenergy, existmat):
    if l < 1:
        return ''

    chararr = []
    curlen = 0
    remainlen = l
    while remainlen > 0:
        exist = False
        
        for typeid in xrange(NUM_TYPE):
            remainenergy = energy
            if curlen > 0:
                remainenergy -= pairwiseenergy(chararr[curlen - 1], TYPE_ARR[typeid])

            if existmat[typeid][remainlen - 1][remainenergy]:
                chararr.append(TYPE_ARR[typeid])
                energy = remainenergy
                curlen += 1
                remainlen -= 1
                exist = True
                break

        if not exist:
            return ''

    return ''.join(chararr)

############################################################################################3
#
def build_free_energy_poly(length, pairwiseenergy):
    polypairlist = []

    if length == 0:
        return polypairlist

    if length == 1:
        dictx = {}
        dicty = {}
        for typea in TYPE_ARR:
            for typeb in TYPE_ARR:
                if typea == typeb:
                    dictx[(length, typea, typeb)] = m_poly.MyPoly([1])
                else:
                    dictx[(length, typea, typeb)] = m_poly.MyPoly()
        for typea in TYPE_ARR:
            for typeb in TYPE_ARR:
                for d1 in TYPE_ARR:
                    for d2 in TYPE_ARR:
                        for d3 in TYPE_ARR:
                            dicty[(length, typea, typeb, d1, d2, d3)] = None
        polypairlist.append((dictx, dicty))

        return polypairlist

    # Recurse
    polypairlist = build_free_energy_poly(length >> 1, pairwiseenergy)

    if (length & 1) == 0:
        # length is even
        dicty = {}
        for typea in TYPE_ARR:
            for typeb in TYPE_ARR:
                for d1 in TYPE_ARR:
                    for d2 in TYPE_ARR:
                        temppoly = polypairlist[-1][0][(length >> 1, typea, d1)] * polypairlist[-1][0][(length >> 1, d2, typeb)]
                        temppoly = temppoly.multiplyBySingleton(pairwiseenergy(d1, d2), 1)
                        dicty[(length, typea, typeb, d1, d2)] = temppoly

        
        dictx = {}
        for typea in TYPE_ARR:
            for typeb in TYPE_ARR:
                poly = m_poly.MyPoly()
                for d1 in TYPE_ARR:
                    for d2 in TYPE_ARR:
                        poly += dicty[(length, typea, typeb, d1, d2)]
                dictx[(length, typea, typeb)] = poly

        polypairlist.append((dictx, dicty))
    else:
        # length is odd
        dicty = {}
        for typea in TYPE_ARR:
            for typeb in TYPE_ARR:
                for d1 in TYPE_ARR:
                    for d2 in TYPE_ARR:
                        for d3 in TYPE_ARR:
                            temppoly = polypairlist[-1][0][(length >> 1, typea, d1)] * polypairlist[-1][0][(length >> 1, d2, typeb)]
                            temppoly = temppoly.multiplyBySingleton(pairwiseenergy(d1, d3) + pairwiseenergy(d3, d2), 1)
                            dicty[(length, typea, typeb, d1, d2, d3)] = temppoly

        dictx = {}
        for typea in TYPE_ARR:
            for typeb in TYPE_ARR:
                poly = m_poly.MyPoly()
                for d1 in TYPE_ARR:
                    for d2 in TYPE_ARR:
                        for d3 in TYPE_ARR:
                            poly += dicty[(length, typea, typeb, d1, d2, d3)]
                dictx[(length, typea, typeb)] = poly
        polypairlist.append((dictx, dicty))    

    return polypairlist

####################################################################################
#
def extract_free_energy_dna_poly(length, energy, pairwiseenergy, polypairlist):
    lastind = len(polypairlist) - 1
    for typea in TYPE_ARR:
        for typeb in TYPE_ARR:
            if has_positive_energy_coeff(polypairlist[-1][0][(length, typea, typeb)], energy):
                chararr = recursive_extract(length, energy, typea, typeb, pairwiseenergy, polypairlist, lastind)
                return ''.join(chararr)

    return ''    # Default: Return empty string

################################################################################################
#
def recursive_extract(length, energy, typea, typeb, pairwiseenergy, polypairlist, lengthpolyind):
    # Base cases:
    if length == 2:
        return [typea, typeb]
    if length == 1:
        return [typea]

    halflength = (length >> 1)
    
    paramtuple = None
    #########################################################
    # CHECK THIS PART AGAIN
    if (length & 1) == 0:
        # length is even
        # Find triple (z, d1, d2)
        paramtuple = findd1d2(halflength, energy, typea, typeb, pairwiseenergy, polypairlist[lengthpolyind - 1][0])
    else:
        # length is odd
        # Find 4-tuple (z, d1, d2, d3)
        paramtuple = findd1d2d3(halflength, energy, typea, typeb, pairwiseenergy, polypairlist[lengthpolyind - 1][0])
    ########################################################3
    if paramtuple is None:
        return []

    chararr = recursive_extract(halflength, paramtuple[0], typea, paramtuple[1], pairwiseenergy, polypairlist, lengthpolyind - 1)

    if chararr == []:
        return []
    
    subenergy = 0
    if (length & 1) == 0:
        # length is even
        subenergy = energy - paramtuple[0] - pairwiseenergy(paramtuple[1], paramtuple[2]) 
    else:
        # length is odd
        subenergy = energy - paramtuple[0] - pairwiseenergy(paramtuple[1], paramtuple[3]) - pairwiseenergy(paramtuple[3], paramtuple[2])
        chararr.append(paramtuple[3])
    chararr.extend(recursive_extract(halflength, subenergy, paramtuple[2], typeb, pairwiseenergy, polypairlist, lengthpolyind - 1))

    return chararr

####################################################################################
#
def has_positive_energy_coeff(poly, energy):
    if energy < 0:
        return False
    if poly.degree < energy:
        return False

    return poly.coeff[energy] > 0

#######################################################################################
#
def findd1d2(length, energy, typea, typeb, pairwiseenergy, polyset):
    for d1 in TYPE_ARR:
        polyd1 =  polyset[(length, typea, d1)]
        maxenergy = polyd1.degree
        for subenergy in xrange(maxenergy + 1):
            if polyd1.coeff[subenergy] < 1:
                continue
            if energy - subenergy < 0:
                break

            for d2 in TYPE_ARR:
                leftenergy = energy - subenergy - pairwiseenergy(d1, d2)
                if has_positive_energy_coeff(polyset[(length, d2, typeb)], leftenergy):
                    return (subenergy, d1, d2)
        
    return None

#######################################################################################
#
def findd1d2d3(length, energy, typea, typeb, pairwiseenergy, polyset):
    for d1 in TYPE_ARR:
        polyd1 = polyset[(length, typea, d1)]
        maxenergy = polyd1.degree

        for subenergy in xrange(maxenergy + 1):
            if polyd1.coeff[subenergy] < 1:
                continue
            if energy - subenergy < 0:
                break

            for d2 in TYPE_ARR:
                for d3 in TYPE_ARR:
                    leftenergy = energy - subenergy - pairwiseenergy(d1, d3) - pairwiseenergy(d3, d2)
                    if has_positive_energy_coeff(polyset[(length, d2, typeb)], leftenergy):
                        return (subenergy, d1, d2, d3)

    return None

##################################################################################################
# construct_bounded_energy_dna_list_poly(triplelist, pairwiseenergy):
#     + Constructs and returns the list of DNA strings with the SAME length and bounded energy 
#      satisfied the conditions specified in triplelist. 
#
# Input specification:
#    + triplelist: a list of triples (l, A, B). Each triple (l, A, B) contains information
#        about the DNA we want to construct:
#            - l: the length of the DNA string
#            - A: the lower bound of the free energy of the DNA string
#            - B: the upper bound of the free energy of the DNA string
#    + pairwiseenergy: a function that takes in two inputs X and Y where
#        X, Y in {'A', 'C', 'G', 'T'} and returns the free energy
#        between the ordered pair (X, Y), which is an integer value.
#
# Assumption:
#     + Unlike construct_bounded_energy_dna_list, in this function, we assume all
#        triples must have the same length information.
#
# Detailed Flow:
#    + 
def construct_bounded_energy_dna_list_poly(triplelist, pairwiseenergy):
    numdna = len(triplelist)
    if numdna == 0:
        return []
    if max([triplelist[i][0] for i in xrange(numdna)]) != min([triplelist[i][0] for i in xrange(numdna)]):
        raise RuntimeError("This function only generates bounded DNA strings of EQUAL length")

    length = triplelist[0][0]
    
    polypairlist = build_free_energy_poly(length, pairwiseenergy)

    dnalist = []
    for i in xrange(numdna):
        lowerenergy = triplelist[i][1]
        upperenergy = triplelist[i][2]
        if lowerenergy > upperenergy:
            lowerenergy, upperenergy = upperenergy, lowerenergy

        dnastr = ''
        for energy in xrange(lowerenergy, upperenergy + 1):
            dnastr = extract_free_energy_dna_poly(length, energy, pairwiseenergy, polypairlist)
            if dnastr != '':
                break
        dnalist.append(dnastr)

    return dnalist
