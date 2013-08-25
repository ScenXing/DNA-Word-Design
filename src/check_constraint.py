# Author: Nguyen Truong Duy
# Email: truongduy134@gmail.com

"""
This module is for testing purpose. For example, it contains various function to check if a given list of DNA words satisfies a subset of constraints.
"""

# Self-written python modules
import helper
import algo_subroutine
import free_energy_routine

# Standard modules
import math

####################################################################################
# checkc1(strlist, k1):
#
# Time complexity: O(n^2 * l) where 
#   + n is the size of the string list, and
#   + l is the length of each string in the list
def checkc1(strlist, k1):
    """
    Check if the input string list satisfies C1 constraint (Basic Hamming Constraint)

    Inputs:
    + strlist: a list of Python strings
    + k1: parameter of C1 constraint

    Output:
    + True if the input string list satisfies C1(k1)
    + False otherwise

    Exception:
    + a RuntimeError will be raised if the input strings do not have equal length
    """

    numstr = len(strlist)
        
    for indone in xrange(numstr):
        for indtwo in xrange(indone + 1, numstr):
            if not checkc1each(strlist[indone], strlist[indtwo], k1):
                return False
    return True

# checkc1each(str1, str2, k1):
#
# Time complexity: O(l) where l is the length of str1, and str2
def checkc1each(str1, str2, k1):
    """
    Check if the input string pair satisfies C1 constraint (Basic Hamming Constraint)

    Inputs:
    + str1, str2: Python strings
    + k1: parameter of C1 constraint

    Output:
    + True if the input string pair satisfies C1(k1), i.e. their Hamming distance is at least k1
    + False otherwise

    Exception:
    + a RuntimeError will be raised if str1 and str2 do not have equal length
    """

    return helper.hamming_dist(str1, str2) >= k1

##########################################################################################
# checkc2(strlist, k2):
#
# Time complexity: O(n^2 * l) where 
#   + n is the size of the string list, and
#   + l is the length of each string in the list
def checkc2(strlist, k2):
    """
    Check if the input string list satisfies C2 constraint (Reverse Complementary Constraint)

    Inputs:
    + strlist: a list of Python strings
    + k2: parameter of C2 constraint

    Output:
    + True if the input string list satisfies C2(k2)
    + False otherwise

    Exception:
    + a RuntimeError will be raised if the input strings do not have equal length
    """

    numstr = len(strlist)

    for indone in xrange(numstr):
        for indtwo in xrange(numstr):
            if indone != indtwo:
                if not checkc2each(strlist[indone], strlist[indtwo], k2):
                    return False
    return True

# checkc2each(str1, str2, k2):
#
# Time complexity: O(l) where l is the length of str1, and str2 
def checkc2each(str1, str2, k2):
    """
    Check if the input string pair satisfies C2 constraint (Reverse Complementary Constraint)

    Inputs:
    + str1, str2: Python strings
    + k2: parameter of C2 constraint

    Output:
    + True if the input string pair satisfies C2(k2), i.e. their Hamming distance between str1 and the reverse of str2 is at least k2
    + False otherwise

    Exception:
    + a RuntimeError will be raised if str1 and str2 do not have equal length
    """

    revcompstr2 = helper.reverse_str(helper.complement_str(str2))
    return helper.hamming_dist(str1, revcompstr2) >= k2

###############################################################################################
# checkc3(strlist, k3):
#
# Time complexity: O(n * l) where 
#       + n is the size of the string list, and
#       + l is the length of each string in the list
def checkc3(strlist, k3):
    """
    Check if the input string list satisfies C3 constraint (Self Reverse Complementary Constraint)

    Inputs:
    + strlist: a list of Python strings
    + k3: parameter of C3 constraint

    Output:
    + True if the input string list satisfies C3(k3)
    + False otherwise
    """

    numstr = len(strlist)

    for ind in xrange(numstr):
        if not checkc3each(strlist[ind], k3):
            return False
    return True

# checkc3each(str1, k3):
#
# Time complexity: O(l) where l is the length of str1
def checkc3each(str1, k3):
    """
    Check if the input string satisfies C3 constraint (Self Reverse Complementary Constraint )

    Inputs:
    + str1: a Python string
    + k3: parameter of C3 constraint

    Output:
    + True if the input string satisfies C3(k3), i.e. the Hamming distance between str1 and is reverse is at least k3
    + False otherwise
    """

    if k3 <= 0:
        return True
    return checkc2each(str1, str1, k3)

################################################################################################
# checkc4(strlist, k4):
#
# Time complexity: O(n^2 * l * k4) where 
#   + n is the size of the string list, and
#   + l is the length of each string in the list
def checkc4(strlist, k4):
    """
    Check if the input string list satisfies C4 constraint (Shifting Hamming Constraint)

    Inputs:
    + strlist: a list of Python strings
    + k4: parameter of C4 constraint

    Output:
    + True if the input string list satisfies C4(k4)
    + False otherwise

    Exception:
    + a RuntimeError will be raised if the input strings do not have equal length
    """

    numstr = len(strlist)

    for indone in xrange(numstr):
        for indtwo in xrange(numstr):
            if indone != indtwo:
                if not checkc4each(strlist[indone], strlist[indtwo], k4):
                    return False
    return True

# checkc4each(str1, str2, k4):
#
# Time complexity: O(l * k4) where l is the length of str1, and of str2
def checkc4each(str1, str2, k4):
    """
    Check if the input string pair satisfies C4 constraint (Shifting Hamming Constraint )

    Inputs:
    + str1, str2: Python strings
    + k4: parameter of C4 constraint

    Output:
    + True if the input string pair satisfies C4(k4)
    + False otherwise

    Exception:
    + a RuntimeError will be raised if str1 and str2 do not have equal length
    """

    len1 = len(str1)
    len2 = len(str2)
    if len1 != len2:
        raise RuntimeError("Input strings should have equal length")
    
    mlen = len1
    # Check corner cases
    if k4 <= 0:
        return True
    if k4 > mlen:
        return False

    for sublen in xrange(mlen - k4, mlen + 1):
        if sublen == 0:
            continue
        hammingdist = helper.hamming_dist_substr(str1, 0, sublen - 1, str2, mlen - sublen, mlen - 1)
        if hammingdist < (k4 - (mlen - sublen)):
            return False
    
    return True

############################################################################################################
# checkc5(strlist, k5):
#
# Time complexity: O(n^2 * l * k5) where 
#    + n is the size of the string list, and
#    + l is the length of each string in the list    
def checkc5(strlist, k5):
    """
    Check if the input string list satisfies C5 constraint (Shifting Reverse Complementary Constraint)

    Inputs:
    + strlist: a list of Python strings
    + k5: parameter of C5 constraint

    Output:
    + True if the input string list satisfies C5(k5)
    + False otherwise

    Exception:
     + a RuntimeError will be raised if the input strings do not have equal length
    """

    numstr = len(strlist)
    
    for indone in xrange(numstr):
        for indtwo in xrange(numstr):
            if indone != indtwo:
                if not checkc5each(strlist[indone], strlist[indtwo], k5):
                    return False
    return True

# checkc5each(str1, str2, k5):
#
# Time complexity: O(l * k5) where l is the length of str1, or str2
def checkc5each(str1, str2, k5):
    """
    Check if the input string pair satisfies C5 constraint (Shifting Reverse Complementary Constraint)

    Inputs:
    + str1, str2: Python strings
    + k5: parameter of C5 constraint

    Output:
    + True if the input string pair satisfies C5(k5)
    + False otherwise

    Exception:
     + a RuntimeError will be raised if str1 and str2 do not have equal length
    """

    len1 = len(str1)
    len2 = len(str2)
    if len1 != len2:
        raise RuntimeError("Input strings should have equal length")
    
    mlen = len1
    # Check corner cases
    if k5 <= 0:
        return True
    if k5 > mlen:
        return False

    revcompstr2 = helper.reverse_str(helper.complement_str(str2))
    for sublen in xrange(mlen - k5, mlen + 1):
        if sublen == 0:
            continue

        hammingdist = helper.hamming_dist_substr(str1, 0, sublen - 1, revcompstr2, 0, sublen - 1)
        if hammingdist < (k5 - (mlen - sublen)):
            return False

        hammingdist = helper.hamming_dist_substr(str1, mlen - sublen, mlen - 1, revcompstr2, mlen - sublen, mlen - 1)
        if hammingdist < (k5 - (mlen - sublen)):
            return False
    
    return True

#########################################################################################################
# checkc6(strlist, k6):
#
# Time complexity: O(n * l * k6) where 
#    + n is the size of the string list, and
#    + l is the length of each string in the list
def checkc6(strlist, k6):
    """
    Check if the input string list satisfies C6 constraint (Shifting Self Reverse Complementary Constraint)

    Inputs:
    + strlist: a list of Python strings
    + k6: parameter of C6 constraint

    Output:
    + True if the input string list satisfies C6(k6)
    + False otherwise
    """

    for mstr in strlist:
        if not checkc6each(mstr, k6):
            return False
    return True

# checkc6each(str1, k6):
#
# Time complexity: O(l * k6) where l is the length of str1
def checkc6each(str1, k6):
    """
    Check if the input string satisfies C6 constraint (Shifting Self Reverse Complementary Constraint)

    Inputs:
    + str1: a Python string
    + k6: parameter of C6 constraint

    Output:
    + True if the input string satisfies C6(k6)
    + False otherwise
    """
    return checkc5each(str1, str1, k6)

##################################################################################################
# checkc7(strlist, ratio):
#
# Time complexity: O(n * l) where 
#    + n is the size of the string list, and
#    + l is the length of each string in the list    
def checkc7(strlist, ratio):
    """
    Check if the input string list satisfies C7 constraint (GC Content Constraint)

    Inputs:
    + strlist: a list of Python strings
    + ratio: parameter of C7 constraint

    Output:
    + True if the input string list satisfies C7(ratio)
    + False otherwise
    """

    for mstr in strlist:
        if not checkc7each(mstr, ratio):
            return False
    return True

# checkc7each(str1, ratio):
#    - Returns true if the input string satisfies GC Content Constraint
#    - Returns false otherwise
def checkc7each(str1, ratio):
    """
    Check if the input string satisfies C7 constraint (GC Content Constraint)

    Inputs:
    + str1: a Python string
    + ratio: parameter of C7 constraint

    Output:
    + True if the input string satisfies C7(ratio)
    + False otherwise
    """
    if ratio < 0 or ratio > 1:
        return False

    len1 = len(str1)
    nummustpresentceil = int(math.ceil(ratio * len1))
    nummustpresentfloor = int(ratio * len1)

    numpresent = 0
    for mchar in str1:
        if mchar == 'G' or mchar == 'C':
            numpresent += 1

    return numpresent == nummustpresentfloor or numpresent == nummustpresentceil

##################################################################################################
# checkc8(strlist, maxlenrun):
#
# Time complexity: O(n * l) where 
#    + n is the size of the string list, and
#    + l is the length of each string in the list
def checkc8(strlist, maxlenrun):
    """
    Check if the input string list satisfies C8 constraint (Consecutive Base Constraint Constraint)

    Inputs:
    + strlist: a list of Python strings
    + maxlenrun: parameter of C8 constraint

    Output:
    + True if the input string list satisfies C8(maxlenrun)
    + False otherwise
    """

    if maxlenrun <= 0:
        return False
    for mstr in strlist:
        if not checkc8each(mstr, maxlenrun):
            return False
    return True

# checkc8each(str1, maxlenrun):
#
# Time complexity: O(l) where l is the length of str1
def checkc8each(str1, maxlenrun):
    """
    Check if the input string satisfies C8 constraint (Consecutive Base Constraint Constraint)

    Inputs:
    + str1: a Python string
    + maxlenrun: parameter of C8 constraint

    Output:
    + True if the input string satisfies C8(maxlenrun)
    + False otherwise
    """
    if maxlenrun <= 0:
        return False

    startind = 0
    len1 = len(str1)
    while startind < len1:
        lenrun = 1
        endind = startind + 1
        while endind < len1 and str1[startind] == str1[endind]:
            lenrun += 1
            endind += 1
        
        if lenrun > maxlenrun:
            return False
        startind = endind        

    return True

########################################################################################
# checkc9(strlist, sigma, pairwise_energy):
#
# Time complexity: O(n * l + n^2) where 
#    + n is the size of the string list, and
#    + l is the length of each string in the list    
def checkc9(strlist, sigma, pairwise_energy):
    """
    Check if the input string list satisfies C9 constraint (Free Energy Constraint)

    Inputs:
    + strlist: a list of Python strings
    + sigma: parameter of C9 constraint
    + pairwise_energy: a function that takes in two inputs X and Y where X, Y in {'A', 'C', 'G', 'T'} and returns a non-negative integer indicating the free energy between the ordered pair (X, Y).

    Output:
    + True if the input string list satisfies C9(sigma) with respect to the pairwise energy function
    + False otherwise
    """

    numstr = len(strlist)

    energylist = [free_energy_routine.compute_free_energy(strlist[ind], pairwise_energy) for ind in xrange(numstr)]

    for indone in xrange(numstr):
        for indtwo in xrange(indone + 1, numstr):
            if math.fabs(energylist[indone] - energylist[indtwo]) > sigma:
                return False
    return True

###########################################################################################
def checkconstraintset(strlist, maptypetoparam):
    """
    Check if the input string list satisfies a set of constraints.

    Inputs:
    + strlist: a list of Python strings
    + maptypetoparam: a dictionary that maps an integer (representing the constraint type) to the parameter corresponding to that constraint type. Note that if maptypetoparam contains 9 as a key, then maptypetoparam[9] should be a pair whose first entry is a parameter for C9 constraint, and second entry is a pairwise free energy function.

    Output:
    + True if the input string list satisfies the set of constraints specified in maptypetoparam
    + False otherwise
    """

    checkfunclist = [checkc1, checkc2, checkc3, checkc4, checkc5, checkc6, checkc7, checkc8, checkc9]
    
    result = True
    for typekey in maptypetoparam:
        if typekey >= 1 and typekey <= 9:
            if typekey != 9:
                result = checkfunclist[typekey - 1](strlist, maptypetoparam[typekey])
            else:
                result = checkfunclist[typekey - 1](strlist, maptypetoparam[typekey][0], maptypetoparam[typekey][1])

            if not result:
                print typekey
                return False

    return True

###########################################################################################
def checkgendnaword1to6and9algo(strlist, maptypetoparam):
    """
    Check if the input string list satisfies C1 through C6 constraint, and C9 constraint

    Inputs:
    + strlist: a list of Python strings
    + maptypetoparam: a dictionary that maps an integer (representing the constraint type) to the parameter corresponding to that constraint type. It must have 1, 2, ..., 6 and 9 as keys. Note that maptypetoparam[9] is a 4 x 4 matrix of non-negative integers where rows and columns are indexed by typeArr = ['A', 'C', 'G', 'T'] in that order, presenting the pairwise free energy values.

    Output:
    + True if the input string list satisfies C1 through C6 constraints and C9(sigma) with respect to the pairwise energy function where sigma = 4 * (maxpairwise - minpairwise) + maxpairwise, and maxpairwise and minpairwise are the largest and smallest entries in maptypetoparam[9] respectively 
    + False otherwise
    """

    pairwise_energy = free_energy_routine.create_pairwise_energy_func(maptypetoparam[9])
    maxpairwise = free_energy_routine.find_extreme_pairwise_energy(pairwise_energy, max)
    minpairwise = free_energy_routine.find_extreme_pairwise_energy(pairwise_energy, min)
    D = maxpairwise - minpairwise
    checkmaptypetoparam = maptypetoparam.copy()
    checkmaptypetoparam[9] = (4 * D + maxpairwise, pairwise_energy)
    return checkconstraintset(strlist, checkmaptypetoparam)

###########################################################################################
def checklength14(l, n, k1, k4):
    """
    Check if the input values l, n, k1, k4 satisfy the following two conditions:
       1) l >= 2 * k where k = max(k1, k4)
       2) l - k * log(e) - k * log(l / k) - 2 * log(n) - 2 * log(k) > 0

    Inputs:
    + l, n, k1, k4: integers

    Output:
    + True if l, n, k1, k1 satisfy the above two conditions
    + False otherwise
    """

    k = max(k1, k4)
    if l < 2 * k:
        return False

    expr = l - k * math.log(math.e, 2) - k * math.log(l * 1.0 / k, 2) - 2 * math.log(n, 2) - 2 * math.log(k, 2)
    if expr > 0:
        return True
    return False

#####################################################################################
def checkexpcountempty(n, l, k1, k4):
    """
    Check if the input values l, n, k1, k4 satisfy that
        ExpCount(M, k1, k4) > n * (n - 1) / 2 * (1 + 2(k4 - 1)) - 1
    where M is a partially assigned matrix of size n x l where every entry is an unkown

    Inputs:
    + l, n, k1, k4: integers

    Output:
    + True if l, n, k1, k1 satisfy the above condition
    + False otherwise
    """
    pascalprefixsum = helper.generate_prefix_sum_pascal_triangle(l)

    expcount = algo_subroutine.compute_expcount_empty(n, l, k1, k4, pascalprefixsum)
    rhs = (n * (n - 1) / 2) * (1 + 2 * (k4 - 1)) - 1
    return expcount > rhs
