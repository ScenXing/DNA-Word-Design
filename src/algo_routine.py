# Author: Nguyen Truong Duy
# Email: truongduy134@gmail.com

#
# Note: All the functions in this module are related to the following papers
#
#   1) "Deterministic Polynomial-Time Algorithms for Designing
#       Short DNA Words" by Kao et al.
#       
#       You can retrieve a copy of this paper at:
#           http://arxiv.org/pdf/1201.6358.pdf 
#
#   2) (For free-energy constraint C9) "Randomized Fast Design of Short DNA
#      Words" by Kao et al.
#
#      You can retreive a copy of this paper at:
#           http://dl.acm.org/citation.cfm?id=1597047
#

"""
This module contains implementation of various algorithms to generate a set of short DNA words satisfying the following constraints:
1) C1 and C4
2) C1 to C6
3) C1 to C7
4) C1, C2, C3, C7 and C8
5) C1 to C8
6) C1 to C6, and C9
"""

# Self-written modules
import helper
import algo_subroutine
import free_energy_routine
import m_fraction

# Builtin modules
import math

##########################################################################################
# genbinword14sub(n, l, k1, k4, pascalprefixsum):
#
# Assumption: (n, l) must satisfy the condition of Lemma 3 in paper (1), i.e.
#            ExpCount(M, k1, k4) > nC2 * (1 + 2 * (k4 - 1)) - 1
#    so that such a matrix M can exist. Otherwise, the function will
#    throw a Runtime Error
#
# The implementation strictly follows Algorithm 1 presented in paper (1)
#
# Detail:
def genbinword14sub(n, l, k1, k4, pascalprefixsum):
    """
    Compute and return a BINARY (k1, k4)-distance matrix M (consisting of characters '0' and '1') of dimension n x l

    Inputs:
     + n: the number of DNA strings to generate
     + l: the length of each DNA string to generate
     + k1: parameter of C1 constraint
     + k4: parameter of C4 constraint
     + pascalprefixsum: a 2D array with at least l + 1 rows. Row i has (i + 1) entries. pascalprefixsum[i][j] = sum ( i Choose h ) for h = 0, ..., j.

    Output:
     + M: a BINARY (k1, k4)-distance matrix M of dimension n x l  
    """

    # 'X' stands for unknown
    M = [['X'] * l for i in xrange(n)]


    numdiff1 = [[0] * (i + 1) for i in xrange(n)]
    numdiff4 = [[[0] * n for i in xrange(n)] for k in xrange(k4 - 1)]

    prevexpcount = algo_subroutine.compute_expcount_empty(n, l, k1, k4, pascalprefixsum)

    guardval = m_fraction.MyFraction(((n * (n - 1)) / 2) * (1 + 2 * (k4 - 1)) - 1)
    for strpos in xrange(l):
        for strid in xrange(n):
            if prevexpcount <= guardval:
                raise RuntimeError("Your Algo or Input (n, l) is wrong")

            expcount0 = prevexpcount + algo_subroutine.compute_change_in_expcount(M, strpos, strid, '0', 'X', numdiff1, numdiff4, n, l, k1, k4, pascalprefixsum)
            expcount1 = m_fraction.MyFraction(2) * prevexpcount - expcount0
            
            if expcount0 >= expcount1:
                M[strid][strpos] = '0'
                prevexpcount = expcount0
            else:
                M[strid][strpos] = '1'
                prevexpcount = expcount1
            algo_subroutine.update_numdiff(M, strpos, strid, 'X', numdiff1, numdiff4, n, l, k4)

    return M

##################################################################################
# genbinword14(n, k1, k4):
#
# Algorithm flow:
#    + Precompute pascalprefixsum which is
#        - a 2D array with at least l + 1 rows. Row i has (i + 1) entries
#        - pascalprefixsum[i][j] = sum ( i Choose h ) for h = 0, ..., j
#      This is to save computation time in other subroutines.
#    + Find the minimum length l (via binary search) that satisfies Lemma 3
#      of paper (1)
#    + Call the routine genbinword14sub(n, l, k1, k4, pascalprefixsum)
#      (implemeted above)
#
# Time complexity: O(n^2 * (max(k1, k4) + log n)^2)
def genbinword14(n, k1, k4):
    """
    Compute and return the (k1, k4)-distance matrix M (consisting of characters '0' and '1'). Each row of M can be viewed as a string, and M can be viewed as a set of binary strings satisfying C1 and C4 constraints.

    Inputs:
     + n: the number of DNA strings to generate
     + k1: parameter of C1 constraint
     + k4: parameter of C4 constraint 

    Output:
     + M: a BINARY (k1, k4)-distance matrix M with n rows 
    """

    initlen = algo_subroutine.initlength14(n, k1, k4)
    # Precompute the prefix sum of computations, which will be used a lot later
    pascalprefixsum = helper.generate_prefix_sum_pascal_triangle(initlen)

    # Find the minimum length that satisfies Lemma 3 of paper (1) by Binary Search
    minL = algo_subroutine.binarysearchlength14(n, k1, k4, initlen, pascalprefixsum)

    M = genbinword14sub(n, minL, k1, k4, pascalprefixsum)

    return M

######################################################################################
# gendnaword14(n, maptypetoparam):
#
# Detail:
#    1) Generate a BINARY (k1, k4) distance matrix M (consisting of characters 
#       '0' and '1'). Each row in M can be viewed as a string that satisfies 
#       C1(k1) and C4(k4)
#    2) For each entry in M, change '0' to 'C' and '1' to 'G'.
#    3) Then the list {W(0), ..., W(n - 1)} (W(i) is a string formed by the i-th row of M)
#       is a list of DNA words satifying C1 to C4
def gendnaword14(n, maptypetoparam):
    """
    Generate and return a set of DNA words satisfying C1 and C4 constraints.

    Inputs:
     + n: the number of strings to generate
     + maptypetoparam: a dictionary that maps an integer (representing the constraint type) to the parameter corresponding to that constraint type. It must have 1 and 4 as keys

    Output:
     + strlist: a Python list of DNA words satisfying C1 and C4 constraints

    Example: to generate a set of 15 DNA words satisfying C1(8) and C4(9), call the function
            gendnaword14(15, {1 : 8, 4 : 9 })
    """

    if n <= 0:
        return []

    M = genbinword14(n, maptypetoparam[1], maptypetoparam[4])
    
    helper.change_char_in_mat(M, range(len(M[0])), {'0': 'C', '1': 'G'})
    return helper.convert_mat_to_strlist(M, n)

######################################################################################
# gendnaword1to6(n, maptypetoparam):
#
# Implementation is based on Lemma 12 in paper (1)
#
# Detail: (Note that k(i) = maptypetoparam[i])
#    1) Generate a BINARY (k1, k4) distance matrix M (consisting of characters 
#       '0' and '1'). Each row in M can be viewed as a string that satisfies 
#       C1(k1) and C4(k4)
#    2) For each entry in M, change '0' to 'A' and '1' to 'T'.
#    3) Let k = max{k2, k3, k5, k6}. Add k copies of 'C' at the beginning of each 
#       word formed by each row i in M, called such word W(i).
#       Then the list {W(0), ..., W(n - 1)} is a list of DNA words satifying C1 to C6 
def gendnaword1to6(n, maptypetoparam):
    """
    Generate and return a set of DNA words satisfying C1 and C6 constraints.

    Inputs:
     + n: the number of strings to generate
     + maptypetoparam: a dictionary that maps an integer (representing the constraint type) to the parameter corresponding to that constraint type. It must have 1, 2, ..., 6 as keys

    Output:
     + strlist: a Python list of DNA words satisfying C1 and C6 constraints

    Example: to generate a set of 25 DNA words satisfying C1(8), C2(4), C3(5), C4(7), C5(8), C6(10), call the function
        gendnaword1to6(25, {1 : 8, 2 : 4, 3 : 5, 4 : 7, 5 : 8, 6 : 10})
    """

    if n <= 0:
        return []

    M = genbinword14(n, maptypetoparam[1], maptypetoparam[4])

    # In M, change '0' to 'A' and '1' to 'T'
    helper.change_char_in_mat(M, range(len(M[0])), {'0': 'A', '1' : 'T'})
    
    k = max(maptypetoparam[2], maptypetoparam[3], maptypetoparam[5], maptypetoparam[6])
    leadingstr = 'C' * k

    # Append k = max(k2, k3, k5, k6) to the start of each row in M
    strlist = []
    for row in xrange(n):
        strlist.append(leadingstr + ''.join(M[row]))
    
    return strlist

######################################################################################
# gendnaword1to7(n, maptypetoparam):
#
# Implementation is based on Lemma 14 in paper (1)
#
# Detail:
#    1) Generate a BINARY (k1, k4) distance matrix M (consisting of characters 
#       '0' and '1'). Each row in M can be viewed as a string that satisfies 
#       C1(k1) and C4(k4)
#    2) Let k = max{k2, k3, k5, k6} (k(i) = maptypetoparam[i]). 
#       For each row of M, add k copies of '1' at the beginning and k copies of '1'
#       at the end.
#    3) Let l be the new length of each row in M now. Let gamma = maptypetoparam[7]. 
#       Assume 0 <= gamma <= 1. Choose randomly a subset of size ceil(gamma * l) of
#       the set {0, ..., l - 1}, let it be P
#    4) For each row of M, for positions in P, change '0' to 'C' and '1' to 'G'.
#       For the remaining positions NOT in P, change '0' to 'A' and '1' to 'T'
#    5) The set of strings formed by rows in the new M satifies C1 to C7.
def gendnaword1to7(n, maptypetoparam):
    """
    Generate and return a set of DNA words satisfying C1 and C7 constraints.

    Inputs:
     + n: the number of strings to generate
     + maptypetoparam: a dictionary that maps an integer (representing the constraint type) to the parameter corresponding to that constraint type. It must have 1, 2, ..., 7 as keys.

    Outputs:
     + strlist: a Python list of DNA words satisfying C1 and C7 constraints

    Example: to generate a set of 25 DNA words satisfying C1(8), C2(4), C3(5), C4(7), C5(8), C6(10), C7(0.7), call the function
        gendnaword1to7(25, {1 : 8, 2 : 4, 3 : 5, 4 : 7, 5 : 8, 6 : 10, 7 : 0.7})
    """

    if n <= 0:
        return []
    if maptypetoparam[7] > 1:
        return []

    M = genbinword14(n, maptypetoparam[1], maptypetoparam[4])
    
    k = max(maptypetoparam[2], maptypetoparam[3], maptypetoparam[5], maptypetoparam[6])
    l = len(M[0]) + k + k

    chosencolumn = []
    if int(maptypetoparam[7]) == 1:
        chosencolumn = range(l)
    else:
        chosencolumn = helper.choose_random_pos_list(l, int(math.ceil(maptypetoparam[7] * l)))
    allcolumn = range(l)

    strlist = []
    for row in xrange(n):
        # Append k instances of '1' at the beginning and the end of each row in M
        newlist = ['1'] * k
        newlist.extend(M[row])
        newlist.extend(['1'] * k)

        helper.change_char_in_mat([newlist], chosencolumn, {'0': 'C', '1': 'G'})
        helper.change_char_in_mat([newlist], allcolumn, {'0': 'A', '1': 'T'})
        strlist.append(''.join(newlist))

    return strlist

######################################################################################
# gendnaword12378(n, maptypetoparam):
#
# Implementation is based on Lemma 16 in paper (1)
#
# Detail: (Note that k(i) = maptypetoparam[i])
#    1) Generate a BINARY (k1, 1) distance matrix M (consisting of characters 
#       '0' and '1'). Each row in M can be viewed as a string that satisfies 
#       C1(k1) and C4(1). This is essentially equivalent to the fact that each
#       row of M satisfy C1(k1) only.
#    2) Let l0 be the number of columns in M. If l0 is odd, we append '0' at the
#       end of each row in M.
#    3) Let k = max(k2, k3). 
#       For each row of M, add k copies of '1' at the beginning and k copies of '1'
#       at the end.
#    4) Let S be a list of strings formed by rows in M. Apply breakrun function to 
#       each string in S so that all strings in S do not have runs of the same
#       characters longer than maptypetoparam[8].
#    5) Let l be the new length of each string in S after Step 4. Let gamma = maptypetoparam[7]. 
#       Assume 0 <= gamma <= 1. Choose randomly a subset of size ceil(gamma * l) of
#       the set {0, ..., l - 1}, let it be P
#    6) For each string in S, for positions in P, change '0' to 'C' and '1' to 'G'.
#       For the remaining positions NOT in P, change '0' to 'A' and '1' to 'T'
#    7) The new list S contains DNA words satisfying C1, C2, C3, C7 and C8 constraints. 
def gendnaword12378(n, maptypetoparam):
    """
    Generate and return a set of DNA words satisfying C1, C2, C3, C7 and C8 constraints.

    Inputs:
     + n: the number of strings to generate
     + maptypetoparam: a dictionary that maps an integer (representing the constraint type) to the parameter corresponding to that constraint type. It must have 1, 2, 3, 7 and 8 as keys.

    Output:
     + strlist: a Python list of DNA words satisfying C1, C2, C3, C7 and C8 constraints

    Example: to generate a set of 50 DNA words satisfying C1(8), C2(4), C3(5), C7(0.7), C8(3), call the function
        gendnaword12378(50, {1 : 8, 2 : 4, 3 : 5, 7 : 0.7, 8 : 3})
    """

    if n <= 0:
        return []
    if maptypetoparam[7] > 1:
        return []

    M = genbinword14(n, maptypetoparam[1], 1)
    l0 = len(M[0])
    if l0 & 1:
        # If l0 is odd, append '0' at the end of every word in M so that
        # the new length is even
        for row in xrange(n):
            M[row].append('0')

    strlist = []
    k = max(maptypetoparam[2], maptypetoparam[3])
    for row in xrange(n):
        newlist = ['1'] * k
        newlist.extend(M[row])
        newlist.extend(['1'] * k)
        strlist.append(''.join(newlist))

    # Break run
    for strid in xrange(n):
        strlist[strid] = algo_subroutine.breakrun(strlist[strid], maptypetoparam[8])

    newlen = len(strlist[0])
    chosencolumn = []
    allcolumn = range(newlen)
    if int(maptypetoparam[7]) == 1:
        chosencolumn = range(newlen)
    else:
        chosencolumn = helper.choose_random_pos_list(newlen, int(math.ceil(maptypetoparam[7] * newlen)))

    for strid in xrange(n):
        curlist = list(strlist[strid])
        helper.change_char_in_mat([curlist], chosencolumn, {'0': 'C', '1': 'G'})
        helper.change_char_in_mat([curlist], allcolumn, {'0': 'A', '1':'T'})
        
        strlist[strid] = ''.join(curlist)

    return strlist

#######################################################################################
# gendnaword1to8(n, maptypetoparam):
#
# Implementation is a combination of 
#    + gendnaword1to8algo2 function (implemented based on Lemma 20 and Theorem 21 in paper (1))
#    + gendnaword1to8algo1 function (implemented based on Lemma 18 and Theorem 19 in paper (1))
#
# Description:
#    + If 1 / (d + 1) <= gamma <= d / (d + 1), we will use gendnaword1to8algo1 to generate
#      the set of DNA words. It is because generally gendnaword1to8algo1 produces shorter words
#    + Otherwise, if d >= 2, we will use gendnaword1to8algo2
#    + Otherwise, a RuntimeError will be thrown
def gendnaword1to8(n, maptypetoparam):
    """
    Generate and return a set of DNA words satisfying C1 through C8 constraints.

    Inputs:
     + n: the number of strings to generate
     + maptypetoparam: a dictionary that maps an integer (representing the constraint type) to the parameter corresponding to that constraint type. It must have 1, 2, ..., 8 as keys.

    Output:
     + strlist: a Python list of DNA words satisfying C1 through C8 constraints

    Example: to generate a set of 50 DNA words satisfying C1(8), C2(4), C3(5), C4(7), C5(8), C6(10), C7(0.7), C8(3), call the function
            gendnaword1to8(50, {1 : 8, 2 : 4, 3 : 5, 4 : 7, 5 : 8, 6 : 10, 7 : 0.7, 8 : 3})
    """

    if n <= 0:
        return []

    gamma = maptypetoparam[7]
    if gamma > 1:
        return []
    
    d = maptypetoparam[8]
    if d < 2:
        raise RuntimeError("gendnaword1to8algo2 only works with maxlenrun >= 2")

    if (1.0 / (d + 1) > gamma) or (gamma > d * 1.0 / (d + 1)):
        return gendnaword1to8algo2(n, maptypetoparam)
    return gendnaword1to8algo1(n, maptypetoparam)

########################################################################################
# gendnaword1to8algo2(n, maptypetoparam):
#
# The implementation is based on Lemma 20 and Theorem 21 in paper (1)
#
# Detail: (Note that k(i) = maptypetoparam[i] and d = maptypetoparam[8])
#    1) Generate a BINARY (max(k1, k4), 1) distance matrix M (consisting of characters 
#       '0' and '1'). Each row in M can be viewed as a string that satisfies 
#       C1(max(k1, k4)) and C4(1). This is essentially equivalent to the fact that each
#       row of M satisfy C1(max(k1, k4)) only.
#    2) For each row in M, after every d - 1 characters, or when we reach the last character
#       of a row, insert a new character which is the complementary ('0' is complementary to '1', 
#       and vice versa) to the immediate previous character.
#    3) Add 1 copy of '1' at the beginning and 1 copy of '0' at the end of each row in M.
#    4) Let k = max{k2, k3, k4, k5, k6}. For each row of M, add
#       ceil(k / d) copies of the length-d string 11...10 at the beginning and the end of each
#       row in M.
#    5) Let l be the new length of each row in M now (or the number of columns). 
#       Let gamma = maptypetoparam[7]. Assume 0 <= gamma <= 1. 
#       Choose randomly a subset of size ceil(gamma * l) of the set {0, ..., l - 1}, let it be P
#    6) For each row of M, for positions in P, change '0' to 'C' and '1' to 'G'.
#       For the remaining positions NOT in P, change '0' to 'A' and '1' to 'T'
#    7) The list of strings formed by rows in the new M satifies C1 through C8.
def gendnaword1to8algo2(n, maptypetoparam):
    """
    Generate and return a set of DNA words satisfying C1 through C8 constraints.

    Inputs:
     + n: the number of strings to generate
     + maptypetoparam: a dictionary that maps an integer (representing the constraint type) to the parameter corresponding to that constraint type. It must have 1, 2, ..., 8 as keys.

    Output:
     + strlist: a Python list of DNA words satisfying C1 through C8 constraints

    Exception:
     + A RuntimeError will be raised if maptypetoparam[8] < 2
    """

    if n <= 0:
        return []
    if maptypetoparam[7] > 1:
        return []

    if maptypetoparam[8] < 2:
        raise RuntimeError("gendnaword1to8algo2 only works with maxlenrun >= 2")

    # Generate the set of strings satisfies C1 constraint
    M = genbinword14(n, max(maptypetoparam[1], maptypetoparam[4]), 1)

    k = max([maptypetoparam[i] for i in xrange(2, 6 + 1)])

    newM = []
    l0 = len(M[0])

    # Prepare the 'string' (list of characters) used later
    baselist = ['1'] * (maptypetoparam[8] - 1)
    baselist.append('0')
    numtime = int(math.ceil(k * 1.0 / maptypetoparam[8]))
    supplist = baselist * numtime

    for row in xrange(n):
        newrow = []
        newrow.extend(supplist)
        newrow.append('1')
        sublen = 0
        for ind in xrange(l0):
            newrow.append(M[row][ind])
            sublen += 1
            if (sublen == maptypetoparam[8] - 1) or (ind == l0 - 1):
                newrow.append(helper.get_complement_letter(M[row][ind]))
                sublen = 0
        newrow.append('0')
        newrow.extend(supplist)
        newM.append(newrow)


    newlen = len(newM[0])
    allcolumn = range(newlen)
    if maptypetoparam[7] == 1:
        chosencolumn = range(newlen)
    else:
        chosencolumn = helper.choose_random_pos_list(newlen, int(math.ceil(maptypetoparam[7] * newlen)))

    helper.change_char_in_mat(newM, chosencolumn, {'0': 'C', '1': 'G'})
    helper.change_char_in_mat(newM, allcolumn, {'0': 'A', '1': 'T'})

    return helper.convert_mat_to_strlist(newM, n)

########################################################################################
# gendnaword1to8algo1(n, maptypetoparam):
#
# Assumption: Let gamma = maptypetoparam[7] and d = maptypetoparam[8]. We must have 
#                1 / (d + 1) <= gamma <= d / (d + 1)
#             Otherwise, an RunTimeError will be raised
#
# The implementation is based on Lemma 18 and Theorem 19 in paper (1)
#
# Detail: (Note that d = maptypetoparam[8])
#    1) Generate a BINARY (max(k1, k4), 1) distance matrix M (consisting of characters 
#       '0' and '1'). Each row in M can be viewed as a string that satisfies 
#       C1(max(k1, k4)) and C4(1). This is essentially equivalent to the fact that each
#       row of M satisfy C1(max(k1, k4)) only.
#    2) Let k = max{k2, k3, k4, k5, k6}. 
#       For each row of M, add k copies of '1' at the beginning and k copies of '1'
#       at the end.
#    3) Let l be the new length of rows in M (or the number of columns in M). 
#       Partition the integer interval [0, l - 1] into subintervals Z(1), Z(2), ..., Z(s)
#       for some s such that
#        (1) Each subinterval consists of at most d integers and at least one integer
#        (2) The total number of integers in the odd-indexed subintervals is 
#            floor(newgamma * l) where newgamma = max(gamma, 1 - gamma)
#
#       The step is not straightforward to implement. So we describe how it is implemented in
#       function
#        + Let newgamma = max(gamma, 1 - gamma)
#        + Since 1 / (d + 1) <= gamma <= d / (d + 1), we have
#                1/2 <= newgamma <= d / (d + 1)
#        + Let numchoose = floor(newgamma * l) and numnotchoose = l - numchoose
#         + Let q = floor(numchoose / numnotchoose) and r = numchoose % numnotchoose
#        + It can be shown that 1 <= q <= d. And when q = d, r = 0
#        + Instead of finding Z(1), ..., Z(s) explicitly, we just need to maintain 2 list
#            - oddlist is the union odd-indexed subintervals Z(1), Z(3), ...
#            - evenlist is the union of even-indexed subintervals Z(2), Z(4), ...
#        + Starting with i = 0, we repeated the following procedure:
#            - Add the next q integers including i (i.e. i, i + 1, ..., i + q - 1) 
#              or add until we reach l - 1 into oddlist
#            - If we reach (add) l - 1, stop the procedure.
#            - Otherwise, if r is positive, add the next integer, i.e. i + q into oddlist,
#              and decrement r
#            - Add the next integer (if it is less than l) into evenlist
#            - Update i for the next iteration
#
#    4) For each row in M:
#           (1) If gamma < 0.5: 
#        + For entries whose position is in an odd-indexed subinterval, change '0' to 'A', 
#          '1' to 'T'.
#        + For entries whose position is in an even-indexed subinterval, change '0' to 'C',
#          '1' to 'G'.
#        (2) If gamma >= 0.5: 
#        + For entries whose position is in an odd-indexed subinterval, change '0' to 'C', 
#          '1' to 'G'.
#        + For entries whose position is in an even-indexed subinterval, change '0' to 'A',
#          '1' to 'T'.
#    5) The list of strings formed by rows in the new M satisfies C1 through C8.
def gendnaword1to8algo1(n, maptypetoparam):
    """
    Generate and return a set of DNA words satisfying C1 through C8 constraints.

    Inputs:
     + n: the number of strings to generate
     + maptypetoparam: a dictionary that maps an integer (representing the constraint type) to the parameter corresponding to that constraint type. It must have 1, 2, ..., 8 as keys.

    Output:
     + strlist: a Python list of DNA words satisfying C1 through C8 constraints

    Exception:
     + A RuntimeError will be raised if the following condition is NOT satified: 
                1 / (d + 1) <= gamma <= d / (d + 1)
     where gamma = maptypetoparam[7] and d = maptypetoparam[8]
    """

    gamma = maptypetoparam[7]
    if n <= 0 or gamma > 1:
        return []

    d = maptypetoparam[8]
    if (1.0 / (d + 1) > gamma) or (gamma > d * 1.0 / (d + 1)):
        raise RuntimeError("gendnaword1to8algo1 works only if 1 / (d + 1) <= gamma <= d / (d + 1)")

        
    # Generate the set of strings satisfies C1 constraint
    M = genbinword14(n, max(maptypetoparam[1], maptypetoparam[4]), 1)

    k = max([maptypetoparam[i] for i in xrange(2, 6 + 1)])

    newM = []
    for row in xrange(n):
        newrow = []
        newrow.extend(['1'] * k)
        newrow.extend(M[row])
        newrow.extend(['1'] * k)
        newM.append(newrow)

    # Find oddlist and evenlist as mentioned in Step 3 (see comments above)
    newlen = len(newM[0])
    newgamma = gamma
    if newgamma < 0.5:
        newgamma = 1 - newgamma
    numchoose = int(newgamma * newlen)
    numnotchoose = newlen - numchoose
    minoddsize = int(numchoose * 1.0 / numnotchoose)
    numleft = numchoose % numnotchoose
    oddlist = []
    evenlist = []
    ind = 0
    oddsize = 0
    while ind < newlen:
        oddlist.append(ind)
        ind += 1
        oddsize += 1
        if ind < newlen and oddsize == minoddsize:
            if numleft != 0:
                oddlist.append(ind)
                numleft -= 1
                ind += 1
            if ind < newlen:
                evenlist.append(ind)
                ind += 1
                oddsize = 0

    # Convert binary words into DNA words
    if gamma < 0.5:
        helper.change_char_in_mat(newM, oddlist, {'0': 'A', '1': 'T'})
        helper.change_char_in_mat(newM, evenlist, {'0': 'C', '1': 'G'})
    else:
        helper.change_char_in_mat(newM, oddlist, {'0': 'C', '1': 'G'})
        helper.change_char_in_mat(newM, evenlist, {'0': 'A', '1': 'T'})

    return helper.convert_mat_to_strlist(newM, n)

##############################################################################################
# gendnaword1to6And9(n, maptypetoparam):
#
#    Refer to gendnaword1to6and9generic comments
#
# The function we use to construct bounded energy DNA strings here is
#    free_energy_routine.construct_bounded_energy_dna_list, which employs dynamic programming
#    approach
#
def gendnaword1to6and9(n, maptypetoparam):
    """
    Generate and return a set of DNA words satisfying C1 to C6 with specified parameters, and C9(4D + Gamma_max) constraints where
      - Gamma_max and Gamma_min are the largest and smallest among 16 entries of maptypetoparam[9].
      - D = Gamma_max - Gamma_min
    In this implementation, we use dynamic programming to solve the bounded energy strand generation problem.

    See gendnaword1to6and9generic for more information
    """
    return gendnaword1to6and9generic(n, maptypetoparam, free_energy_routine.construct_bounded_energy_dna_list)

##############################################################################################
# gendnaword1to6and9poly(n, maptypetoparam)
#
#    Refer to gendnaword1to6and9generic comments
#
# The function we use to construct bounded energy DNA strings here is
#    construct_bounded_energy_dna_list_poly, which employs polynomial multiplication and
#    generating function approach
#
def gendnaword1to6and9poly(n, maptypetoparam):
    """
    Generate and return a set of DNA words satisfying C1 to C6 with specified parameters, and C9(4D + Gamma_max) constraints where
      - Gamma_max and Gamma_min are the largest and smallest among 16 entries of maptypetoparam[9].
      - D = Gamma_max - Gamma_min
    In this implementation, we use an approach of polynomial multiplication and generating function to solve the bounded energy strand generation problem. 

    See gendnaword1to6and9generic for more information
    """
    return gendnaword1to6and9generic(n, maptypetoparam, free_energy_routine.construct_bounded_energy_dna_list_poly)

##############################################################################################
# gendnaword1to6and9generic(n, maptypetoparam, construct_bounded_energy_dna):
#   + Generate and return a set of DNA words satisfying C1 to C6, and C9 constraints.
#
# Input specfication:
#   + n: the number of strings to generate
#   + maptypetoparam: a dictionary that maps an integer (representing
#       the constraint type) to the parameter corresponding to that
#       constraint type.
#     In this function, maptypetoparam must contain
#       maptypetoparam[i] = k(i) (an integer for C(i) constraint)
#     for i = 1, ..., 6
#
#       maptypetoparam[9]: is a 4 x 4 matrix where rows and columns are indexed
#     by the list ['A', 'C', 'G', 'T'] in that order, presenting the pairwise free energy values.
#     For example, if
#       maptypetoparam[9] = M = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 16, 16]]
#     Then the free energy between 'T' and 'T' is M[3][3] = 16, and
#          the free energy between 'A' and 'G' is M[0][2] = 3, and
#          the free energy between 'G' and 'A' is M[2][0] = 9, etc.
#     All 16 entries of maptypetoparam[9] should be non-negative integers!
#
#   + construct_bounded_energy_dna: a function to generate a list of string with specified length
#       and that has free energy bounded between the specified interval.
#
# Implementation is based on Lemma 5, 6, 7, 8 of paper (2)
#
# Detail:
#   1) Generate a list S of DNA words satisfying C1 through C6 constraints by calling gendnaword1to6
#   2) Let femax and femin be the largest and smallest free energy of DNA words in S.
#   3) If femax - femin <= 4 * D + Gamma_max, return S.
#   4) Otherwise, let l the length of strings in S. 
#      Find feminstr which is the minimum free energy among all DNA words of length 2l
#   5) Let alpha = femax + feminstr, and beta = alpha + 2 * D. For each string S(i) in S, construct
#      a DNA word T(i) with length 2l and has free energy between alpha - FE(S(i)) and beta - FE(S(i))
#      (FE(W) denotes the free energy of a word W). This can be done by calling 
#      construct_bounded_energy_dna
#   6) Let U(i) be the string constructed by concatenating T(i)[0 ... l - 1] (the first half of T(i))
#      and S(i) and T(i)[l ... 2l - 1] (the second half of T(i))
#   7) The list of strings U(0), ..., U(n - 1) satisfying C1 through C6 constraints, and
#      C9(4 * D + Gamma_max)
def gendnaword1to6and9generic(n, maptypetoparam, construct_bounded_energy_dna):
    """
    Generate and return a set of DNA words satisfying C1 to C6 with specified parameters, and C9(4D + Gamma_max) constraints where
      - Gamma_max and Gamma_min are the largest and smallest among 16 entries of maptypetoparam[9].
      - D = Gamma_max - Gamma_min

    Inputs:
     + n: the number of strings to generate
     + maptypetoparam: a dictionary that maps an integer (representing the constraint type) to the parameter corresponding to that constraint type. It must have 1, 2, ..., 6 and 9 as keys. Note that maptypetoparam[9] is a 4 x 4 matrix of non-negative integers where rows and columns are indexed by typeArr = ['A', 'C', 'G', 'T'] in that order, presenting the pairwise free energy values.
     + construct_bounded_energy_dna: a function to solve the bounded energy strand generation problem. It generates a list of string with specified length and that has free energy bounded between the specified interval.

    Output:
     + strlist: a Python list of DNA words satisfying C1 through C6 constraints, and C9(4D + Gamma_max)

    Example:
    + The input pairwise free energy matrix M is given by 
            M = [[5, 4, 6, 1], [2, 10, 3, 4], [6, 11, 5, 8], [1, 3, 4, 8]]. 
    + To generate a list of 95 DNA words satisfying C1(3), C2(6), C3(2), C4(7), C5(8), C6(9) and C9 constraint, call the function:
            strlist = gen_dna_words(95, {1: 3, 2: 6, 3: 2, 4: 7, 5: 8, 6: 9, 9: M})
    """

    # Generate a set of words satisfying C1 to C6 constraint
    strlist16 = gendnaword1to6(n, maptypetoparam)

    # Create the pairwise energy function from the 4 x 4 matrices
    pairwiseenergy = free_energy_routine.create_pairwise_energy_func(maptypetoparam[9])
    maxpairwise = free_energy_routine.find_extreme_pairwise_energy(pairwiseenergy, max)
    minpairwise = free_energy_routine.find_extreme_pairwise_energy(pairwiseenergy, min)
    D = maxpairwise - minpairwise

    # Compute the free energy of each string
    freeenergylist = [free_energy_routine.compute_free_energy(strlist16[i], pairwiseenergy) for i in xrange(n)]
    femax = max(freeenergylist)
    femin = min(freeenergylist)

    if femax - femin <= 4 * D + maxpairwise:
        return strlist16

    newlen = len(strlist16[0])
    newlen += newlen

    alpha = femax + free_energy_routine.compute_min_free_energy(newlen, pairwiseenergy)
    beta = alpha + D + D
    triplelist =  [(newlen, alpha - freeenergylist[i], beta - freeenergylist[i]) for i in xrange(n)]

    # Find the list of bounded energy DNA strings
    boundedestrlist = construct_bounded_energy_dna(triplelist, pairwiseenergy)

    strlist = []
    halflen = newlen / 2
    for ind in xrange(n):
        newstr = boundedestrlist[ind][0 : halflen] + strlist16[ind] + boundedestrlist[ind][halflen : newlen]
        strlist.append(newstr)

    return strlist 
