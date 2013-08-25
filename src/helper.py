# Author: Nguyen Truong Duy
# Email: truongduy134@gmail.com

"""
This module constains implementation of various helper functions.
"""

import random

############################################################################
def hamming_dist(str1, str2):
    """
    Compute the hamming distance between str1 and str2

    Inputs:
     + str1, str2: strings
    
    Output:
     + an integer which is the Hamming distance between str1 and str2

    Exception:
     + RuntimeError if str1 and str2 do not have equal length

    Example:
     hamming_dist("ABC", "DEF") returns 3
     hamming_dist("ABC", "ABE") returns 1
     hamming_dist("A", "ABCDEF") raises a RuntimeError
    """
    len1 = len(str1)
    len2 = len(str2)
    
    if len1 != len2:
        raise RuntimeError("Input strings should have equal length")

    numdiff = 0
    for ind in xrange(len1):
        if str1[ind] != str2[ind]:
            numdiff += 1

    return numdiff

######################################################################################
def hamming_dist_substr(str1, start1, end1, str2, start2, end2):
    """
    Compute the hamming distance between substrings str1[start1 ... end1] and str2[start2 ... end2]

    Inputs:
     + str1, str2: Python strings
     + start1, end1: integers indicating the start and end of the substring of str1
     + start2, end2: integers indicating the start and end of the substring of str2

    Output:
     + an integer which is the Hamming distance between str1[start1 ... end1] and str2[start2 ... end2]

    Exception:
     + RuntimeError if the input substrings do have have equal length

    Assumption: 0 <= start1 <= end1 < len(str1) and 0 <= start2 <= end2 < len(str2)
    """
    lensubstr1 = end1 - start1 + 1
    lensubstr2 = end2 - start2 + 1
        
    if lensubstr1 != lensubstr2:
        raise RuntimeError("Input substrings should have equal length")
    
    numdiff = 0
    for offset in xrange(lensubstr1):
        if str1[start1 + offset] != str2[start2 + offset]:
            numdiff += 1

    return numdiff

#########################################################################
# reverse_str(str1):
#
# Time complexity: O(L) where L is the length of the input string
def reverse_str(str1):
    """
    Return a string which is the reverse of the input string

    Input:
     + str1: a Python string
    
    Output:
     + str2: a reverse string of str1 
    """
    return str1[::-1]

############################################################################
# complement_str(str1):
#
# Time complexity: O(L) where L is the length of the input string
def complement_str(str1):
    """
    Find the complement string of the input string

    Input:
     + str1: a Python string

    Output:
     + str2: a complement string of str1 (obtained from str1 by taking the complement of each character in str1)
    """
    chararr = [get_complement_letter(str1[i]) for i in xrange(len(str1))]
    return ''.join(chararr)

############################################################################
def get_complement_letter(mchar):
    """
    Return the complement of the input character

    Input:
     + mchar: input character
    
    Output:
     + the character which is the complement of the input character

    Complementary rules:
     + Complement of '0' is '1' and vice versa
     + Complement of 'A' is 'T' and vice versa
     + Complement of 'G' is 'C' and vice versa
     + For the rest of the characters, their complements are themselves.
    """

    if mchar == '0':
        return '1'
    elif mchar == '1':
        return '0'
    elif mchar == 'A':
        return 'T'
    elif mchar == 'T':
        return 'A'
    elif mchar == 'G':
        return 'C'
    elif mchar == 'C':
        return 'G'
    else:
        return mchar

##############################################################################
# generate_pascal_triangle(l):
#
# Time complexity: O(l^2) (ignoring the complexity of arithmetic)
def generate_pascal_triangle(l):
    """
    Generate a pascal triangle PT with l + 1 rows where PT[r][c] is the value of r C c (r choose c) (0 <= c <= r <= l)

    Input:
     + l: the input which indicates the pascal triangle will have l + 1 rows

    Output
     + PT: the pascal triangle with l + 1 rows. 
    """
    pascal = []
    pascal.append([1])
    prevlen = 1

    for currow in xrange(1, l + 1):
        curlist = [1]
        for ind in xrange(1, prevlen):
            curlist.append(pascal[currow - 1][ind] + pascal[currow - 1][ind - 1])
        curlist.append(1)
        pascal.append(curlist)
        prevlen += 1

    return pascal 

###############################################################################
def get_prefix_sum_pascal_triangle(pascal):
    """
    Modify the input pascal triangle such that pascal[row][col] (new value) = sum of pascal[row][i] (old) for i = 0, ..., col - 1

    Input:
     + pascal: a pascal triangle satisfying the i-th row has exactly (i + 1) column (i >= 0)

    Assumption: the input is correct. That is, pascal is a pascal triangle, or it must at least satisfy: The i-th row has exactly (i + 1) column (i >= 0)
    """

    numrow = len(pascal)

    numcol = 1
    for row in xrange(numrow):
        for col in xrange(1, numcol):
            pascal[row][col] += pascal[row][col - 1]
        numcol += 1

#############################################################################
def generate_prefix_sum_pascal_triangle(l):
    """
    Generate the prefix sum of a pascal triangle

    Input:
     + l: the input which indicates the pascal triangle will have l + 1 rows

    Output
     + PT: a triangle with l + 1 rows such that  PT[r][c] is the sum of r C i (r choose c) for i = 0, ..., c - 1
    """

    pascalprefixsum = generate_pascal_triangle(l)
    get_prefix_sum_pascal_triangle(pascalprefixsum)
    return pascalprefixsum

##############################################################################
# fastpower(a, b):
#
# Assumption: b >= 0 and is an integer
#
# Time complexity: O(log b) (ignoring the time complexity of arithmetics)
def fastpower(a, b):
    """
    Compute and return a**b in logarithmic time

    Inputs:
     + a: an integer representing the base
     + b: a non-negative integer representing the power

    Output:
     + the result of a**b
    """

    if b == 0:
        return 1
    if b == 1:
        return a

    halfpower = fastpower(a, b >> 1)
    result = halfpower * halfpower
    if (b & 1):
        # b is odd
        result *= a
    return result

###############################################################################
def choose_random_pos_list(n, k):
    """
    Choose and return k distinct elements from the set {0, ..., n - 1}

    Inputs:
     + n: an integer indicating the size of the set
     + k: an integer indicating the number of distinct elements to choose

    Output:
     + L: a Python list containing k distinct integers chosen randomly from the set {0, ..., n - 1}
    """

    if n <= 0 or k <= 0 or n < k:
        return []
    if n == k:
        return range(n)

    randomlist = random.sample(range(n), k)
    return randomlist

##############################################################################
# convert_mat_to_strlist(M, n):
#    - Convert a 2D matrix of characters (where each row can be viewed as
#        as a string) into a list of strings
#
# Input description:
#    M: a 2D matrix of characters that has n rows
def convert_mat_to_strlist(M, n):
    """
    Convert a 2D matrix of characters (where each row can be viewed as as a string) into a list of strings

    Inputs:
     + n: an integer indicating the number of rows in the matrix
     + M: a 2D matrix (A Python list of Python list) of characters that has n rows

    Output:
     + strlist: a list of Python strings, each formed by a row of M
    """

    strlist = []
    
    for ind in xrange(n):
        strlist.append(''.join(M[ind]))

    return strlist

##############################################################################
def count_assigned_pos(str1, start1, end1, str2, start2, end2, unkownchar):
    """
    Compute the pair of value (s, t) as decribed below

    Inputs:
     + str1, str2: Python strings
     + start1, end1: integers indicating the start and end of the substring of str1
     + start2, end2: integers indicating the start and end of the substring of str2
     + unknownchar: a character indicating the unknown character.

    Output:
     + a pair of integers (s, t) where 
        - s = the number of positions k in the substrings str1[start1 ... end1] and str2[start2 ... end2] such that str1[k] and str2[k] are not unknown
        - t = the number of positions k in the substrings str1[start1 ... end1] and str2[start2 ... end2] such that str1[k] and str2[k] are not unknown and str1[k] != str2[k]

    Exception:
     + RuntimeError if the input substrings do have have equal length

    Assumption: 0 <= start1 <= end1 < len(str1) and 0 <= start2 <= end2 < len(str2)
    """

    lensub1 = end1 - start1 + 1
    lensub2 = end2 - start2 + 1
    if lensub1 != lensub2:
        raise RuntimeError("Input substrings must have equal length")

    numassigned = 0
    numdiff = 0
    for i in xrange(lensub1):
        if str1[start1 + i] != unkownchar and str2[start2 + i] != unkownchar:
            numassigned += 1
            if str1[start1 + i] != str2[start2 + i]:
                numdiff += 1

    return (numassigned, numdiff)

#######################################################################################
# change_char_in_mat(M, colarr, mapoldtonewchardict):
#
# Time complexity: O(n * m) where 
#    - n is the number of rows in M
#    - m is the size of colarr
def change_char_in_mat(M, colarr, mapoldtonewchardict):
    """
    Replace characters in the columns specified in colarr in the matrix M of characters by new characters. The replacement scheme is specified in the input dictionary mapoldtonewchardict

    Inputs:
     + M: a 2D matrix (A Python list of Python list) of characters
     + colarr: a Python list of integers indicating the column indices in M
     + mapoldtonewchardict: a dictionary that maps old character to new character

    Output:
     + a new matrix in which the characters in the columns specified in colarr in M are replaced by new characters specified in mapoldtonewchardict. If no mapping is available, the characters are unchanged.

    Example:
     To replace all characters '0' by 'A', '1' by 'T' in column 1, 2, 3, 8 of M, call the function:
            change_char_in_mat(M, [1, 2, 3, 8], {'0': 'A', '1': 'T'})
    """

    numrow = len(M)
    if numrow == 0:
        return
    numcol = len(M[0])

    def changechar(mchar):
        if mchar in mapoldtonewchardict:
            return mapoldtonewchardict[mchar]
        return mchar

    for row in xrange(numrow):
        for col in colarr:
            if col >= 0 and col < numcol:
                M[row][col] = changechar(M[row][col])
