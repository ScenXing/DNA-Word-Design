# Author: Nguyen Truong Duy
# Email: truongduy134@gmail.com

import random

############################################################################
# hamming_dist(A, B):
#	- Returns the hamming distance between A and B
#	- Throws an exception if the input strings do have have equal length
def hamming_dist(strA, strB):
	lenA = len(strA)
	lenB = len(strB)
	
	if lenA != lenB:
		raise RuntimeError("Input strings should have equal length")

	numDiff = 0
	for ind in xrange(lenA):
		if strA[ind] != strB[ind]:
			numDiff += 1

	return numDiff

######################################################################################
# hamming_dist_substr(strA, startA, endA, strB, startB, endB):
#	- Returns the hamming distance between the substrings A[startA ... endA] and
#		B[startB ... endB]
#	- Throws an exception if the input substrings do have have equal length
#
# Assumption: [startA ... endA] and [startB ... endB] are valid ranges with
#	- 0 <= startA <= endA < len(strA)
#	- 0 <= startB <= endB < len(strB)
def hamming_dist_substr(strA, startA, endA, strB, startB, endB):
	lenSubstrA = endA - startA + 1
	lenSubstrB = endB - startB + 1
		
	if lenSubstrA != lenSubstrB:
		raise RuntimeError("Input substrings should have equal length")
	
	numDiff = 0
	for offset in xrange(lenSubstrA):
		if strA[startA + offset] != strB[startB + offset]:
			numDiff += 1

	return numDiff

#########################################################################
# reverse_str(strA):
#	- Reverses the input string and returns a new string
#
# Time complexity: O(L) where L is the length of the input string
def reverse_str(strA):
	return strA[::-1]

############################################################################
# complement_str(strA):
#	- Returns the new string that is a complement of the input string
#
# Time complexity: O(L) where L is the length of the input string
def complement_str(strA):
	charArr = [getComplementLetter(strA[i]) for i in xrange(len(strA))]
	return ''.join(charArr)

############################################################################
# getComplementLetter(mChar):
#	- Returns the complement of the input character
#
# Complement of '0' is '1' and vice versa
# Complement of 'A' is 'T' and vice versa
# Complement of 'G' is 'C' and vice versa
# For programming ease, we assume, for the rest of the characters, their
#	complements are themselves.
def getComplementLetter(mChar):
	if mChar == '0':
		return '1'
	elif mChar == '1':
		return '0'
	elif mChar == 'A':
		return 'T'
	elif mChar == 'T':
		return 'A'
	elif mChar == 'G':
		return 'C'
	elif mChar == 'C':
		return 'G'
	else:
		return mChar

##############################################################################
# generatePascalTriangle(l):
#	- Generates a pascal triangle PT with l + 1 rows where
#		PT[r][c] is the value of r C c (r choose c) (0 <= c <= r <= l)
#
# Time complexity: O(l^2) (ignoring the complexity of arithmetic)
def generatePascalTriangle(l):
	pascal = []
	pascal.append([1])
	prevLen = 1

	for curRow in xrange(1, l + 1):
		curList = [1]
		for ind in xrange(1, prevLen):
			curList.append(pascal[curRow - 1][ind] + pascal[curRow - 1][ind - 1])
		curList.append(1)
		pascal.append(curList)
		prevLen += 1

	return pascal 

###############################################################################
# getPrefixSumPascalTriangle(pascal):
#	- Modifies the input pascal triangle such that
#		pascal[row][col] (new value) = sum of pascal[row][i] (old) 
#	  for i = 0, ..., col - 1
#
# Assumption: the input is correct. That is, pascal is a pascal triangle, or it
# 	must at least satisfy: The i-th row has exactly (i + 1) column (i >= 0)
def getPrefixSumPascalTriangle(pascal):
	numRow = len(pascal)

	numCol = 1
	for row in xrange(numRow):
		for col in xrange(1, numCol):
			pascal[row][col] += pascal[row][col - 1]
		numCol += 1

##############################################################################
# fastPower(a, b):
#	- Returns a**b
#
# Assumption: b >= 0 and is an integer
#
# Time complexity: O(log b) (ignoring the time complexity of arithmetics)
def fastPower(a, b):
	if b == 0:
		return 1
	if b == 1:
		return a

	halfPower = fastPower(a, b >> 1)
	result = halfPower * halfPower
	if (b & 1):
		# b is odd
		result *= a
	return result

###############################################################################
# chooseRandomPosList(n, k):
#	- Chooses k distinct elements from the set {0, ..., n - 1}
#	- Returns the chosen list
def chooseRandomPosList(n, k):
	if n <= 0 or k <= 0 or n < k:
		return []
	if n == k:
		return range(n)

	randomList = random.sample(range(n), k)
	return randomList

##############################################################################
# convertMatToStrSet(M, n):
#	- Converts a 2D matrix of characters (where each row can be viewed as
#		as a string) into a set of strings
#
# Input description:
#	M: a 2D matrix of characters that has n rows
def convertMatToStrSet(M, n):
	strSet = []
	
	for ind in xrange(n):
		strSet.append(''.join(M[ind]))

	return strSet

##############################################################################
# countAssignedPos(strA, startA, endA, strB, startB, endB, unknownChar):
#	+ Return the pair of value (s, t) where
#		s = the number of positions k in the substrings strA[startA ... endA] and
#	            strB[startB ... endB] such that strA[k] and strB[k] are not unknown
#		t = the number of positions k in the substrings strA[startA ... endA] and
#	            strB[startB ... endB] such that strA[k] and strB[k] are not unknown and
#		    strA[k] != strB[k]
#	+ Throws an exception if the substrings strA[startA ... endA] and strB[startB ... endB]
#	  do not have equal length
#
def countAssignedPos(strA, startA, endA, strB, startB, endB, unknownChar):
	lenSubA = endA - startA + 1
	lenSubB = endB - startB + 1
	if lenSubA != lenSubB:
		raise RuntimeError("Input substrings must have equal length")

	numAssigned = 0
	numDiff = 0
	for i in xrange(lenSubA):
		if strA[startA + i] != unknownChar and strB[startB + i] != unknownChar:
			numAssigned += 1
			if strA[startA + i] != strB[startB + i]:
				numDiff += 1

	return (numAssigned, numDiff)

#######################################################################################
# changeCharInMat(M, colArr, mapOldToNewCharDict):
#	+ Replaces characters in those columns in colArr in the matrix of characters
#	  by new characters. The replacement scheme is specified in the input dictionary
#	  mapOldToNewCharDict
#
# Input description:
#	+ M: a 2D matrix of characters
#	+ mapOldToNewCharDict: a dictionary that maps old character to new character
#
# Example: To replace all characters '0' by 'A', '1' by 'T' in column 1, 2, 3, 8 of M,
#	call the function:
#		changeCharInMat(M, [1, 2, 3, 8], {'0': 'A', '1': 'T'})  
def changeCharInMat(M, colArr, mapOldToNewCharDict):
	numRow = len(M)
	if numRow == 0:
		return
	numCol = len(M[0])

	def changeChar(mChar):
		if mChar in mapOldToNewCharDict:
			return mapOldToNewCharDict[mChar]
		return mChar

	for row in xrange(numRow):
		for col in colArr:
			if col >= 0 and col < numCol:
				M[row][col] = changeChar(M[row][col])
