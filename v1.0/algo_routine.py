# Author: Nguyen Truong Duy
# Email: truongduy134@gmail.com

#
# Note: All the functions in this module are related to the following papers
#
#	1) "Deterministic Polynomial-Time Algorithms for Designing
#	    Short DNA Words" by Kao et al.
#	    
#	    You can retrieve a copy of this paper at:
#			http://arxiv.org/pdf/1201.6358.pdf 
#
#	2) (For free-energy constraint C9) "Randomized Fast Design of Short DNA
#	   Words" by Kao et al.
#
#	   You can retreive a copy of this paper at:
#			http://dl.acm.org/citation.cfm?id=1597047
#

# Self-written modules
import helper
import algo_subroutine
import free_energy_routine
import m_fraction

# Builtin modules
import math

############################################################################
# genBinWord14Sub(n, l, k1, k4, pascalPrefixSum):
#	+ Computes and returns the BINARY (k1, k4)-distance matrix M
#	  of dimension n x l
#
# Input description:
#	+ pascalPrefixSum: is a 2D array with at least l + 1 rows. Row i has (i + 1) entries
#	  pascalPrefixSum[i][j] = sum ( iCh ) for h = 0, ..., j
# Note: iCh denotes i Choose h
#
# Assumption: (n, l) must satisfy the condition of Lemma 3 in paper (1), i.e.
#			ExpCount(M, k1, k4) > nC2 * (1 + 2 * (k4 - 1)) - 1
#	so that such a matrix M can exist. Otherwise, the function will
#	throw a Runtime Error
#
# The implementation strictly follows Algorithm 1 presented in paper (1)
def genBinWord14Sub(n, l, k1, k4, pascalPrefixSum):
	# 'X' stands for unknown
	M = [['X'] * l for i in xrange(n)]


	numDiff1 = [[0] * (i + 1) for i in xrange(n)]
	numDiff4 = [[[0] * n for i in xrange(n)] for k in xrange(k4 - 1)]

	prevExpCount = algo_subroutine.computeExpCountEmpty(n, l, k1, k4, pascalPrefixSum)

	guardVal = m_fraction.MyFraction(((n * (n - 1)) / 2) * (1 + 2 * (k4 - 1)) - 1)
	for strPos in xrange(l):
		for strId in xrange(n):
			if prevExpCount <= guardVal:
				raise RuntimeError("Your Algo or Input (n, l) is wrong")

			expCount0 = prevExpCount + algo_subroutine.computeChangeInExpCount(M, strPos, strId, '0', 'X', numDiff1, numDiff4, n, l, k1, k4, pascalPrefixSum)
			expCount1 = m_fraction.MyFraction(2) * prevExpCount - expCount0
			
			if expCount0 >= expCount1:
				M[strId][strPos] = '0'
				prevExpCount = expCount0
			else:
				M[strId][strPos] = '1'
				prevExpCount = expCount1
			algo_subroutine.updateNumDiff(M, strPos, strId, 'X', numDiff1, numDiff4, n, l, k4)

	return M

##################################################################################
# genBinWord14(n, k1, k4):
#	+ Computes and returns the BINARY (k1, k4)-distance matrix M.
#	  Each row of M can be viewed as a string, and M can be viewed as a set
#	  of binary strings satisfying C1 and C4 constraints.
#
# Algorithm flow:
#	+ Precompute pascalPrefixSum which is
#		- a 2D array with at least l + 1 rows. Row i has (i + 1) entries
#		- pascalPrefixSum[i][j] = sum ( iCh ) for h = 0, ..., j
#	  Note: iCh denotes i Choose h
#	  This is to save computation time in other subroutines.
#	+ Find the minimum length l (via binary search) that satisfies Lemma 3
#	  of paper (1)
#	+ Call the routine genBinWord14Sub(n, l, k1, k4, pascalPrefixSum)
#	  (implemeted above)
#
# Time complexity: O(n^2 * (max(k1, k4) + log n)^2)
def genBinWord14(n, k1, k4):
	initL = algo_subroutine.initLength14(n, k1, k4)
	# Precompute the prefix sum of computations, which will be used a lot later
	pascalPrefixSum = helper.generatePascalTriangle(initL)
	helper.getPrefixSumPascalTriangle(pascalPrefixSum)

	# Find the minimum length that satisfies Lemma 3 of paper (1) by Binary Search
	minL = algo_subroutine.binarySearchLength14(n, k1, k4, initL, pascalPrefixSum)

	M = genBinWord14Sub(n, minL, k1, k4, pascalPrefixSum)

	return M

######################################################################################
# genDnaWord14(n, mapTypeToParam):
#	+ Generates and returns a set of DNA words satisfying C1 and C4 constraints.
#
# Input specfication:
#	+ n: the number of strings in the set
#	+ mapTypeToParam: a dictionary that maps an integer (representing
#		the constraint type) to the parameter corresponding to that
#		constraint type.
#	  In this function, mapTypeToParam must contain
#		mapTypeToParam[1] = k1 (an integer for C1 constraint)
#		mapTypeToParam[4] = k4 (an integer for C4 constraint)
#
# For example, to generate a set of 15 DNA words satisfying C1(8) and C4(9), we
#	call the function
#		genDnaWord14(15, {1 : 8, 4 : 9 })
def genDnaWord14(n, mapTypeToParam):
	if n <= 0:
		return []

	M = genBinWord14(n, mapTypeToParam[1], mapTypeToParam[4])
	
	helper.changeCharInMat(M, range(len(M[0])), {'0': 'C', '1': 'G'})
	return helper.convertMatToStrSet(M, n)

######################################################################################
# genDnaWord1To6(n, mapTypeToParam):
#	+ Generates and returns a set of DNA words satisfying C1 to C6 constraints.
#
# Input specfication:
#	+ n: the number of strings in the set
#	+ mapTypeToParam: a dictionary that maps an integer (representing
#		the constraint type) to the parameter corresponding to that
#		constraint type.
#	  In this function, mapTypeToParam must contain
#		mapTypeToParam[i] = k(i) (an integer for C(i) constraint)
#	  for i = 1, ..., 6
#
# For example, to generate a set of 25 DNA words satisfying C1(8), C2(4), C3(5),
#	C4(7), C5(8), C6(10), call the function
#		genDnaWord1To6(25, {1 : 8, 2 : 4, 3 : 5, 4 : 7, 5 : 8, 6 : 10})
#
# Implementation is based on Lemma 12 in paper (1)
#
# Detail:
#	+ 
def genDnaWord1To6(n, mapTypeToParam):
	if n <= 0:
		return []

	M = genBinWord14(n, mapTypeToParam[1], mapTypeToParam[4])

	# In M, change '0' to 'A' and '1' to 'T'
	helper.changeCharInMat(M, range(len(M[0])), {'0': 'A', '1' : 'T'})
	
	k = max(mapTypeToParam[2], mapTypeToParam[3], mapTypeToParam[5], mapTypeToParam[6])
	leadingStr = 'C' * k

	# Append k = max(k2, k3, k5, k6) to the start of each row in M
	strSet = []
	for row in xrange(n):
		strSet.append(leadingStr + ''.join(M[row]))
	
	return strSet

######################################################################################
# genDnaWord1To7(n, mapTypeToParam):
#	+ Generates and returns a set of DNA words satisfying C1 to C7 constraints.
#
# Input specfication:
#	+ n: the number of strings in the set
#	+ mapTypeToParam: a dictionary that maps an integer (representing
#		the constraint type) to the parameter corresponding to that
#		constraint type.
#	  In this function, mapTypeToParam must contain
#		mapTypeToParam[i] = k(i) (an integer for C(i) constraint)
#	  for i = 1, ..., 6
#		mapTypeToParam[7] = gamma (a real number in [0, 1] that indicating
#	  the amount of GC content, i.e. exactly ceiling(gamma * l) characters in
#	  the generated string are 'G' or 'C') (C7 constraint)
#
# For example, to generate a set of 25 DNA words satisfying C1(8), C2(4), C3(5),
#	C4(7), C5(8), C6(10), C7(0.7) call the function
#		genDnaWord1To7(25, {1 : 8, 2 : 4, 3 : 5, 4 : 7, 5 : 8, 6 : 10, 7 : 0.7})
#
# Implementation is based on Lemma 14 in paper (1)
#
# Detail:
#	+
def genDnaWord1To7(n, mapTypeToParam):
	if n <= 0:
		return []
	if mapTypeToParam[7] > 1:
		return []

	M = genBinWord14(n, mapTypeToParam[1], mapTypeToParam[4])
	
	k = max(mapTypeToParam[2], mapTypeToParam[3], mapTypeToParam[5], mapTypeToParam[6])
	l = len(M[0]) + k + k

	chosenColumn = []
	if int(mapTypeToParam[7]) == 1:
		chosenColumn = range(l)
	else:
		chosenColumn = helper.chooseRandomPosList(l, int(math.ceil(mapTypeToParam[7] * l)))
 	allColumn = range(l)

	strSet = []
	for row in xrange(n):
		# Append k instances of '1' at the beginning and the end of each row in M
		newList = ['1'] * k
		newList.extend(M[row])
		newList.extend(['1'] * k)

		helper.changeCharInMat([newList], chosenColumn, {'0': 'C', '1': 'G'})
		helper.changeCharInMat([newList], allColumn, {'0': 'A', '1': 'T'})
		strSet.append(''.join(newList))

	return strSet

######################################################################################
# genDnaWord12378(n, mapTypeToParam):
#	+ Generates and returns a set of DNA words satisfying C1, C2, C3
#		C7 and C8 constraints.
#
# Input specfication:
#	+ n: the number of strings in the set
#	+ mapTypeToParam: a dictionary that maps an integer (representing
#		the constraint type) to the parameter corresponding to that
#		constraint type.
#	  In this function, mapTypeToParam must contain
#		mapTypeToParam[i] = k(i) (an integer for C(i) constraint)
#	  for i = 1, 2, 3
#
#		mapTypeToParam[7] = gamma (a real number in [0, 1] that indicating
#	  the amount of GC content, i.e. exactly ceiling(gamma * l) characters in
#	  the generated string are 'G' or 'C') (C7 constraint)
#
#		mapTypeToParam[8] = d (an integer for C8 constraint)
#
# For example, to generate a set of 50 DNA words satisfying C1(8), C2(4), C3(5),
#	C7(0.7), C8(3) call the function
#		genDnaWord12378(25, {1 : 8, 2 : 4, 3 : 5, 7 : 0.7, 8 : 3})
#
# Implementation is based on Lemma 16 in paper (1)
#
# Detail:
#	+
def genDnaWord12378(n, mapTypeToParam):
	if n <= 0:
		return []
	if mapTypeToParam[7] > 1:
		return []

	M = genBinWord14(n, mapTypeToParam[1], 1)
	l0 = len(M[0])
	if l0 & 1:
		# If l0 is odd, append '0' at the end of every word in M so that
		# the new length is even
		for row in xrange(n):
			M[row].append('0')

	strSet = []
	k = max(mapTypeToParam[2], mapTypeToParam[3])
	for row in xrange(n):
		newList = ['1'] * k
		newList.extend(M[row])
		newList.extend(['1'] * k)
		strSet.append(''.join(newList))

	# Break run
	for strId in xrange(n):
		strSet[strId] = algo_subroutine.breakRun(strSet[strId], mapTypeToParam[8])

	newLen = len(strSet[0])
	chosenColumn = []
	allColumn = range(newLen)
	if int(mapTypeToParam[7]) == 1:
		chosenColumn = range(newLen)
	else:
		chosenColumn = helper.chooseRandomPosList(newLen, int(math.ceil(mapTypeToParam[7] * newLen)))

	for strId in xrange(n):
		curList = list(strSet[strId])
		helper.changeCharInMat([curList], chosenColumn, {'0': 'C', '1': 'G'})
		helper.changeCharInMat([curList], allColumn, {'0': 'A', '1':'T'})
		
		strSet[strId] = ''.join(curList)

	return strSet

#######################################################################################
# genDnaWord1To8(n, mapTypeToParam):
#	+ Generates and returns a set of DNA words satisfying C1 to C8 constraints.
#
# Input specfication:
#	+ n: the number of strings in the set
#	+ mapTypeToParam: a dictionary that maps an integer (representing
#		the constraint type) to the parameter corresponding to that
#		constraint type.
#	  In this function, mapTypeToParam must contain
#		mapTypeToParam[i] = k(i) (an integer for C(i) constraint)
#	  for i = 1, ..., 6
#
#		mapTypeToParam[7] = gamma (a real number in [0, 1] that indicating
#	  the amount of GC content, i.e. exactly ceiling(gamma * l) characters in
#	  the generated string are 'G' or 'C') (C7 constraint)
#
#		mapTypeToParam[8] = d (an integer for C8 constraint)
#
# For example, to generate a set of 50 DNA words satisfying C1(8), C2(4), C3(5),
#	C4(7), C5(8), C6(10), C7(0.7), C8(3) call the function
#		genDnaWord12378(25, {1 : 8, 2 : 4, 3 : 5, 4 : 7, 5 : 8, 6 : 10, 7 : 0.7, 8 : 3})
#
# Implementation is a combination of 
#	+ genDnaWord1To8Algo2 function (implemented based on Lemma 20 and Theorem 21 in paper (1))
#	+ genDnaWord1To8Algo1 function (implemented based on Lemma 18 and Theorem 19 in paper (1))
#
# Description:
#	+ If 1 / (d + 1) <= gamma <= d / (d + 1), we will use genDnaWord1To8Algo1 to generate
#	  the set of DNA words. It is because generally genDnaWord1To8Algo1 produces shorter words
#	+ Otherwise, if d >= 3, we will use genDnaWord1To8Algo2
#	+ Otherwise, a RuntimeError will be thrown
#
def genDnaWord1To8(n, mapTypeToParam):
	if n <= 0:
		return []

	gamma = mapTypeToParam[7]
	if gamma > 1:
		return []
	
	d = mapTypeToParam[8]
	if (1.0 / (d + 1) > gamma) or (gamma > d * 1.0 / (d + 1)):
		return genDnaWord1To8Algo2(n, mapTypeToParam)
	return genDnaWord1To8Algo1(n, mapTypeToParam)

########################################################################################
# genDnaWord1To8Algo2(n, mapTypeToParam):
#	+ Generates and returns a set of DNA words satisfying C1 to C8 constraints.
#
# The implementation is based on Lemma 20 and Theorem 21 in paper (1)
#
# Detail:
def genDnaWord1To8Algo2(n, mapTypeToParam):
	if n <= 0:
		return []
	if mapTypeToParam[7] > 1:
		return []

	if mapTypeToParam[8] < 3:
		raise RuntimeError("genDnaWord1To8Algo2 only works with maxLenRun >= 3")

	# Generate the set of strings satisfies C1 constraint
	M = genBinWord14(n, max(mapTypeToParam[1], mapTypeToParam[4]), 1)

	k = max([mapTypeToParam[i] for i in xrange(2, 6 + 1)])

	newM = []
	l0 = len(M[0])

	# Prepare the 'string' (list of characters) used later
	baseList = ['1'] * (mapTypeToParam[8] - 1)
	baseList.append('0')
	numTime = int(math.ceil(k * 1.0 / (mapTypeToParam[8] - 2)))
	suppList = baseList * numTime

	for row in xrange(n):
		newRow = []
		newRow.extend(suppList)
		newRow.append('1')
		subLen = 0
		for ind in xrange(l0):
			newRow.append(M[row][ind])
			subLen += 1
			if (subLen == mapTypeToParam[8] - 1) or (ind == l0 - 1):
				newRow.append(helper.getComplementLetter(M[row][ind]))
				subLen = 0
		newRow.append('0')
		newRow.extend(suppList)
		newM.append(newRow)


	newLen = len(newM[0])
	allColumn = range(newLen)
	if mapTypeToParam[7] == 1:
		chosenColumn = range(newLen)
	else:
		chosenColumn = helper.chooseRandomPosList(newLen, int(math.ceil(mapTypeToParam[7] * newLen)))

	helper.changeCharInMat(newM, chosenColumn, {'0': 'C', '1': 'G'})
	helper.changeCharInMat(newM, allColumn, {'0': 'A', '1': 'T'})

	return helper.convertMatToStrSet(newM, n)

########################################################################################
# genDnaWord1To8Algo1(n, mapTypeToParam):
#	+ Generates and returns a set of DNA words satisfying C1 to C8 constraints.
#
# The implementation is based on Lemma 18 and Theorem 19 in paper (1)
#
# Detail:
def genDnaWord1To8Algo1(n, mapTypeToParam):
	gamma = mapTypeToParam[7]
	if n <= 0 or gamma > 1:
		return []

	d = mapTypeToParam[8]
	if (1.0 / (d + 1) > gamma) or (gamma > d * 1.0 / (d + 1)):
		raise RuntimeError("genDnaWord1To8Algo1 works only if 1 / (d + 1) <= gamma <= d / (d + 1)")

		
	# Generate the set of strings satisfies C1 constraint
	M = genBinWord14(n, max(mapTypeToParam[1], mapTypeToParam[4]), 1)

	k = max([mapTypeToParam[i] for i in xrange(2, 6 + 1)])

	newM = []
	for row in xrange(n):
		newRow = []
		newRow.extend(['1'] * k)
		newRow.extend(M[row])
		newRow.extend(['1'] * k)
		newM.append(newRow)

	newLen = len(newM[0])
	newGamma = gamma
	if newGamma < 0.5:
		newGamma = 1 - newGamma
	numChoose = int(newGamma * newLen)
	numNotChoose = newLen - numChoose
	minOddSize = int(numChoose * 1.0 / numNotChoose)
	numLeft = numChoose % numNotChoose
	oddList = []
	evenList = []
	ind = 0
	oddSize = 0
	while ind < newLen:
		oddList.append(ind)
		ind += 1
		oddSize += 1
		if oddSize == minOddSize:
			if numLeft != 0:
				oddList.append(ind)
				numLeft -= 1
				ind += 1
			evenList.append(ind)
			ind += 1
			oddSize = 0

	if gamma < 0.5:
		helper.changeCharInMat(newM, oddList, {'0': 'A', '1': 'T'})
		helper.changeCharInMat(newM, evenList, {'0': 'C', '1': 'G'})
	else:
		helper.changeCharInMat(newM, oddList, {'0': 'C', '1': 'G'})
		helper.changeCharInMat(newM, evenList, {'0': 'A', '1': 'T'})

	return helper.convertMatToStrSet(newM, n)

##############################################################################################
# genDnaWord1To6And9(n, mapTypeToParam):
#
#	Refer to genDnaWord1To6And9Generic comments
#
# The function we use to construct bounded energy DNA strings here is
#	free_energy_routine.constructBoundedEnergyDnaList, which employs dynamic programming
#	approach
#
def genDnaWord1To6And9(n, mapTypeToParam):
	return genDnaWord1To6And9Generic(n, mapTypeToParam, free_energy_routine.constructBoundedEnergyDnaList)

##############################################################################################
# genDnaWord1To6And9Poly(n, mapTypeToParam)
#
#	Refer to genDnaWord1To6And9Generic comments
#
# The function we use to construct bounded energy DNA strings here is
#	constructBoundedEnergyDnaListPoly, which employs polynomial multiplication and
#	generating function approach
#
def genDnaWord1To6And9Poly(n, mapTypeToParam):
	return genDnaWord1To6And9Generic(n, mapTypeToParam, free_energy_routine.constructBoundedEnergyDnaListPoly)

##############################################################################################
# genDnaWord1To6And9Generic(n, mapTypeToParam, constructBoundedEDnaFunc):
#	+ Generates and returns a set of DNA words satisfying C1 to C6, and C9 constraints.
#
# Input specfication:
#	+ n: the number of strings in the set
#	+ mapTypeToParam: a dictionary that maps an integer (representing
#		the constraint type) to the parameter corresponding to that
#		constraint type.
#	  In this function, mapTypeToParam must contain
#		mapTypeToParam[i] = k(i) (an integer for C(i) constraint)
#	  for i = 1, ..., 6
#
#		mapTypeToParam[9]: is a 4 x 4 matrix where rows and columns are indexed
#	  by typeArr = ['A', 'C', 'G', 'T'] in that order, presenting the pairwise free energy values.
#	  For example, if
#		mapTypeToParam[9] = M = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 16, 16]]
#	  Then the free energy between 'T' and 'T' is M[3][3] = 16, and
#	       the free energy between 'A' and 'G' is M[0][2] = 3, and
#	       the free energy between 'G' and 'A' is M[2][0] = 9, etc.
#	  All 16 entries of mapTypeToParam[9] should be non-negative integers!
#
#	+ constructBoundedEDnaFunc: the function to generate 
#
# Implementation is based on Lemma 5, 6, 7, 8 of paper (2)
#
# Detail:
def genDnaWord1To6And9Generic(n, mapTypeToParam, constructBoundedEDnaFunc):
	# Generate a set of words satisfying C1 to C6 constraint
	strSet16 = genDnaWord1To6(n, mapTypeToParam)

	# Create the pairwise energy function from the 4 x 4 matrices
	pairwiseEnergy = free_energy_routine.createPairwiseEnergyFunc(mapTypeToParam[9])
	maxPairwise = free_energy_routine.findExtremePairwiseEnergy(pairwiseEnergy, max)
	minPairwise = free_energy_routine.findExtremePairwiseEnergy(pairwiseEnergy, min)
	D = maxPairwise - minPairwise

	# Compute the free energy of each string
	freeEnergyList = [free_energy_routine.computeFreeEnergy(strSet16[i], pairwiseEnergy) for i in xrange(n)]
	feMax = max(freeEnergyList)
	feMin = min(freeEnergyList)

	if feMax - feMin <= 3 * D:
		return strSet16

	newLen = len(strSet16[0])
	newLen += newLen

	alpha = feMax + free_energy_routine.computeMinFreeEnergy(newLen, pairwiseEnergy)
	beta = alpha + D + D
	tripleList =  [(newLen, alpha - freeEnergyList[i], beta - freeEnergyList[i]) for i in xrange(n)]
	
	# Find the list of bounded energy DNA strings
	boundedEStrList = constructBoundedEDnaFunc(tripleList, pairwiseEnergy)

	strSet = []
	halfLen = newLen / 2
	for ind in xrange(n):
		newStr = boundedEStrList[ind][0 : halfLen] + strSet16[ind] + boundedEStrList[ind][halfLen : newLen]
		strSet.append(newStr)

	return strSet 
