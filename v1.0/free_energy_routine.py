# Author: Nguyen Truong Duy
# Email: truongduy134@gmail.com

# Self-written modules
import m_poly

# Global variables
typeArr = ['A', 'C', 'G', 'T']
numType = 4

#############################################################################
# computeFreeEnergy(strA, pairwiseEnergy):
#	- Returns the approximated free energy of the input DNA string.
#
# pairwiseEnergy is the input function in which
#	pairwiseEnergy(X, Y) = some real number where X, Y = {A, G, T, C}
#
# The way we compute the approximated free energy of a DNA string X is specified
#	Section 4 of the paper.
# freeEnergy(X) = sum pairwiseEnergy(X(i), X(i + 1)) for i = 0, ..., L - 2 
def computeFreeEnergy(strA, pairwiseEnergy):
	freeEnergy = 0
	lenA = len(strA)
	for ind in xrange(lenA - 1):
		freeEnergy += pairwiseEnergy(strA[ind], strA[ind + 1])
	return freeEnergy

##############################################################################
# findExtremePairwiseEnergy(pairwiseEnergy, extremeFunc):
#	+ Finds an extreme value (min / max) among all 16 possible entries
#		returned from pairwiseEnergy
#
# Input specification:
#	+ extremeFunc: a binary function min or max
#	+ pairwiseEnergy: a function that takes in two inputs X and Y where
#		X, Y in {'A', 'C', 'G', 'T'} and returns the free energy
#		between the ordered pair (X, Y), which is an integer value. 
def findExtremePairwiseEnergy(pairwiseEnergy, extremeFunc):
	global typeArr
	global numType

	extremeVal = pairwiseEnergy('A', 'A')
	for typeOne in xrange(numType):
		for typeTwo in xrange(numType):
			extremeVal = extremeFunc(extremeVal, pairwiseEnergy(typeArr[typeOne], typeArr[typeTwo]))

	return extremeVal

################################################################################
# createPairwiseEnergyFunc(pairwiseEnergyMat, indArr = ['A', 'C', 'G', 'T']):
#	+ Creates a pairwise free energy function from the input 
def createPairwiseEnergyFunc(pairwiseEnergyMat, indArr = ['A', 'C', 'G', 'T']):
	global numType

	mapTypeToInd = {}
	for ind in xrange(numType):
		mapTypeToInd[indArr[ind]] = ind
	
	def pairwiseFunc(typeOne, typeTwo):
		return pairwiseEnergyMat[mapTypeToInd[typeOne]][mapTypeToInd[typeTwo]]
	return pairwiseFunc

#################################################################################
# computeMinFreeEnergy(l, pairwiseEnergy):
#	+ Returns the minimum free energy among all possible values of all DNA
#	  strings of length l with respect to the given pairwiseEnergy function.
#
# Input specification:
#	+ pairwiseEnergy: a function that takes in two inputs X and Y where
#		X, Y in {'A', 'C', 'G', 'T'} and returns the free energy
#		between the ordered pair (X, Y), which is an integer value.
#
# Methodology:
#	+ We use dynamic programming to compute the value.
# 	+ Define:
#		f(l, c) = the minimum free energy among all possible values of all
#		DNA strings of length l and the first character is c
#	+ Then the value we want to compute is:
#		min{f(l, 'A'), f(l, 'C'), f(l, 'G'), f(l, 'T')}
#
#	+ Recurrence relation:
#		f(l, c) = min{pairwiseEnergy(c, d) + f(l - 1, d)} where d in {'A', 'C', 'G', 'T'}
#	+ Base cases:
#		f(1, c) = 0 for all c
#		f(0, c) = 0 for all c
#
#	+ We implment computeMinFreeEnergyDP with recursion and memoization to compute f
#
# Time complexity: O(l)
#
# Note on the implementation:
#	+ We use integers to represent the type of nucleotide, i.e. i represents typeArr[i]
#	  (typeArr is a global variable) for i = 0, 1, 2, 3
def computeMinFreeEnergy(l, pairwiseEnergy):
	global numType
	memoTable = [[-1] * l for i in xrange(numType)]

	minEnergy = computeMinFreeEnergyDP(l, 0, pairwiseEnergy, memoTable)
	for typeId in xrange(1, numType):
		minEnergy = min(minEnergy, computeMinFreeEnergyDP(l, typeId, pairwiseEnergy, memoTable))

	return minEnergy

#
# computeMinFreeEnergyDP(l, typeId, pairwiseEnergy, memoTable):
#
#	Refer to the documentation of computeMinFreeEnergy
#
def computeMinFreeEnergyDP(l, typeId, pairwiseEnergy, memoTable):
	global typeArr
	global numType

	if l <= 1:
		return 0
	if memoTable[typeId][l - 1] >= 0:
		return memoTable[typeId][l - 1]

	
	minEnergy = -1
	for nextTypeId in xrange(numType):
		energy = pairwiseEnergy(typeArr[typeId], typeArr[nextTypeId]) + computeMinFreeEnergyDP(l - 1, nextTypeId, pairwiseEnergy, memoTable)
		if minEnergy < 0 or energy < minEnergy:
			minEnergy = energy

	memoTable[typeId][l - 1] = minEnergy
	return minEnergy

###############################################################################################
# computeExistFreeEnergyDnaMat(maxLen, pairwiseEnergy, maxEnergy = -1):
#	+ Computes and returns a binary matrix M of dimension numType x maxLen x maxE 
#	  (where numType = 4) such that
#		M[c][L][E] = 1 means there exists a DNA string of length (L + 1) that 
#			starts with the character c and has free energy E
#			with respect to the given pairwiseEnergy function.
#		M[c][L][E] = 0 otherwise 
#
# Input specification:
#	+ maxLen: an integer indicating the maximum length of a DNA string
#	+ pairwiseEnergy: a function that takes in two inputs X and Y where
#		X, Y in {'A', 'C', 'G', 'T'} and returns the free energy
#		between the ordered pair (X, Y), which is an integer value.
#	+ maxEnergy: an integer indicating the maximum energy we are considering.
#		maxEnergy = -1 if it is not specified.
#
# This is how we determine maxE (Note that our resulted matrix is of size maxLen x maxE):
#	+ Let W be maximum value among 16 values of pairwiseEnergy
#	+ Let maxPossibleE = W * (maxLen - 1). maxPossibleE is the maximum free energy possible
#		for all DNA strings of length maxLen
#	+ If maxEnergy > 0 and maxEnergy < maxPossibleE, then maxE = maxEnergy
#	+ Otherwise, maxE = maxPossibleE
#
# Methodology to compute the matrix:
#	+ We use dynamic programming to solve this problem.
#	+ Define the function
#		f(c, l, E) = 1 if there exists a DNA string of length l that 
#			starts with the character c and has free energy E
#		f(c, l, E) = 0 otherwise
#	+ We see that M[c][L][E] = f(c, L + 1, E)
#	+ Recurrence relation:
#		f(c, l, E) = OR{f(d, l - 1, E - pairwiseEnergy(c, d)} for d in {'A', 'C', 'G', 'T'}
#	+ Base cases:
#		f(c, l, 0) = 1 if l <= 1 for all c
#		f(c, l, E) = 0 if l <= 1 and E > 0
#		f(c, l, E) = 0 if E < 0
#
#	+ We implment existFreeEnergyDnaDP with recursion and memoization to compute f
#
# Time complexity: O(L^2) (Note that maxE = O(L)
#
# Note on the implementation:
#	+ We use integers to represent the type of nucleotide, i.e. i represents typeArr[i]
#	  (typeArr is a global variable) for i = 0, 1, 2, 3
#
# Important note: M[c][L][E] = f(c, L + 1, E)
#
def computeExistFreeEnergyDnaMat(maxLen, pairwiseEnergy, maxEnergy = -1):
	global numType
	
	# Find maximum possible energy which is equal to
	#	(maxLen - 1) * maxPairwiseEnergy
	# where maxPairwiseEnergy = max(pairwiseEnergy(X, Y)) (X, Y can be 'A', 'C', 'G', 'T')
	maxPairwiseEnergy = findExtremePairwiseEnergy(pairwiseEnergy, max)
	maxPossibleEnergy = (maxLen - 1) * maxPairwiseEnergy
	if maxEnergy < 0 or maxEnergy > maxPossibleEnergy:
		maxEnergy = maxPossibleEnergy

	# Declare and initialize existence matrix
	existMat = [[[-1] * (maxEnergy + 1) for l in xrange(maxLen)] for typeId in xrange(numType)]

	# Repeatedly call the recursive routine (with memoization) to compute all entries of 
	# the existence matrix
	for typeId in xrange(numType):
		for l in xrange(maxLen):
			for energy in xrange(maxEnergy + 1):
				if existMat[typeId][l][energy] < 0:
					existMat[typeId][l][energy] = existFreeEnergyDnaDP(l + 1, typeId, energy, pairwiseEnergy, existMat)

	return existMat

#
# existFreeEnergyDnaDP(l, typeId, freeEnergy, pairwiseEnergy, memoTable):
#
# 	Refer to the documentation of computeExistFreeEnergyDnaMat
#
def existFreeEnergyDnaDP(l, typeId, freeEnergy, pairwiseEnergy, memoTable):
	global typeArr
	global numType

	if freeEnergy < 0:
		return 0
	if l <= 1:
		if freeEnergy == 0:
			return 1
		return 0

	if memoTable[typeId][l - 1][freeEnergy] >= 0:
		return memoTable[typeId][l - 1][freeEnergy]

	existFlg = 0
	for nextTypeId in xrange(numType):
		existFlg = existFlg or existFreeEnergyDnaDP(l - 1, nextTypeId, freeEnergy - pairwiseEnergy(typeArr[typeId], typeArr[nextTypeId]), pairwiseEnergy, memoTable)
		if existFlg:
			break

	memoTable[typeId][l - 1][freeEnergy] = existFlg
	return existFlg

##############################################################################################
# constructBoundedEnergyDnaList(tripleList, pairwiseEnergy):
# 	+ Constructs and returns the list of DNA strings with length and bounded energy 
#	  satisfied the conditions specified in tripleList. 
#
# Input specification:
#	+ tripleList: a list of triples (l, A, B). Each triple (l, A, B) contains information
#		about the DNA we want to construct:
#			- l: the length of the DNA string
#			- A: the lower bound of the free energy of the DNA string
#			- B: the upper bound of the free energy of the DNA string
#	+ pairwiseEnergy: a function that takes in two inputs X and Y where
#		X, Y in {'A', 'C', 'G', 'T'} and returns the free energy
#		between the ordered pair (X, Y), which is an integer value.
#
# Detailed Flow:
#	+ Let L be the maximum length of all DNA strings we want to construct
#	  Let E be the maximum possible free energy 
#	+ Construct the existence matrix M of dimension numType x L x E (using the function
#		computeExistFreeEnergyDnaMat) (numType = 4 typically) such that:
#			M[c][l][e] = 1 means there exists a DNA string of length (l + 1) that 
#				starts with the character c and has free energy e
#				with respect to the given pairwiseEnergy function.
#			M[c][l][e] = 0 otherwise
#	+ For each triple (l, A, B), we search if there exists some energy e BETWEEN A and B
#		such that M[c][l][e] = 1 for some type c
#	  If such e exists, we call the function constructFreeEnergyDna to construct
#		the DNA string of length l, starting with character c and has free energy e
#		which is bounded between A and B
#
# Time complexity: O(L^2 + n * L) since (E = O(L)) where
#	- O(L^2) is due to constructing the existence matrix
#	- O(n * L) is due to constructing n bounded free energy DNA strings of length at most L
#
def constructBoundedEnergyDnaList(tripleList, pairwiseEnergy):
	global numType

	numDna = len(tripleList)
	if numDna == 0:
		return []

	maxLen = max([tripleList[i][0] for i in xrange(numDna)])
	maxEnergy = max([tripleList[i][2] for i in xrange(numDna)])
	maxEnergy = max(maxEnergy, max([tripleList[i][1] for i in xrange(numDna)]))

	existMat = computeExistFreeEnergyDnaMat(maxLen, pairwiseEnergy, maxEnergy)
	maxPossibleE = len(existMat[0][0]) - 1

	boundedEStrSet = []

	# Assumption: l >= 1, 0 <= energy <= maxPossibleE
	def existEnergyDna(l, energy):
		for typeId in xrange(numType):
			if existMat[typeId][l - 1][energy]:
				return True
		return False

	for ind in xrange(numDna):
		# Preprocess / Standardize input data
		minE = tripleList[ind][1]
		if minE < 0:
			minE = 0
		maxE = tripleList[ind][2]
		if maxE > maxPossibleE:
			maxE = maxPossibleE
		if minE > maxE:
			temp = minE
			minE = maxE
			maxE = temp
		l = tripleList[ind][0]

		if l < 1:
			boundedEStrSet.append('')
			continue

		found = False
		for energy in xrange(minE, maxE + 1):
			if existEnergyDna(l, energy):
				found = True
				boundedEStrSet.append(constructFreeEnergyDna(l, energy, pairwiseEnergy, existMat))
				break

		if not found:
			boundedEStrSet.append('')

	return boundedEStrSet		
			
	
##############################################################################################
# constructFreeEnergyDna(l, energy, pairwiseEnergy, existMat):
#	+ Constructs and returns a DNA string of length l and having free energy as specified
#		in the second parameter with respect to the input pairwiseEnergy function
#
# Input specfication:
#	+ l: an integer indicating the length of the DNA string to be constructed
#	+ energy: the free energy of the DNA string to be constructed
#	+ pairwiseEnergy: a function that takes in two inputs X and Y where
#		X, Y in {'A', 'C', 'G', 'T'} and returns the free energy
#		between the ordered pair (X, Y), which is an integer value.
#	+ existMat: a 3D array of dimension numType x L x E (using the function
#		computeExistFreeEnergyDnaMat) (numType = 4 typically) such that:
#			- L >= l
#			- E > energy
#			- M[c][l][e] = 1 means there exists a DNA string of length (l + 1) that 
#				starts with the character c and has free energy e
#				with respect to the given pairwiseEnergy function.
#			- M[c][l][e] = 0 otherwise
#
# Detail:
#	+ Clearly, if l < 1, then we return the empty string
#	+ For the first character c1 of the DNA, we choose any c1 among {'A', 'C', 'G', 'T'} such that
#		existMat[c1][l - 1][energy] = 1
#	+ For the second character c2 of the DNA, our problem now becomes finding a DNA string of length
#		l - 1, starting with c2 and have energy equals to energy - pairwiseEnergy(c1, c2). This #		means we choose any c2 among {'A', 'C', 'G', 'T'} such that
#			existMat[c2][l - 2][energy - pairwiseEnergy(c1, c2)] = 1
#	+ Continue in this manner until we find all n characters
#
# Time complexity: O(l)
def constructFreeEnergyDna(l, energy, pairwiseEnergy, existMat):
	global numType
	global typeArr

	if l < 1:
		return ''

	charArr = []
	curLen = 0
	remainL = l
	while remainL > 0:
		exist = False
		
		for typeId in xrange(numType):
			remainEnergy = energy
			if curLen > 0:
				remainEnergy -= pairwiseEnergy(charArr[curLen - 1], typeArr[typeId])

			if existMat[typeId][remainL - 1][remainEnergy]:
				charArr.append(typeArr[typeId])
				energy = remainEnergy
				curLen += 1
				remainL -= 1
				exist = True
				break

		if not exist:
			return ''

	return ''.join(charArr)

############################################################################################3
#
def buildFreeEnergyPoly(length, pairwiseEnergy):
	global typeArr

	polyPairList = []

	if length == 0:
		return polyPairList

	if length == 1:
		dictX = {}
		dictY = {}
		for typeA in typeArr:
			for typeB in typeArr:
				if typeA == typeB:
					dictX[(length, typeA, typeB)] = m_poly.MyPoly([1])
				else:
					dictX[(length, typeA, typeB)] = m_poly.MyPoly()
		for typeA in typeArr:
			for typeB in typeArr:
				for d1 in typeArr:
					for d2 in typeArr:
						for d3 in typeArr:
							dictY[(length, typeA, typeB, d1, d2, d3)] = None
		polyPairList.append((dictX, dictY))

		return polyPairList

	# Recurse
	polyPairList = buildFreeEnergyPoly(length >> 1, pairwiseEnergy)

	if (length & 1) == 0:
		# length is even
		dictY = {}
		for typeA in typeArr:
			for typeB in typeArr:
				for d1 in typeArr:
					for d2 in typeArr:
						tempPoly = polyPairList[-1][0][(length >> 1, typeA, d1)] * polyPairList[-1][0][(length >> 1, d2, typeB)]
						tempPoly = tempPoly.multiplyBySingleton(pairwiseEnergy(d1, d2), 1)
						dictY[(length, typeA, typeB, d1, d2)] = tempPoly

		
		dictX = {}
		for typeA in typeArr:
			for typeB in typeArr:
				poly = m_poly.MyPoly()
				for d1 in typeArr:
					for d2 in typeArr:
						poly += dictY[(length, typeA, typeB, d1, d2)]
				dictX[(length, typeA, typeB)] = poly

		polyPairList.append((dictX, dictY))
	else:
		# length is odd
		dictY = {}
		for typeA in typeArr:
			for typeB in typeArr:
				for d1 in typeArr:
					for d2 in typeArr:
						for d3 in typeArr:
							tempPoly = polyPairList[-1][0][(length >> 1, typeA, d1)] * polyPairList[-1][0][(length >> 1, d2, typeB)]
							tempPoly = tempPoly.multiplyBySingleton(pairwiseEnergy(d1, d3) + pairwiseEnergy(d3, d2), 1)
							dictY[(length, typeA, typeB, d1, d2, d3)] = tempPoly

		dictX = {}
		for typeA in typeArr:
			for typeB in typeArr:
				poly = m_poly.MyPoly()
				for d1 in typeArr:
					for d2 in typeArr:
						for d3 in typeArr:
							poly += dictY[(length, typeA, typeB, d1, d2, d3)]
				dictX[(length, typeA, typeB)] = poly
		polyPairList.append((dictX, dictY))	

	return polyPairList

####################################################################################
#
def extractFreeEnergyDnaPoly(length, energy, pairwiseEnergy, polyPairList):
	global typeArr
	
	lastInd = len(polyPairList) - 1
	for typeA in typeArr:
		for typeB in typeArr:
			if hasPositiveEnergyCoeff(polyPairList[-1][0][(length, typeA, typeB)], energy):
				charArr = recursiveExtract(length, energy, typeA, typeB, pairwiseEnergy, polyPairList, lastInd)
				return ''.join(charArr)

	return ''	# Default: Return empty string

################################################################################################
#
def recursiveExtract(length, energy, typeA, typeB, pairwiseEnergy, polyPairList, lengthPolyInd):
	global typeArr

	# Base cases:
	if length == 2:
		return [typeA, typeB]
	if length == 1:
		return [typeA]

	halfLength = (length >> 1)
	
	paramTuple = None
	#########################################################
	# CHECK THIS PART AGAIN
	if (length & 1) == 0:
		# length is even
		# Find triple (z, d1, d2)
		paramTuple = findD1D2(halfLength, energy, typeA, typeB, pairwiseEnergy, polyPairList[lengthPolyInd - 1][0])
	else:
		# length is odd
		# Find 4-tuple (z, d1, d2, d3)
		paramTuple = findD1D2D3(halfLength, energy, typeA, typeB, pairwiseEnergy, polyPairList[lengthPolyInd - 1][0])
	########################################################3
	if paramTuple is None:
		return []

	charArr = recursiveExtract(halfLength, paramTuple[0], typeA, paramTuple[1], pairwiseEnergy, polyPairList, lengthPolyInd - 1)

	if charArr == []:
		return []
	
	subE = 0
	if (length & 1) == 0:
		# length is even
		subE = energy - paramTuple[0] - pairwiseEnergy(paramTuple[1], paramTuple[2]) 
	else:
		# length is odd
		subE = energy - paramTuple[0] - pairwiseEnergy(paramTuple[1], paramTuple[3]) - pairwiseEnergy(paramTuple[3], paramTuple[2])
		charArr.append(paramTuple[3])
	charArr.extend(recursiveExtract(halfLength, subE, paramTuple[2], typeB, pairwiseEnergy, polyPairList, lengthPolyInd - 1))

	return charArr

####################################################################################
#
def hasPositiveEnergyCoeff(poly, energy):
	if energy < 0:
		return False
	if poly.degree < energy:
		return False

	return poly.coeff[energy] > 0

#######################################################################################
#
def findD1D2(length, energy, typeA, typeB, pairwiseEnergy, polySet):
	global typeArr

	for d1 in typeArr:
		polyD1 =  polySet[(length, typeA, d1)]
		maxE = polyD1.degree
		for subE in xrange(maxE + 1):
			if polyD1.coeff[subE] < 1:
				continue
			if energy - subE < 0:
				break

			for d2 in typeArr:
				leftE = energy - subE - pairwiseEnergy(d1, d2)
				if hasPositiveEnergyCoeff(polySet[(length, d2, typeB)], leftE):
					return (subE, d1, d2)
		
	return None

#######################################################################################
#
def findD1D2D3(length, energy, typeA, typeB, pairwiseEnergy, polySet):
	global typeArr

	for d1 in typeArr:
		polyD1 = polySet[(length, typeA, d1)]
		maxE = polyD1.degree

		for subE in xrange(maxE + 1):
			if polyD1.coeff[subE] < 1:
				continue
			if energy - subE < 0:
				break

			for d2 in typeArr:
				for d3 in typeArr:
					leftE = energy - subE - pairwiseEnergy(d1, d3) - pairwiseEnergy(d3, d2)
					if hasPositiveEnergyCoeff(polySet[(length, d2, typeB)], leftE):
						return (subE, d1, d2, d3)

	return None

##################################################################################################
# constructBoundedEnergyDnaListPoly(tripleList, pairwiseEnergy):
# 	+ Constructs and returns the list of DNA strings with the SAME length and bounded energy 
#	  satisfied the conditions specified in tripleList. 
#
# Input specification:
#	+ tripleList: a list of triples (l, A, B). Each triple (l, A, B) contains information
#		about the DNA we want to construct:
#			- l: the length of the DNA string
#			- A: the lower bound of the free energy of the DNA string
#			- B: the upper bound of the free energy of the DNA string
#	+ pairwiseEnergy: a function that takes in two inputs X and Y where
#		X, Y in {'A', 'C', 'G', 'T'} and returns the free energy
#		between the ordered pair (X, Y), which is an integer value.
#
# Assumption:
# 	+ Unlike constructBoundedEnergyDnaList, in this function, we assume all
#		triples must have the same length information.
#
# Detailed Flow:
#	+ 
def constructBoundedEnergyDnaListPoly(tripleList, pairwiseEnergy):
	numDna = len(tripleList)
	if numDna == 0:
		return []
	if max([tripleList[i][0] for i in xrange(numDna)]) != min([tripleList[i][0] for i in xrange(numDna)]):
		raise RuntimeError("This function only generates bounded DNA strings of EQUAL length")

	length = tripleList[0][0]
	
	polyPairList = buildFreeEnergyPoly(length, pairwiseEnergy)

	dnaList = []
	for i in xrange(numDna):
		lowerE = tripleList[i][1]
		upperE = tripleList[i][2]
		if lowerE > upperE:
			lowerE, upperE = upperE, lowerE

		dnaStr = ''
		for energy in xrange(lowerE, upperE + 1):
			dnaStr = extractFreeEnergyDnaPoly(length, energy, pairwiseEnergy, polyPairList)
			if dnaStr != '':
				break
		dnaList.append(dnaStr)

	return dnaList
