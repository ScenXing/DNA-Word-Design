# Author: Nguyen Truong Duy
# Email: truongduy134@gmail.com

# Self-written modules
import m_poly

# Global variables
TYPE_ARR = ['A', 'C', 'G', 'T']
NUM_TYPE = 4

#############################################################################
# compute_free_energy(strA, pairwiseenergy):
#	- Returns the approximated free energy of the input DNA string.
#
# pairwiseenergy is the input function in which
#	pairwiseenergy(X, Y) = some real number where X, Y = {A, G, T, C}
#
# The way we compute the approximated free energy of a DNA string X is specified
#	Section 4 of the paper.
# freeenergy(X) = sum pairwiseenergy(X(i), X(i + 1)) for i = 0, ..., L - 2 
def compute_free_energy(strA, pairwiseenergy):
	"""
	"""

	freeenergy = 0
	lenA = len(strA)
	for ind in xrange(lenA - 1):
		freeenergy += pairwiseenergy(strA[ind], strA[ind + 1])
	return freeenergy

##############################################################################
# find_extreme_pairwise_energy(pairwiseenergy, extremeFunc):
#	+ Finds an extreme value (min / max) among all 16 possible entries
#		returned from pairwiseenergy
#
# Input specification:
#	+ extremeFunc: a binary function min or max
#	+ pairwiseenergy: a function that takes in two inputs X and Y where
#		X, Y in {'A', 'C', 'G', 'T'} and returns the free energy
#		between the ordered pair (X, Y), which is an integer value. 
def find_extreme_pairwise_energy(pairwiseenergy, extremeFunc):
	global TYPE_ARR
	global NUM_TYPE

	extremeVal = pairwiseenergy('A', 'A')
	for typeone in xrange(NUM_TYPE):
		for typetwo in xrange(NUM_TYPE):
			extremeVal = extremeFunc(extremeVal, pairwiseenergy(TYPE_ARR[typeone], TYPE_ARR[typetwo]))

	return extremeVal

################################################################################
# create_pairwise_energy_func(pairwiseenergyMat, indArr = ['A', 'C', 'G', 'T']):
#	+ Creates a pairwise free energy function from the input 
def create_pairwise_energy_func(pairwiseenergyMat, indArr = ['A', 'C', 'G', 'T']):
	global NUM_TYPE

	mapTypeToInd = {}
	for ind in xrange(NUM_TYPE):
		mapTypeToInd[indArr[ind]] = ind
	
	def pairwiseFunc(typeone, typetwo):
		return pairwiseenergyMat[mapTypeToInd[typeone]][mapTypeToInd[typetwo]]
	return pairwiseFunc

#################################################################################
# compute_min_free_energy(l, pairwiseenergy):
#	+ Returns the minimum free energy among all possible values of all DNA
#	  strings of length l with respect to the given pairwiseenergy function.
#
# Input specification:
#	+ pairwiseenergy: a function that takes in two inputs X and Y where
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
#		f(l, c) = min{pairwiseenergy(c, d) + f(l - 1, d)} where d in {'A', 'C', 'G', 'T'}
#	+ Base cases:
#		f(1, c) = 0 for all c
#		f(0, c) = 0 for all c
#
#	+ We implment compute_min_free_energy_dp with recursion and memoization to compute f
#
# Time complexity: O(l)
#
# Note on the implementation:
#	+ We use integers to represent the type of nucleotide, i.e. i represents TYPE_ARR[i]
#	  (TYPE_ARR is a global variable) for i = 0, 1, 2, 3
def compute_min_free_energy(l, pairwiseenergy):
	global NUM_TYPE
	memotable = [[-1] * l for i in xrange(NUM_TYPE)]

	minEnergy = compute_min_free_energy_dp(l, 0, pairwiseenergy, memotable)
	for typeId in xrange(1, NUM_TYPE):
		minEnergy = min(minEnergy, compute_min_free_energy_dp(l, typeId, pairwiseenergy, memotable))

	return minEnergy

#
# compute_min_free_energy_dp(l, typeId, pairwiseenergy, memotable):
#
#	Refer to the documentation of compute_min_free_energy
#
def compute_min_free_energy_dp(l, typeId, pairwiseenergy, memotable):
	global TYPE_ARR
	global NUM_TYPE

	if l <= 1:
		return 0
	if memotable[typeId][l - 1] >= 0:
		return memotable[typeId][l - 1]

	
	minEnergy = -1
	for nextTypeId in xrange(NUM_TYPE):
		energy = pairwiseenergy(TYPE_ARR[typeId], TYPE_ARR[nextTypeId]) + compute_min_free_energy_dp(l - 1, nextTypeId, pairwiseenergy, memotable)
		if minEnergy < 0 or energy < minEnergy:
			minEnergy = energy

	memotable[typeId][l - 1] = minEnergy
	return minEnergy

###############################################################################################
# compute_exist_free_energy_dna_mat(maxLen, pairwiseenergy, maxEnergy = -1):
#	+ Computes and returns a binary matrix M of dimension NUM_TYPE x maxLen x maxE 
#	  (where NUM_TYPE = 4) such that
#		M[c][L][E] = 1 means there exists a DNA string of length (L + 1) that 
#			starts with the character c and has free energy E
#			with respect to the given pairwiseenergy function.
#		M[c][L][E] = 0 otherwise 
#
# Input specification:
#	+ maxLen: an integer indicating the maximum length of a DNA string
#	+ pairwiseenergy: a function that takes in two inputs X and Y where
#		X, Y in {'A', 'C', 'G', 'T'} and returns the free energy
#		between the ordered pair (X, Y), which is an integer value.
#	+ maxEnergy: an integer indicating the maximum energy we are considering.
#		maxEnergy = -1 if it is not specified.
#
# This is how we determine maxE (Note that our resulted matrix is of size maxLen x maxE):
#	+ Let W be maximum value among 16 values of pairwiseenergy
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
#		f(c, l, E) = OR{f(d, l - 1, E - pairwiseenergy(c, d)} for d in {'A', 'C', 'G', 'T'}
#	+ Base cases:
#		f(c, l, 0) = 1 if l <= 1 for all c
#		f(c, l, E) = 0 if l <= 1 and E > 0
#		f(c, l, E) = 0 if E < 0
#
#	+ We implment exist_free_energy_dna_dp with recursion and memoization to compute f
#
# Time complexity: O(L^2) (Note that maxE = O(L)
#
# Note on the implementation:
#	+ We use integers to represent the type of nucleotide, i.e. i represents TYPE_ARR[i]
#	  (TYPE_ARR is a global variable) for i = 0, 1, 2, 3
#
# Important note: M[c][L][E] = f(c, L + 1, E)
#
def compute_exist_free_energy_dna_mat(maxLen, pairwiseenergy, maxEnergy = -1):
	global NUM_TYPE
	
	# Find maximum possible energy which is equal to
	#	(maxLen - 1) * maxPairwiseEnergy
	# where maxPairwiseEnergy = max(pairwiseenergy(X, Y)) (X, Y can be 'A', 'C', 'G', 'T')
	maxPairwiseEnergy = find_extreme_pairwise_energy(pairwiseenergy, max)
	maxPossibleEnergy = (maxLen - 1) * maxPairwiseEnergy
	if maxEnergy < 0 or maxEnergy > maxPossibleEnergy:
		maxEnergy = maxPossibleEnergy

	# Declare and initialize existence matrix
	existmat = [[[-1] * (maxEnergy + 1) for l in xrange(maxLen)] for typeId in xrange(NUM_TYPE)]

	# Repeatedly call the recursive routine (with memoization) to compute all entries of 
	# the existence matrix
	for typeId in xrange(NUM_TYPE):
		for l in xrange(maxLen):
			for energy in xrange(maxEnergy + 1):
				if existmat[typeId][l][energy] < 0:
					existmat[typeId][l][energy] = exist_free_energy_dna_dp(l + 1, typeId, energy, pairwiseenergy, existmat)

	return existmat

#
# exist_free_energy_dna_dp(l, typeId, freeenergy, pairwiseenergy, memotable):
#
# 	Refer to the documentation of compute_exist_free_energy_dna_mat
#
def exist_free_energy_dna_dp(l, typeId, freeenergy, pairwiseenergy, memotable):
	global TYPE_ARR
	global NUM_TYPE

	if freeenergy < 0:
		return 0
	if l <= 1:
		if freeenergy == 0:
			return 1
		return 0

	if memotable[typeId][l - 1][freeenergy] >= 0:
		return memotable[typeId][l - 1][freeenergy]

	existFlg = 0
	for nextTypeId in xrange(NUM_TYPE):
		existFlg = existFlg or exist_free_energy_dna_dp(l - 1, nextTypeId, freeenergy - pairwiseenergy(TYPE_ARR[typeId], TYPE_ARR[nextTypeId]), pairwiseenergy, memotable)
		if existFlg:
			break

	memotable[typeId][l - 1][freeenergy] = existFlg
	return existFlg

##############################################################################################
# construct_bounded_energy_dna_list(triplelist, pairwiseenergy):
# 	+ Constructs and returns the list of DNA strings with length and bounded energy 
#	  satisfied the conditions specified in triplelist. 
#
# Input specification:
#	+ triplelist: a list of triples (l, A, B). Each triple (l, A, B) contains information
#		about the DNA we want to construct:
#			- l: the length of the DNA string
#			- A: the lower bound of the free energy of the DNA string
#			- B: the upper bound of the free energy of the DNA string
#	+ pairwiseenergy: a function that takes in two inputs X and Y where
#		X, Y in {'A', 'C', 'G', 'T'} and returns the free energy
#		between the ordered pair (X, Y), which is an integer value.
#
# Detailed Flow:
#	+ Let L be the maximum length of all DNA strings we want to construct
#	  Let E be the maximum possible free energy 
#	+ Construct the existence matrix M of dimension NUM_TYPE x L x E (using the function
#		compute_exist_free_energy_dna_mat) (NUM_TYPE = 4 typically) such that:
#			M[c][l][e] = 1 means there exists a DNA string of length (l + 1) that 
#				starts with the character c and has free energy e
#				with respect to the given pairwiseenergy function.
#			M[c][l][e] = 0 otherwise
#	+ For each triple (l, A, B), we search if there exists some energy e BETWEEN A and B
#		such that M[c][l][e] = 1 for some type c
#	  If such e exists, we call the function construct_free_energy_dna to construct
#		the DNA string of length l, starting with character c and has free energy e
#		which is bounded between A and B
#
# Time complexity: O(L^2 + n * L) since (E = O(L)) where
#	- O(L^2) is due to constructing the existence matrix
#	- O(n * L) is due to constructing n bounded free energy DNA strings of length at most L
#
def construct_bounded_energy_dna_list(triplelist, pairwiseenergy):
	global NUM_TYPE

	numdna = len(triplelist)
	if numdna == 0:
		return []

	maxLen = max([triplelist[i][0] for i in xrange(numdna)])
	maxEnergy = max([triplelist[i][2] for i in xrange(numdna)])
	maxEnergy = max(maxEnergy, max([triplelist[i][1] for i in xrange(numdna)]))

	existmat = compute_exist_free_energy_dna_mat(maxLen, pairwiseenergy, maxEnergy)
	maxPossibleE = len(existmat[0][0]) - 1

	boundedEStrSet = []

	# Assumption: l >= 1, 0 <= energy <= maxPossibleE
	def existEnergyDna(l, energy):
		for typeId in xrange(NUM_TYPE):
			if existmat[typeId][l - 1][energy]:
				return True
		return False

	for ind in xrange(numdna):
		# Preprocess / Standardize input data
		minE = triplelist[ind][1]
		if minE < 0:
			minE = 0
		maxE = triplelist[ind][2]
		if maxE > maxPossibleE:
			maxE = maxPossibleE
		if minE > maxE:
			temp = minE
			minE = maxE
			maxE = temp
		l = triplelist[ind][0]

		if l < 1:
			boundedEStrSet.append('')
			continue

		found = False
		for energy in xrange(minE, maxE + 1):
			if existEnergyDna(l, energy):
				found = True
				boundedEStrSet.append(construct_free_energy_dna(l, energy, pairwiseenergy, existmat))
				break

		if not found:
			boundedEStrSet.append('')

	return boundedEStrSet		
			
	
##############################################################################################
# construct_free_energy_dna(l, energy, pairwiseenergy, existmat):
#	+ Constructs and returns a DNA string of length l and having free energy as specified
#		in the second parameter with respect to the input pairwiseenergy function
#
# Input specfication:
#	+ l: an integer indicating the length of the DNA string to be constructed
#	+ energy: the free energy of the DNA string to be constructed
#	+ pairwiseenergy: a function that takes in two inputs X and Y where
#		X, Y in {'A', 'C', 'G', 'T'} and returns the free energy
#		between the ordered pair (X, Y), which is an integer value.
#	+ existmat: a 3D array of dimension NUM_TYPE x L x E (using the function
#		compute_exist_free_energy_dna_mat) (NUM_TYPE = 4 typically) such that:
#			- L >= l
#			- E > energy
#			- M[c][l][e] = 1 means there exists a DNA string of length (l + 1) that 
#				starts with the character c and has free energy e
#				with respect to the given pairwiseenergy function.
#			- M[c][l][e] = 0 otherwise
#
# Detail:
#	+ Clearly, if l < 1, then we return the empty string
#	+ For the first character c1 of the DNA, we choose any c1 among {'A', 'C', 'G', 'T'} such that
#		existmat[c1][l - 1][energy] = 1
#	+ For the second character c2 of the DNA, our problem now becomes finding a DNA string of length
#		l - 1, starting with c2 and have energy equals to energy - pairwiseenergy(c1, c2). This #		means we choose any c2 among {'A', 'C', 'G', 'T'} such that
#			existmat[c2][l - 2][energy - pairwiseenergy(c1, c2)] = 1
#	+ Continue in this manner until we find all n characters
#
# Time complexity: O(l)
def construct_free_energy_dna(l, energy, pairwiseenergy, existmat):
	global NUM_TYPE
	global TYPE_ARR

	if l < 1:
		return ''

	charArr = []
	curLen = 0
	remainL = l
	while remainL > 0:
		exist = False
		
		for typeId in xrange(NUM_TYPE):
			remainEnergy = energy
			if curLen > 0:
				remainEnergy -= pairwiseenergy(charArr[curLen - 1], TYPE_ARR[typeId])

			if existmat[typeId][remainL - 1][remainEnergy]:
				charArr.append(TYPE_ARR[typeId])
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
def build_free_energy_poly(length, pairwiseenergy):
	global TYPE_ARR

	polypairlist = []

	if length == 0:
		return polypairlist

	if length == 1:
		dictX = {}
		dictY = {}
		for typeA in TYPE_ARR:
			for typeB in TYPE_ARR:
				if typeA == typeB:
					dictX[(length, typeA, typeB)] = m_poly.MyPoly([1])
				else:
					dictX[(length, typeA, typeB)] = m_poly.MyPoly()
		for typeA in TYPE_ARR:
			for typeB in TYPE_ARR:
				for d1 in TYPE_ARR:
					for d2 in TYPE_ARR:
						for d3 in TYPE_ARR:
							dictY[(length, typeA, typeB, d1, d2, d3)] = None
		polypairlist.append((dictX, dictY))

		return polypairlist

	# Recurse
	polypairlist = build_free_energy_poly(length >> 1, pairwiseenergy)

	if (length & 1) == 0:
		# length is even
		dictY = {}
		for typeA in TYPE_ARR:
			for typeB in TYPE_ARR:
				for d1 in TYPE_ARR:
					for d2 in TYPE_ARR:
						tempPoly = polypairlist[-1][0][(length >> 1, typeA, d1)] * polypairlist[-1][0][(length >> 1, d2, typeB)]
						tempPoly = tempPoly.multiplyBySingleton(pairwiseenergy(d1, d2), 1)
						dictY[(length, typeA, typeB, d1, d2)] = tempPoly

		
		dictX = {}
		for typeA in TYPE_ARR:
			for typeB in TYPE_ARR:
				poly = m_poly.MyPoly()
				for d1 in TYPE_ARR:
					for d2 in TYPE_ARR:
						poly += dictY[(length, typeA, typeB, d1, d2)]
				dictX[(length, typeA, typeB)] = poly

		polypairlist.append((dictX, dictY))
	else:
		# length is odd
		dictY = {}
		for typeA in TYPE_ARR:
			for typeB in TYPE_ARR:
				for d1 in TYPE_ARR:
					for d2 in TYPE_ARR:
						for d3 in TYPE_ARR:
							tempPoly = polypairlist[-1][0][(length >> 1, typeA, d1)] * polypairlist[-1][0][(length >> 1, d2, typeB)]
							tempPoly = tempPoly.multiplyBySingleton(pairwiseenergy(d1, d3) + pairwiseenergy(d3, d2), 1)
							dictY[(length, typeA, typeB, d1, d2, d3)] = tempPoly

		dictX = {}
		for typeA in TYPE_ARR:
			for typeB in TYPE_ARR:
				poly = m_poly.MyPoly()
				for d1 in TYPE_ARR:
					for d2 in TYPE_ARR:
						for d3 in TYPE_ARR:
							poly += dictY[(length, typeA, typeB, d1, d2, d3)]
				dictX[(length, typeA, typeB)] = poly
		polypairlist.append((dictX, dictY))	

	return polypairlist

####################################################################################
#
def extract_free_energy_dna_poly(length, energy, pairwiseenergy, polypairlist):
	global TYPE_ARR
	
	lastInd = len(polypairlist) - 1
	for typeA in TYPE_ARR:
		for typeB in TYPE_ARR:
			if has_positive_energy_coeff(polypairlist[-1][0][(length, typeA, typeB)], energy):
				charArr = recursive_extract(length, energy, typeA, typeB, pairwiseenergy, polypairlist, lastInd)
				return ''.join(charArr)

	return ''	# Default: Return empty string

################################################################################################
#
def recursive_extract(length, energy, typeA, typeB, pairwiseenergy, polypairlist, lengthPolyInd):
	global TYPE_ARR

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
		paramTuple = findd1d2(halfLength, energy, typeA, typeB, pairwiseenergy, polypairlist[lengthPolyInd - 1][0])
	else:
		# length is odd
		# Find 4-tuple (z, d1, d2, d3)
		paramTuple = findd1d2d3(halfLength, energy, typeA, typeB, pairwiseenergy, polypairlist[lengthPolyInd - 1][0])
	########################################################3
	if paramTuple is None:
		return []

	charArr = recursive_extract(halfLength, paramTuple[0], typeA, paramTuple[1], pairwiseenergy, polypairlist, lengthPolyInd - 1)

	if charArr == []:
		return []
	
	subE = 0
	if (length & 1) == 0:
		# length is even
		subE = energy - paramTuple[0] - pairwiseenergy(paramTuple[1], paramTuple[2]) 
	else:
		# length is odd
		subE = energy - paramTuple[0] - pairwiseenergy(paramTuple[1], paramTuple[3]) - pairwiseenergy(paramTuple[3], paramTuple[2])
		charArr.append(paramTuple[3])
	charArr.extend(recursive_extract(halfLength, subE, paramTuple[2], typeB, pairwiseenergy, polypairlist, lengthPolyInd - 1))

	return charArr

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
def findd1d2(length, energy, typeA, typeB, pairwiseenergy, polySet):
	global TYPE_ARR

	for d1 in TYPE_ARR:
		polyD1 =  polySet[(length, typeA, d1)]
		maxE = polyD1.degree
		for subE in xrange(maxE + 1):
			if polyD1.coeff[subE] < 1:
				continue
			if energy - subE < 0:
				break

			for d2 in TYPE_ARR:
				leftE = energy - subE - pairwiseenergy(d1, d2)
				if has_positive_energy_coeff(polySet[(length, d2, typeB)], leftE):
					return (subE, d1, d2)
		
	return None

#######################################################################################
#
def findd1d2d3(length, energy, typeA, typeB, pairwiseenergy, polySet):
	global TYPE_ARR

	for d1 in TYPE_ARR:
		polyD1 = polySet[(length, typeA, d1)]
		maxE = polyD1.degree

		for subE in xrange(maxE + 1):
			if polyD1.coeff[subE] < 1:
				continue
			if energy - subE < 0:
				break

			for d2 in TYPE_ARR:
				for d3 in TYPE_ARR:
					leftE = energy - subE - pairwiseenergy(d1, d3) - pairwiseenergy(d3, d2)
					if has_positive_energy_coeff(polySet[(length, d2, typeB)], leftE):
						return (subE, d1, d2, d3)

	return None

##################################################################################################
# construct_bounded_energy_dna_list_poly(triplelist, pairwiseenergy):
# 	+ Constructs and returns the list of DNA strings with the SAME length and bounded energy 
#	  satisfied the conditions specified in triplelist. 
#
# Input specification:
#	+ triplelist: a list of triples (l, A, B). Each triple (l, A, B) contains information
#		about the DNA we want to construct:
#			- l: the length of the DNA string
#			- A: the lower bound of the free energy of the DNA string
#			- B: the upper bound of the free energy of the DNA string
#	+ pairwiseenergy: a function that takes in two inputs X and Y where
#		X, Y in {'A', 'C', 'G', 'T'} and returns the free energy
#		between the ordered pair (X, Y), which is an integer value.
#
# Assumption:
# 	+ Unlike construct_bounded_energy_dna_list, in this function, we assume all
#		triples must have the same length information.
#
# Detailed Flow:
#	+ 
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
		lowerE = triplelist[i][1]
		upperE = triplelist[i][2]
		if lowerE > upperE:
			lowerE, upperE = upperE, lowerE

		dnastr = ''
		for energy in xrange(lowerE, upperE + 1):
			dnastr = extract_free_energy_dna_poly(length, energy, pairwiseenergy, polypairlist)
			if dnastr != '':
				break
		dnalist.append(dnastr)

	return dnalist
