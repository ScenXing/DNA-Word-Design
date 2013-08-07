# Author: Nguyen Truong Duy
# Email: truongduy134@gmail.com

#
# Note: All the functions in this module are related to the paper
#	"Deterministic Polynomial-Time Algorithms for Designing
#	 Short DNA Words" by Kao et al.
#
#	You can retrieve a copy of this paper at:
#		http://arxiv.org/pdf/1201.6358.pdf 
#

# Self-written python modules
import helper
import m_fraction

# Standard modules
import math

###############################################################################
# initLength14(n, k1, k4, delta):
#	- Return an integer l indicating the length that satisfies Lemma 6.
#
# The implementation is based on Lemma 7.
#
# Assumption:
#	n >= 2, max{k1, k4} >= 1, and delta is a positive real number
#
# If delta is not specified, the default value 0.1 is used.
def initLength14(n, k1, k4, delta = 0.1):
	c1 = 2 + delta
	c2 = c1 / 2.0 * (math.log(c1 / ((c1 - 2) * math.log(2)), 2) + 2.5 - 1.0 / math.log(2))
	k = max(k1, k4)

	return int(math.ceil(c1 * math.log(n, 2) + c2 * k))

##############################################################################################
# binarySearchLength14(n, k1, k4, initLength, pascalPrefixSum):
#	- Returns an integer which is the minimum length l such that the
#	  distance matrix M of size n x l satisfies Lemma 3, i.e.
#		ExpCount(M, k1, k4) > nC2 * (1 + 2 * (k4 - 1)) - 1
#
# The implementation is based on Lemma 10 and Lemma 11 in the paper. It uses
#	binary search to find the minimum length based on the initial length given
#	as an input.
#
# Input description:
#	- pascalPrefixSum: is a 2D array with at least l + 1 rows. Row i has (i + 1) entries
#	- pascalPrefixSum[i][j] = sum ( iCh ) for h = 0, ..., j
# Note: iCh denotes i Choose h
#
#	- initLength: an integer such that the distance matrix M of size n x initLength
#		satisfies Lemma 3. For example, initLength can be a result returned
#		by the function initLength14 (also implemented in this module)
#
# Time complexity: O(k4 * log(initLength)) 
#	(This does not take into account the complexity of arithmetics)
def binarySearchLength14(n, k1, k4, initLength, pascalPrefixSum):
	k = max(k1, k4)

	######################################################################
	# compute(l):
	#	- Computes and returns the value of ExpCount(M, k1, k4)
	#	  for the EMPTY distance matrix M of size n x l
 	#
	# The formula to compute such value is given in Lemma 10
	#
	# Time complexity: O(k4) due to our pre-computation of some values.
	#	This does not take into account the complexity of arithmetics
	def compute(l):
		result = m_fraction.MyFraction(pascalPrefixSum[l][k - 1], helper.fastPower(2, l))
		
		powTwo = helper.fastPower(2, l - k4 + 1)
		sumTemp = m_fraction.MyFraction(0, 1)
		for i in xrange(l - k4 + 1, l):
			sumTemp += m_fraction.MyFraction(pascalPrefixSum[i][k4 - (l - i) - 1], powTwo)
			powTwo *= 2
 		sumTemp *= m_fraction.MyFraction(2)

		result += sumTemp
		result *= m_fraction.MyFraction((n * (n - 1)) / 2)
		return result

	leftL = 1
	rightL = initLength
	minFound = initLength

	while leftL <= rightL:
		midL = leftL + ((rightL - leftL) >> 1)
		val = compute(midL)
		if val.numer < val.denom:
			minFound = midL
			rightL = midL - 1
		else:
			leftL = midL + 1

	return minFound

#####################################################################################
# computeExpCountEmpty(n, l, k1, k4, pascalPrefixSum):
#	- Computes and returns the value of ExpCount(M, k1, k4) when the matrix
#		M of size n x l is empty. 
#
# The formula is mentioned in the proof of Lemma 10 in the paper.
#
# Input description:
#	- pascalPrefixSum: is a 2D array with at least l + 1 rows. Row i has (i + 1) entries
#	- pascalPrefixSum[i][j] = sum ( iCh ) for h = 0, ..., j
# Note: iCh denotes i Choose h   
def computeExpCountEmpty(n, l, k1, k4, pascalPrefixSum):
	k = max(k1, k4)	
	firstTerm = m_fraction.MyFraction(((n * (n - 1)) / 2) * (1 + 2 * (k4 - 1)), 1)

	secondTerm = m_fraction.MyFraction(pascalPrefixSum[l][k - 1], helper.fastPower(2, l))
		
	powTwo = helper.fastPower(2, l - k4 + 1)
	sumTemp = m_fraction.MyFraction(0, 1)
	for i in xrange(l - k4 + 1, l):
		sumTemp += m_fraction.MyFraction(pascalPrefixSum[i][k4 - (l - i) - 1], powTwo)			
		powTwo *= 2
 	sumTemp *= m_fraction.MyFraction(2)
	secondTerm += sumTemp
	secondTerm *= m_fraction.MyFraction((n * (n - 1)) / 2)

	return firstTerm - secondTerm

#####################################################################################
# breakRun(strA, maxLenRun):
#	- Breaks run in strA to obtain strB such that strB has no more than
#		maxLenRun consecutive bases.
#	  Returns strB
#	- Throws an exception if the input string has odd length or maxLenRun <= 1
#
# Requirement: 
#	- strA must have even length
#	- maxLenRun >= 2
# The implementation is based on Algorithm 2 of the paper
def breakRun(strA, maxLenRun):
	lenA = len(strA)
	if lenA & 1:
		# lenA is odd
		raise RuntimeError("This algorithm requires the input string to have even length")
	if maxLenRun <= 1:
		raise RuntimeError("maxLenRun must be at least 2")

	# Break run in the first half
	leftCharList = []
	lenRun = 0
	for ind in xrange(lenA >> 1):
		leftCharList.append(strA[ind])
		lenRun += 1
		if (lenRun == maxLenRun - 1) or (ind == (lenA >> 1) - 1):
			leftCharList.append(helper.getComplementLetter(strA[ind]))
			lenRun = 0

	# Break run in the second half
	rightCharList = []
	for ind in xrange(lenA - 1, (lenA >> 1) - 1, -1):
		rightCharList.append(strA[ind])
		lenRun += 1
		if (lenRun == maxLenRun - 1) or (ind == (lenA >> 1)):
			rightCharList.append(helper.getComplementLetter(strA[ind]))
			lenRun = 0
	rightCharList.reverse()

	# Return the resulted string
	leftCharList.extend(rightCharList)
	return ''.join(leftCharList)

#########################################################################################
# updateNumDiff(M, strPos, strId, unknownChar, numDiff1, numDiff4, n, l, k4):
#	- Update the difference matrices when the newest fill-in entry is M[strId][strPos]
#		assuming the entries in M are filled in the order top to bottom, 
#		then left to right.
#
# Input description:
#	M: a partially assigned k1, k4-distance matrix of size n x l
#	unknownChar: the character denoting unassigned entries in M
#	strPos, strId: M[strId][strPos] is the newest assigned entry in M
#	k4: the parameter of C4 constraint
#
#	numDiff1: a 2D array that has n rows. Row i has exactly (i + 1) columns (i = 0, ..., n - 1)
#		  numDiff1[a][b] = the number of positions k in strings M[a] and M[b] (b <= a) 
#			such that
#			M[a][k] and M[b][k] are not unknown, and
#			M[a][k] != M[b][k]
#
#	numDiff4: a 3D array of dimension (k4 - 1) x n x n
#		  numDiff4[x][a][b] = the number of positions k in substrings 
#			M[a][0 ... l - k4 + x] and M[b][l - (l - k4 + 1 + x) ... l - 1] such that
#				M[a][k] and M[b][l - (l - k4 + 1 + x) + k] are not unknown, and
#				M[a][k] != M[b][l - (l - k4 + 1 + x) + k]
#
# Time complexity: O(n * k)
def updateNumDiff(M, strPos, strId, unknownChar, numDiff1, numDiff4, n, l, k4):
	# Update numDiff1
	for curStrId in xrange(0, strId):
		if M[curStrId][strPos] != M[strId][strPos]:
			numDiff1[strId][curStrId] += 1

	# Update numDiff4
	for subLen in xrange(l - k4 + 1, l):
		for curStrId in xrange(0, n):
			if curStrId != strId:
				# Y = M[strId] and X = M[curStrId]
				# Compare M[strId][strPos] and M[curStrId][l - subLen + strPos]
				comparePos = l - subLen + strPos
				if comparePos < l and M[curStrId][comparePos] != unknownChar and M[curStrId][comparePos] != M[strId][strPos]:
					numDiff4[subLen - (l - k4 + 1)][strId][curStrId] += 1

				
				# Y = M[curStrId] and X = M[strId]
				# Compare M[strId][strPos] and M[curStrId][subLen + strPos - l]
				comparePos = subLen + strPos - l
				if comparePos >= 0 and M[curStrId][comparePos] != unknownChar and M[curStrId][comparePos] != M[strId][strPos]:
					numDiff4[subLen - (l - k4 + 1)][curStrId][strId] += 1 
					
################################################################################
#
def computeChangeInExpCount(M, strPos, strId, newVal, unknownChar, numDiff1, numDiff4, n, l, k1, k4, pascalPrefixSum):
	k = max(k1, k4)
	
	totalDiff = m_fraction.MyFraction(0, 1)
	# Compute the change contributed by the C1 constraint
	for curStrId in xrange(strId):
		if M[curStrId][strPos] == unknownChar:
			continue

		oldNumDiff = numDiff1[strId][curStrId]
		
		# Testing
		#if (strPos, oldNumDiff) != helper.countAssignedPos(M[curStrId], 0, l - 1, M[strId], 0, l - 1, unknownChar):
		#	raise RuntimeError("Your formula for (s, t) for C1 is WRONG")

		totalDiff += computeChangeInProb(l, k, strPos, oldNumDiff, newVal != M[curStrId][strPos], pascalPrefixSum)
	#if strPos == 1 and strId == 1:
	#	print "C1: totalDiff(1, 0) = " + str(totalDiff)

	# Compute the change contributed by the C4 constraint
	for subLen in xrange(l - k4 + 1, l):
		for curStrId in xrange(n):
			if curStrId == strId:
				continue

			# Consider M[strId][0 ... subLen - 1] and M[curStrId][l - subLen ... l - 1]
			if l - subLen + strPos < l and M[curStrId][l - subLen + strPos] != unknownChar:
				numAssigned = max(0, strPos - l + subLen)
				oldNumDiff = numDiff4[subLen - (l - k4 + 1)][strId][curStrId]

				# Testing
				#if (numAssigned, oldNumDiff) != helper.countAssignedPos(M[strId], 0, subLen - 1, M[curStrId], l - subLen, l - 1, unknownChar):
				#	raise RuntimeError("Your formula for (s, t) for C4 is WRONG")

				totalDiff += computeChangeInProb(subLen, k4 - (l - subLen), numAssigned, oldNumDiff, newVal != M[curStrId][l - subLen + strPos], pascalPrefixSum)
				#if strPos == 1 and strId == 1:
				#	print "C4(1): totalDiff(1, 0) = " + str(totalDiff) + " with subLen = " + str(subLen)
			# Consider M[strId][l - subLen ... l - 1] and M[curStrId][0 ... subLen - 1]
			if subLen + strPos - l >= 0 and M[curStrId][subLen + strPos - l] != unknownChar:
				numAssigned = max(0, strPos - l + subLen)
				oldNumDiff = numDiff4[subLen - (l - k4 + 1)][curStrId][strId]
				
				# Testing
				#if (numAssigned, oldNumDiff) != helper.countAssignedPos(M[curStrId], 0, subLen - 1, M[strId], l - subLen, l - 1, unknownChar):
				#	raise RuntimeError("Your formula for (s, t) for C4 is WRONG")

				totalDiff += computeChangeInProb(subLen, k4 - (l - subLen), numAssigned, oldNumDiff, newVal != M[curStrId][subLen + strPos - l], pascalPrefixSum)
				#if strPos == 1 and strId == 1:
				#	print "C4(2): totalDiff(1, 0) = " + str(totalDiff) + " with subLen = " + str(subLen)
	#if strPos == 1 and strId == 1:
	#	print "totalDiff(1, 0) = " + str(totalDiff)
	return totalDiff
			

################################################################################
# computeChangeInProb(l, k, s, t, diffIncrease, pascalPrefixSum):
#
#
# Input description:
def computeChangeInProb(l, k, s, t, diffIncrease, pascalPrefixSum):
	if t > k - 1:
		return m_fraction.MyFraction(0, 1)
	
	denom = helper.fastPower(2, l - s)
	numer = 1
	if t < k - 1:
		# Compute (l - s - 1) Choose (k - 1 - t)
		if k - 1 - t > l - s - 1:
			numer = 0
		else:
			numer = pascalPrefixSum[l - s - 1][k - 1 - t] - pascalPrefixSum[l - s - 1][k - 2 - t]
	if not diffIncrease:
		numer = -numer
	return m_fraction.MyFraction(numer, denom) 				
