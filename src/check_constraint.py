# Author: Nguyen Truong Duy
# Email: truongduy134@gmail.com

# Self-written python modules
import helper	
import algo_subroutine
import free_energy_routine

# Standard modules
import math

####################################################################################
# checkc1(strset, k1):
#	- Returns true if the input set of strings satisfies Basic Hamming Constraint
#	- Returns false otherwise
def checkc1(strset, k1):
	numstr = len(strset)
	
	for indone in xrange(numstr):
		for indtwo in xrange(indone + 1, numstr):
			if not checkc1each(strset[indone], strset[indtwo], k1):
				return False
	return True

# checkc1each(str1, str2, k1):
#	- Returns true if the two input strings satisfy Basic Hamming Constraint 
#	- Returns false otherwise
def checkc1each(str1, str2, k1):
	if k1 <= 0:
		return True
	return helper.hamming_dist(str1, str2) >= k1

##########################################################################################
# checkc2(strset, k2):
#	- Returns true if the input set of strings satisfies Reverse Complementary Constraint
def checkc2(strset, k2):
	numstr = len(strset)

	for indone in xrange(numstr):
		for indtwo in xrange(numstr):
			if indone != indtwo:
				if not checkc2each(strset[indone], strset[indtwo], k2):
					return False
	return True

# checkc2each(str1, str2, k2):
#	- Returns true if the two input strings satisfy Reverse Complementary Constraint 
#	- Returns false otherwise
def checkc2each(str1, str2, k2):
	if k2 <= 0:
		return True
	revcompstr2 = helper.reverse_str(helper.complement_str(str2))
	return helper.hamming_dist(str1, revcompstr2) >= k2

###############################################################################################
# checkc3(strset, k3):
#	- Returns true if the input set of strings satisfies Self Reverse Complementary Constraint
#	- Returns false otherwise 
def checkc3(strset, k3):
	numstr = len(strset)

	for ind in xrange(numstr):
		if not checkc3each(strset[ind], k3):
			return False
	return True

# checkc3each(str1, k3):
#	- Returns true if the input string satisfies Self Reverse Complementary Constraint 
#	- Returns false otherwise
def checkc3each(str1, k3):
	if k3 <= 0:
		return True
	return checkc2each(str1, str1, k3)

################################################################################################
# checkc4(strset, k4):
#	- Returns true if the input set of strings satisfies Shifting Hamming Constraint
#	- Returns false otherwise 
def checkc4(strset, k4):
	numstr = len(strset)

	for indone in xrange(numstr):
		for indtwo in xrange(numstr):
			if indone != indtwo:
				if not checkc4each(strset[indone], strset[indtwo], k4):
					print "False with indone = " + str(indone) + " and indtwo = " + str(indtwo)
					return False
	return True

# checkc4each(str1, str2, k4):
#	- Returns true if the two input strings satisfy Shifting Hamming Constraint 
#	- Returns false otherwise
def checkc4each(str1, str2, k4):
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
			print "False at i = " + str(sublen)
			return False
	
	return True

############################################################################################################
# checkc5(strset, k5):
#	- Returns true if the input set of strings satisfies Shifting Reverse Complementary Constraint
#	- Returns false otherwise 
def checkc5(strset, k5):
	numstr = len(strset)
	
	for indone in xrange(numstr):
		for indtwo in xrange(numstr):
			if indone != indtwo:
				if not checkc5each(strset[indone], strset[indtwo], k5):
					return False
	return True

# checkc5each(str1, str2, k5):
#	- Returns true if the two input strings satisfy Shifting Reverse Complementary Constraint
#	- Returns false otherwise
def checkc5each(str1, str2, k5):
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
# checkc6(strset, k6):
#	- Returns true if the input set of strings satisfies Shifting Self Reverse Complementary Constraint
#	- Returns false otherwise 
def checkc6(strset, k6):
	for mStr in strset:
		if not checkc6each(mStr, k6):
			return False
	return True

# checkc6each(str1, k6):
#	- Returns true if the input string satisfies Shifting Self Reverse Complementary Constraint
#	- Returns false otherwise
def checkc6each(str1, k6):
	return checkc5each(str1, str1, k6)

##################################################################################################
# checkc7(strset, ratio):
#	- Returns true if the input set of strings satisfies GC Content Constraint
#	- Returns false otherwise 
def checkc7(strset, ratio):
	for mStr in strset:
		if not checkc7each(mStr, ratio):
			return False
	return True

# checkc7each(str1, ratio):
#	- Returns true if the input string satisfies GC Content Constraint
#	- Returns false otherwise
def checkc7each(str1, ratio):
	if ratio < 0 or ratio > 1:
		return False

	len1 = len(str1)
	nummustpresentCeil = int(math.ceil(ratio * len1))
	nummustpresentFloor = int(ratio * len1)

	numpresent = 0
	for mchar in str1:
		if mchar == 'G' or mchar == 'C':
			numpresent += 1

	return numpresent == nummustpresentFloor or numpresent == nummustpresentCeil

##################################################################################################
# checkc8(strset, maxlenrun):
#	- Returns true if the input set of strings satisfyies Consecutive Base Constraint Constraint
#	- Returns false otherwise
def checkc8(strset, maxlenrun):
	if maxlenrun <= 0:
		return False
	for mStr in strset:
		if not checkc8each(mStr, maxlenrun):
			return False
	return True

# checkc8each(str1, maxlenrun):
#	- Returns true if the input string satisfies Consecutive Base Constraint Constraint
#	- Returns false otherwise
def checkc8each(str1, maxlenrun):
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
# checkc9(strset, sigma, pairwise_energy):
#	- Returns true if the input set of strings satisfies Free Energy Constraint
#	- Returns false otherwise
def checkc9(strset, sigma, pairwise_energy):
	numstr = len(strset)

	for indone in xrange(numstr):
		for indtwo in xrange(indone + 1, numstr):
			if not checkc9each(strset[indone], strset[indtwo], sigma, pairwise_energy):
				return False
	return True

# checkc9each(str1, str2, sigma, pairwise_energy):
#	- Returns true if the two input strings str1 and str2 satisfy the Free Energy Constraint
#	- Returns false otherwise
def checkc9each(str1, str2, sigma, pairwise_energy):
	free_energy1 = free_energy_routine.computeFreeEnergy(str1, pairwise_energy)
	free_energy2 = free_energy_routine.computeFreeEnergy(str2, pairwise_energy)

	return not (math.fabs(free_energy1 - free_energy2) > sigma)

###########################################################################################
#
def checkconstraintset(strset, mapTypeToParam):
	checkfunclist = [checkc1, checkc2, checkc3, checkc4, checkc5, checkc6, checkc7, checkc8, checkc9]
	
	result = True
	for typekey in mapTypeToParam:
		if typekey >= 1 and typekey <= 9:
			if typekey != 9:
				result = checkfunclist[typekey - 1](strset, mapTypeToParam[typekey])
			else:
				result = checkfunclist[typekey - 1](strset, mapTypeToParam[typekey][0], mapTypeToParam[typekey][1])

			if not result:
				print typekey
				return False

	return True

###########################################################################################
#
def checkGenDnaWord1To6And9Algo(strset, mapTypeToParam):
	pairwise_energy = free_energy_routine.createPairwiseEnergyFunc(mapTypeToParam[9])
	maxPairwise = free_energy_routine.findExtremePairwiseEnergy(pairwise_energy, max)
	minPairwise = free_energy_routine.findExtremePairwiseEnergy(pairwise_energy, min)
	D = maxPairwise - minPairwise
	checkMapTypeToParam = mapTypeToParam.copy()
	checkMapTypeToParam[9] = (4 * D + maxPairwise, pairwise_energy)
	return checkConstraintSet(strset, checkMapTypeToParam)

###########################################################################################
# checklength14(l, n, k1, k4):
#	- Returns true if l satisfies Lemma 6 with inputs n, k1, k4 (This is to guarantee
#		Lemma 3 holds)
#	- Returns false otherwise
#
# This is implemented based on Lemma 6
def checklength14(l, n, k1, k4):
	k = max(k1, k4)
	if l < 2 * k:
		return False

	expr = l - k * math.log(math.e, 2) - k * math.log(l * 1.0 / k, 2) - 2 * math.log(n, 2) - 2 * math.log(k, 2)
	if expr > 0:
		return True
	return False

#####################################################################################
# checkExpCountEmpty(n, l, k1, k4):
#	- Returns true if the input instance (n, l, k1, k4) satisfies the condition
#		in Lemma 3, i.e.
#			ExpCount(M, k1, k4) > nC2 * (1 + 2(k4 - 1)) - 1
#	- Returns false otherwise
def checkExpCountEmpty(n, l, k1, k4):
	pascalPrefixSum = helper.generatePascalTriangle(l)
	helper.getPrefixSumPascalTriangle(pascalPrefixSum)

	expCount = algo_subroutine.computeExpCountEmpty(n, l, k1, k4, pascalPrefixSum)
	rhs = (n * (n - 1) / 2) * (1 + 2 * (k4 - 1)) - 1
	return expCount > rhs
