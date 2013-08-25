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
# initlength14(n, k1, k4, delta):
#
# The implementation is based on Lemma 7.
#
def initlength14(n, k1, k4, delta = 0.1):
	"""
	Return an integer l such that there exists a partially matrix M of size n x l such that ExpCount(M, k1, k4) > n * (n - 1) / 2 * (1 + (2 * (k4 - 1)) - 1; hence we can obtain a (k1, k4) distance matrix of size n x l from M

	Inputs:
	+ n, k1, k4: integers
	+ delta: a real number. By default, delta = 0.1

	Output:
	+ l: an integer satisfying the above condition.

	Assumption:
	n >= 2 and k = max(k1, k4) >= 1, delta > 0. Otherwise, the correctness of the function is not guaranteed

	Methodology
	+ Let c1 = 2 + delta
	+ Let c2 = 0.5 * c1 * {log(c1 / ((c1 - 2) * ln(2))) + 2.5 - 1 / ln(2)} 
	+ Set l = ceil(c1 * log(n) + c2 * k) 
	"""

	c1 = 2 + delta
	c2 = c1 / 2.0 * (math.log(c1 / ((c1 - 2) * math.log(2)), 2) + 2.5 - 1.0 / math.log(2))
	k = max(k1, k4)

	return int(math.ceil(c1 * math.log(n, 2) + c2 * k))

##############################################################################################
# binarysearchlength14(n, k1, k4, initlength, pascalprefixsum):
#	- Returns an integer which is the minimum length l such that the
#	  distance matrix M of size n x l satisfies Lemma 3, i.e.
#		ExpCount(M, k1, k4) > nC2 * (1 + 2 * (k4 - 1)) - 1
#
# The implementation is based on Lemma 10 and Lemma 11 in the paper. It uses
#	binary search to find the minimum length based on the initial length given
#	as an input.
#
# Time complexity: O(k4 * log(initlength)) 
#	(This does not take into account the complexity of arithmetics)
def binarysearchlength14(n, k1, k4, initlength, pascalprefixsum):
	"""
	Find an integer l which is the minimum length such that there exists a partially assigned matrix M of size n x l satisfying ExpCount(M, k1, k4) > n * (n - 1) / 2 * (1 + (2 * (k4 - 1)) - 1

	Inputs:
	+ n, k1, k4: integers
	+ initlength: a positive integer such that there exists a partially assigned matrix M of size n x initlength satisfying the above condition (e.g. initlength can be the result returned by the function initlength14)
	+ pascalprefixsum: a 2D array (which is essentially a Python list of lists of integers) with at least initlength + 1 rows. Row i has (i + 1) entries. And pascalprefixsum[i][j] = sum ( i Choose h ) for h = 0, ..., j

	Output:
	+ l: an integer indicating the minimum length satisfying the above condition 
	"""
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
		"""
		"""
		result = m_fraction.MyFraction(pascalprefixsum[l][k - 1], helper.fastpower(2, l))
		
		powtwo = helper.fastpower(2, l - k4 + 1)
		sumtemp = m_fraction.MyFraction(0, 1)
		for i in xrange(l - k4 + 1, l):
			sumtemp += m_fraction.MyFraction(pascalprefixsum[i][k4 - (l - i) - 1], powtwo)
			powtwo *= 2
 		sumtemp *= m_fraction.MyFraction(2)

		result += sumtemp
		result *= m_fraction.MyFraction((n * (n - 1)) / 2)
		return result

	leftlen = 1
	rightlen = initlength
	minfound = initlength

	while leftlen <= rightlen:
		midlen = leftlen + ((rightlen - leftlen) >> 1)
		val = compute(midlen)
		if val.numer < val.denom:
			minfound = midlen
			rightlen = midlen - 1
		else:
			leftlen = midlen + 1

	return minfound

#####################################################################################
# compute_expcount_empty(n, l, k1, k4, pascalprefixsum):
#	- Computes and returns the value of ExpCount(M, k1, k4) when the matrix
#		M of size n x l is empty. 
#
# The formula is mentioned in the proof of Lemma 10 in the paper.
#
# Input description:
#	- pascalprefixsum: is a 2D array with at least l + 1 rows. Row i has (i + 1) entries
#	- pascalprefixsum[i][j] = sum ( iCh ) for h = 0, ..., j
# Note: iCh denotes i Choose h   
def compute_expcount_empty(n, l, k1, k4, pascalprefixsum):
	"""
	Compute the value of ExpCount(M, k1, k4) where M is an empty partially assigned matrix of size n x l (i.e. every entry is an unknown).

	Inputs:
	+ n, l: the number of rows and columns of M respectively
	+ k1, k4: parameters of C1 and C4 constraint respectively
	+ pascalprefixsum: a 2D array (which is essentially a Python list of lists of integers) with at least initlength + 1 rows. Row i has (i + 1) entries. And pascalprefixsum[i][j] = sum ( i Choose h ) for h = 0, ..., j

	Output:
	+  
	"""

	k = max(k1, k4)	
	firstterm = m_fraction.MyFraction(((n * (n - 1)) / 2) * (1 + 2 * (k4 - 1)), 1)

	secondterm = m_fraction.MyFraction(pascalprefixsum[l][k - 1], helper.fastpower(2, l))
		
	powtwo = helper.fastpower(2, l - k4 + 1)
	sumtemp = m_fraction.MyFraction(0, 1)
	for i in xrange(l - k4 + 1, l):
		sumtemp += m_fraction.MyFraction(pascalprefixsum[i][k4 - (l - i) - 1], powtwo)
		powtwo *= 2
 	sumtemp *= m_fraction.MyFraction(2)
	secondterm += sumtemp
	secondterm *= m_fraction.MyFraction((n * (n - 1)) / 2)

	return firstterm - secondterm

#####################################################################################
# breakrun(str1, maxlenrun):
#
# Requirement: 
#	- maxlenrun >= 2
#	- str1 contains characters in {'0', '1', 'A', 'C', 'G', 'T'} only. Otherwise,
#	  the correctness of the algorithm is not guaranteed.
#
# The implementation is originally based on Algorithm 2 of the paper. However, it behaves
#	slightly diffrent from the algorithm described in the paper. Also, it works even
#	if the input string has odd length.
#
# Detail:
# 	1) Let L be the length of str1. Split the input string into two halves:
#			X = str1[0 ... floor(L / 2) - 1] and
#			Y = str2[floor(L / 2) ... L - 1]
#	2) Do the following procedure to X and the reverse of Y: After every maxlenrun - 1
#	   characters, or after we reach the last character of the string, insert a new character 
#          which is complementary to the last character of the run we are considering.
#	   Let the resulted string be newX and newRY respectively.
#	3) Concatenate newX and the reverse of newRY to form str2
#
# Time complexity: O(L)
def breakrun(str1, maxlenrun):
	"""
	Break run in str1 to obtain a new string str2 such that str2 has no more than maxlenrun consecutive bases

	Inputs:
	 + str1: a Python string. It should contain only characters in {'0', '1', 'A', 'C', 'G', 'T'}
	 + maxlenrun: an integer indicating the maximum length of a run of the same characters

	Output:
	 + str2: a Python string derived from str1 and that has no more than maxlenrun consecutive same characters.

	Exception:
	 + Throws an exception if maxlenrun <= 1

	Note: The function does not work correctly if the input string contains any letter different from '0', '1', 'A', 'C', 'G', 'T'
	"""

	lenstr = len(str1)
	if maxlenrun <= 1:
		raise RuntimeError("maxlenrun must be at least 2")

	# Break run in the first half
	leftcharlist = []
	lenrun = 0
	for ind in xrange(lenstr >> 1):
		leftcharlist.append(str1[ind])
		lenrun += 1
		if (lenrun == maxlenrun - 1) or (ind == (lenstr >> 1) - 1):
			leftcharlist.append(helper.get_complement_letter(str1[ind]))
			lenrun = 0

	# Break run in the second half
	rightcharlist = []
	for ind in xrange(lenstr - 1, (lenstr >> 1) - 1, -1):
		rightcharlist.append(str1[ind])
		lenrun += 1
		if (lenrun == maxlenrun - 1) or (ind == (lenstr >> 1)):
			rightcharlist.append(helper.get_complement_letter(str1[ind]))
			lenrun = 0
	rightcharlist.reverse()

	# Return the resulted string
	leftcharlist.extend(rightcharlist)
	return ''.join(leftcharlist)

#########################################################################################
# update_numdiff(M, strpos, strid, unknownchar, numdiff1, numdiff4, n, l, k4):
#	- Update the difference matrices when the newest fill-in entry is M[strid][strpos]
#		assuming the entries in M are filled in the order top to bottom, 
#		then left to right.
#
# Input description:
#	M: a partially assigned k1, k4-distance matrix of size n x l
#	unknownchar: the character denoting unassigned entries in M
#	strpos, strid: M[strid][strpos] is the newest assigned entry in M
#	k4: the parameter of C4 constraint
#
#	numdiff1: a 2D array that has n rows. Row i has exactly (i + 1) columns (i = 0, ..., n - 1)
#		  numdiff1[a][b] = the number of positions k in strings M[a] and M[b] (b <= a) 
#			such that
#			M[a][k] and M[b][k] are not unknown, and
#			M[a][k] != M[b][k]
#
#	numdiff4: a 3D array of dimension (k4 - 1) x n x n
#		  numdiff4[x][a][b] = the number of positions k in substrings 
#			M[a][0 ... l - k4 + x] and M[b][l - (l - k4 + 1 + x) ... l - 1] such that
#				M[a][k] and M[b][l - (l - k4 + 1 + x) + k] are not unknown, and
#				M[a][k] != M[b][l - (l - k4 + 1 + x) + k]
#
# Time complexity: O(n * k)
def update_numdiff(M, strpos, strid, unknownchar, numdiff1, numdiff4, n, l, k4):
	"""
	Update the difference matrices when the newest fill-in entry is M[strid][strpos] assuming the entries in M are filled in the order top to bottom, then left to right.

	Inputs:
	 + n: the number of rows in M
	 + l: the number of columns in M
	 + M: a partially assigned k1, k4-distance matrix of size n x l
	 + 
	"""

	# Update numdiff1
	for curstrid in xrange(0, strid):
		if M[curstrid][strpos] != M[strid][strpos]:
			numdiff1[strid][curstrid] += 1

	# Update numdiff4
	for sublen in xrange(l - k4 + 1, l):
		for curstrid in xrange(0, n):
			if curstrid != strid:
				# Y = M[strid] and X = M[curstrid]
				# Compare M[strid][strpos] and M[curstrid][l - sublen + strpos]
				comparepos = l - sublen + strpos
				if comparepos < l and M[curstrid][comparepos] != unknownchar and M[curstrid][comparepos] != M[strid][strpos]:
					numdiff4[sublen - (l - k4 + 1)][strid][curstrid] += 1

				
				# Y = M[curstrid] and X = M[strid]
				# Compare M[strid][strpos] and M[curstrid][sublen + strpos - l]
				comparepos = sublen + strpos - l
				if comparepos >= 0 and M[curstrid][comparepos] != unknownchar and M[curstrid][comparepos] != M[strid][strpos]:
					numdiff4[sublen - (l - k4 + 1)][curstrid][strid] += 1 
					
################################################################################
#
def compute_change_in_expcount(M, strpos, strid, newval, unknownchar, numdiff1, numdiff4, n, l, k1, k4, pascalprefixsum):
	k = max(k1, k4)
	
	totaldiff = m_fraction.MyFraction(0, 1)
	# Compute the change contributed by the C1 constraint
	for curstrid in xrange(strid):
		if M[curstrid][strpos] == unknownchar:
			continue

		oldnumdiff = numdiff1[strid][curstrid]

		totaldiff += compute_change_in_prob(l, k, strpos, oldnumdiff, newval != M[curstrid][strpos], pascalprefixsum)

	# Compute the change contributed by the C4 constraint
	for sublen in xrange(l - k4 + 1, l):
		for curstrid in xrange(n):
			if curstrid == strid:
				continue

			# Consider M[strid][0 ... sublen - 1] and M[curstrid][l - sublen ... l - 1]
			if l - sublen + strpos < l and M[curstrid][l - sublen + strpos] != unknownchar:
				numassigned = max(0, strpos - l + sublen)
				oldnumdiff = numdiff4[sublen - (l - k4 + 1)][strid][curstrid]

				totaldiff += compute_change_in_prob(sublen, k4 - (l - sublen), numassigned, oldnumdiff, newval != M[curstrid][l - sublen + strpos], pascalprefixsum)
				
			# Consider M[strid][l - sublen ... l - 1] and M[curstrid][0 ... sublen - 1]
			if sublen + strpos - l >= 0 and M[curstrid][sublen + strpos - l] != unknownchar:
				numassigned = max(0, strpos - l + sublen)
				oldnumdiff = numdiff4[sublen - (l - k4 + 1)][curstrid][strid]

				totaldiff += compute_change_in_prob(sublen, k4 - (l - sublen), numassigned, oldnumdiff, newval != M[curstrid][sublen + strpos - l], pascalprefixsum)
				
	return totaldiff
			

################################################################################
# compute_change_in_prob(l, k, s, t, diffincrease, pascalprefixsum):
#
#
# Input description:
def compute_change_in_prob(l, k, s, t, diffincrease, pascalprefixsum):
	if t > k - 1:
		return m_fraction.MyFraction(0, 1)
	
	denom = helper.fastpower(2, l - s)
	numer = 1
	if t < k - 1:
		# Compute (l - s - 1) Choose (k - 1 - t)
		if k - 1 - t > l - s - 1:
			numer = 0
		else:
			numer = pascalprefixsum[l - s - 1][k - 1 - t] - pascalprefixsum[l - s - 1][k - 2 - t]
	if not diffincrease:
		numer = -numer
	return m_fraction.MyFraction(numer, denom)
