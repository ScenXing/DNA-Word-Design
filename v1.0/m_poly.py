# Author: Nguyen Truong Duy
# Contact: truongduy134@gmail.com

##########################################################################
# This module implements a light-weight class of polynomial, called
# MyPoly

class MyPoly:
	# Constructor:
	def __init__(self, coeffArr = []):
		numCoeff = len(coeffArr)
		if numCoeff == 0:
			# Create a zero polynomial
			coeffArr = [0]
			numCoeff = 1
		self.degree = numCoeff - 1
		self.coeff = list(coeffArr)
		self.removeLeadingZeroes()

	# Remove leading zero coefficients of the polynomial
	# except for the zero polynomial
	def removeLeadingZeroes(self):
		while self.degree > 0 and self.coeff[self.degree] == 0:
			self.degree -= 1
			self.coeff.pop()

	# print function
	def __str__(self):
		return self.coeff.__str__()

	def __repr__(self):
		return "Polynomial of degree = " + str(self.degree) + " and coeffcient = " + self.coeff.__str__()


	# Operator +
	def __add__(self, other):
		smallCoeff = []
		largeCoeff = []
		minDegree = min(self.degree, other.degree)
		if self.degree < other.degree:
			smallCoeff = self.coeff
			largeCoeff = other.coeff
		else:
			smallCoeff = other.coeff
			largeCoeff = self.coeff

		newCoeff = list(largeCoeff)
		for ind in xrange(minDegree + 1):
			newCoeff[ind] += smallCoeff[ind]

		return MyPoly(newCoeff)

	# Operator +=
	def __iadd__(self, other):
		newSum = self + other
		self.coeff = list(newSum.coeff)
		self.degree = newSum.degree
		return self

	# Constant multiplication
	def multiplyConst(self, constant):
		newList = [x * constant for x in self.coeff]
		return MyPoly(newList)

	# Operator -
	def __sub__(self, other):
		negatePoly = other.multiplyConst(-1)
		return self + negatePoly

	# Operator -=
	def __isub__(self, other):
		subPoly = self - other
		self.coeff = list(subPoly.coeff)
		self.degree = subPoly.degree
		return self

	# Operator *
	def __mul__(self, other):
		# Currently, we use an O(n^2) multiplication algorithm
		# for polynomial multiplication
		return self.multiplyQuadratic(other)
	
	# Operator *=
	def __imul__(self, other):
		mulPoly = self * other
		self.coeff = list(mulPoly.coeff)
		self.degree = mulPoly.degree
		return self

	# Implmenent quadratic polynomial multiplication:
	#
	#
	def multiplyQuadratic(self, other):
		if self.isZeroPoly() or other.isZeroPoly():
			# Any polynomial multiplied by 0 is 0
			return MyPoly()

		newDegree = self.degree * other.degree
		newList = [0] * (newDegree + 1)
		for power in xrange(newDegree + 1):
			maxFirstPow = min(power, self.degree)
			for firstPow in xrange(maxFirstPow + 1):
				secondPow = power - firstPow
				if secondPow <= other.degree:
					newList[power] += self.coeff[firstPow] * other.coeff[secondPow]

		return MyPoly(newList)

	# Return the result by multiplying this polynomial with a singleton of the form
	#		coeffSingleton * X^powSingleton where X is a variable
	def multiplyBySingleton(self, powSingleton, coeffSingleton):
		newCoeff = [0] * powSingleton
		newCoeff.extend(self.coeff)
		for ind in xrange(len(newCoeff)):
			newCoeff[ind] *= coeffSingleton
		return MyPoly(newCoeff)

	# Check if this polynomial is zero
	def isZeroPoly(self):
		return self.degree == 0 and self.coeff[0] == 0 

	# Evaluate the polynomial at a given value
	def evaluate(self, val):
		result = 0
		powVal = 1
		for deg in xrange(self.degree + 1):
			result += self.coeff[deg] * powVal
			powVal *= val
		return result
