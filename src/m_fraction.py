# Author: Nguyen Truong Duy
# Contact: truongduy134@gmail.com

##########################################################################
# This module implements a light-weight class of Fraction, called
# MyFraction
#
# This is because Python Fraction built-in class is too general and slow.
#
# MyFraction class currently supports the following functions
#	1) Constructor: Assume the inputs are integers
# 	   For example, MyFraction(2, 3) creates the fraction 2 / 3,
#	   and MyFraction(2) creates the fraction 2 / 1
#
#	2) Arithmetic operators: +, +=, *, *=, -
#
#	3) Comparison operators: ==, !=, >=, <=, <, >
#
#	4) Outputing
#
# Invariants about MyFraction objects:
#	1) The numerator and denominator are kept to the lowest terms.
#	   That is, their greatest common divisor is 1
#	2) The denominator is always positive
#

def gcd(a, b):
	while b:
		a, b = b, a % b
	return a

class MyFraction:
	# Constructor
	def __init__(self, numer = 0, denom = 1):
		if denom == 0:
			raise ZeroDivisionError("Denominator cannot be zero")

		self.numer = numer
		self.denom = denom
		self.simplify()

	# Simplify fraction
	def simplify(self):
		inumer = self.numer
		if inumer < 0:
			inumer = -inumer
		idenom = self.denom
		if idenom < 0:
			idenom = -idenom
		mgcd = gcd(inumer, idenom)

		if mgcd > 1:
			self.numer /= mgcd
			self.denom /= mgcd
		if self.denom < 0:
			self.numer = -self.numer
			self.denom = -self.denom

	# Add (+)
	def __add__(self, other):
		newnumer = self.numer * other.denom + self.denom * other.numer
		newdenom = self.denom * other.denom
		return MyFraction(newnumer, newdenom)

	# Self add (+=)
	def __iadd__(self, other):
		self.numer = self.numer * other.denom + self.denom * other.numer
		self.denom *= other.denom
		self.simplify()
		return self

	# Mul operator (*)
	def __mul__(self, other):
		newnumer = self.numer * other.numer
		newdenom = self.denom * other.denom
		return MyFraction(newnumer, newdenom)

	# Self mul operator (*=)
	def __imul__(self, other):
		self.numer *= other.numer
		self.denom *= other.denom
		self.simplify()
		return self

	# Subtract operator (-)
	def __sub__(self, other):
		newnumer = self.numer * other.denom - self.denom * other.numer
		newdenom = self.denom * other.denom
		return MyFraction(newnumer, newdenom)

	# Comparision operator
	def __eq__(self, other):
		return self.numer == other.numer and self.denom == other.denom
	
	def __neq__(self, other):
		return not self == other

	def __le__(self, other):
		return (self.numer * other.denom - self.denom * other.numer) <= 0
	
	def __gt__(self, other):
		return not self <= other

	def __ge__(self, other):
		return (self.numer * other.denom - self.denom * other.numer) >= 0

	def __lt__(self, other):
		return not self >= other

	# print function
	def __str__(self):
		return str(self.numer) + "/" + str(self.denom)
	def __repr__(self):
		return self.__str__()	
