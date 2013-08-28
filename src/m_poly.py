# Author: Nguyen Truong Duy
# Contact: truongduy134@gmail.com

##########################################################################
# This module implements a light-weight class of polynomial, called
# MyPoly

class MyPoly:
    # Constructor:
    def __init__(self, coeffarr = None):
        coeffarr = coeffarr or []
        numcoeff = len(coeffarr)
        if numcoeff == 0:
            # Create a zero polynomial
            coeffarr = [0]
            numcoeff = 1
        self.degree = numcoeff - 1
        self.coeff = list(coeffarr)
        self.remove_leading_zeroes()

    # Remove leading zero coefficients of the polynomial
    # except for the zero polynomial
    def remove_leading_zeroes(self):
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
        smallcoeff = []
        largecoeff = []
        minDegree = min(self.degree, other.degree)
        if self.degree < other.degree:
            smallcoeff = self.coeff
            largecoeff = other.coeff
        else:
            smallcoeff = other.coeff
            largecoeff = self.coeff

        newcoeff = list(largecoeff)
        for ind in xrange(minDegree + 1):
            newcoeff[ind] += smallcoeff[ind]

        return MyPoly(newcoeff)

    # Operator +=
    def __iadd__(self, other):
        newsum = self + other
        self.coeff = list(newsum.coeff)
        self.degree = newsum.degree
        return self

    # Constant multiplication
    def multiplyconst(self, constant):
        newlist = [x * constant for x in self.coeff]
        return MyPoly(newlist)

    # Operator -
    def __sub__(self, other):
        negatepoly = other.multiplyconst(-1)
        return self + negatepoly

    # Operator -=
    def __isub__(self, other):
        subpoly = self - other
        self.coeff = list(subpoly.coeff)
        self.degree = subpoly.degree
        return self

    # Operator *
    def __mul__(self, other):
        # Currently, we use an O(n^2) multiplication algorithm
        # for polynomial multiplication
        return self.multiply_quadratic(other)
    
    # Operator *=
    def __imul__(self, other):
        mulpoly = self * other
        self.coeff = list(mulpoly.coeff)
        self.degree = mulpoly.degree
        return self

    # Implmenent quadratic polynomial multiplication:
    #
    #
    def multiply_quadratic(self, other):
        if self.iszeropoly() or other.iszeropoly():
            # Any polynomial multiplied by 0 is 0
            return MyPoly()

        newdegree = self.degree * other.degree
        newlist = [0] * (newdegree + 1)
        for power in xrange(newdegree + 1):
            maxfirstpow = min(power, self.degree)
            for firstpow in xrange(maxfirstpow + 1):
                secondpow = power - firstpow
                if secondpow <= other.degree:
                    newlist[power] += self.coeff[firstpow] * other.coeff[secondpow]

        return MyPoly(newlist)

    # Return the result by multiplying this polynomial with a singleton of the form
    #        coeffsingleton * X^powsingleton where X is a variable
    def multiply_by_singleton(self, powsingleton, coeffsingleton):
        newcoeff = [0] * powsingleton
        newcoeff.extend(self.coeff)
        for ind in xrange(len(newcoeff)):
            newcoeff[ind] *= coeffsingleton
        return MyPoly(newcoeff)

    # Check if this polynomial is zero
    def iszeropoly(self):
        return self.degree == 0 and self.coeff[0] == 0 

    # Evaluate the polynomial at a given value
    def evaluate(self, val):
        result = 0
        powval = 1
        for deg in xrange(self.degree + 1):
            result += self.coeff[deg] * powval
            powval *= val
        return result
