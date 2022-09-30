import unittest
import math
import numpy as np
import sympy as sp
from Lagrange import createLagrangeBasis1D


def computeNewtonCotesQuadrature( fun, num_points, degree ):
    calcMat = np.zeros([degree+1,degree+1])
    for i in range(degree+1):
        expr = sp.Poly(sp.expand(sp.simplify(createLagrangeBasis1D( degree = degree, basis_idx = i))))
        coeffs = expr.all_coeffs()
        for j in range(degree+1):
            calcMat[i,j] = coeffs[j]

    calcMat = np.transpose(calcMat)

    xBounds = np.linspace(-1,1,num_points+1)
    finals = 0
    zeta = sp.symbols('zeta')

    for i in range(num_points):
        eqExpr = calcSingleNewtonCotes(fun = fun, degree = degree, lowerBound = xBounds[i], upperBound = xBounds[i+1], calcMat = calcMat)
        finals += sp.integrate(eqExpr, (zeta, xBounds[i], xBounds[i+1]))

    return finals


def calcSingleNewtonCotes( fun, degree, lowerBound, upperBound, calcMat ):
    xValues = np.linspace(lowerBound,upperBound,degree+1)
    fValues = np.zeros(degree+1)
    for i in range(degree+1):
        fValues[i] = fun.subs(x,xValues[i])
    
    iMat = np.zeros([degree+1, degree+1])
    for i in range(degree + 1):
        for j in range(degree + 1):
            iMat[i, j] = xValues[i]**(degree - j)

    polyValues = np.linalg.solve(iMat,fValues)
    lValues = np.linalg.solve(calcMat,polyValues)

    zeta = sp.symbols('zeta')
    lExpr = zeta*0
    for i in range(degree+1):
        lExpr += lValues[i] * createLagrangeBasis1D( degree = degree, basis_idx = i)

    sExpr = (sp.simplify(lExpr))

    return sExpr





x = sp.symbols('x')
testExpr = 4 * x**4 - x + 2
print(computeNewtonCotesQuadrature(fun = testExpr, num_points = 5, degree = 2))

class Test_computeNewtonCotesQuadrature( unittest.TestCase ):
    def test_integrate_constant_one( self ):
        constant_one = 1 * x**0
        for degree in range( 1, 6 ):
            num_points = degree + 1
            self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = constant_one, num_points = num_points, degree = degree ), second = 2.0, delta = 1e-12 )

    def test_exact_poly_int( self ):
        for degree in range( 1, 6 ):
            num_points = degree + 1
            poly_fun = ( x + 1.0 ) ** degree
            indef_int = ( ( x + 1 ) ** ( degree + 1) ) / ( degree + 1 )
            def_int = indef_int.subs(x,1) - indef_int.subs(x,-1)
            self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = poly_fun, num_points = num_points, degree = degree ), second = def_int, delta = 1e-12 )

    def test_integrate_sin( self ):
        sin = sp.sin(x)
        for num_points in range( 1, 7 ):
            self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = sin, num_points = num_points, degree = 5 ), second = 0.0, delta = 1e-12 )

    def test_integrate_cos( self ):
        cos = sp.cos(x)
        self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = cos, num_points = 6, degree = 5 ), second = 2*math.sin(1), delta = 1e-4 )