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

    #print(calcMat)

    xValues = np.linspace(-1,1,degree+1)
    fValues = np.zeros(degree+1)
    for i in range(degree+1):
        fValues[i] = fun.subs(x,xValues[i])

    #print(fValues)

    iMat = np.zeros([degree+1, degree+1])
    for i in range(degree + 1):
        for j in range(degree + 1):
            iMat[i, j] = xValues[i]**(degree - j)

    polyValues = (np.linalg.solve(iMat,fValues))
    print(np.linalg.solve(calcMat,polyValues))

    fun = sp.Poly(fun)
    funDegree = sp.degree(fun) + 1
    sValues = np.zeros([funDegree, 1])

    #print(sValues)



    return

x = sp.symbols('x')
computeNewtonCotesQuadrature( fun = 2*x + 4, num_points = 1, degree = 2)

class Test_computeNewtonCotesQuadrature( unittest.TestCase ):
    def test_integrate_constant_one( self ):
        constant_one = lambda x : 1 * x**0
        for degree in range( 1, 6 ):
            num_points = degree + 1
            self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = constant_one, num_points = num_points ), second = 2.0, delta = 1e-12 )

    def test_exact_poly_int( self ):
        for degree in range( 1, 6 ):
            num_points = degree + 1
            poly_fun = lambda x : ( x + 1.0 ) ** degree
            indef_int = lambda x : ( ( x + 1 ) ** ( degree + 1) ) / ( degree + 1 )
            def_int = indef_int(1.0) - indef_int(-1.0)
            self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = poly_fun, num_points = num_points ), second = def_int, delta = 1e-12 )

    def test_integrate_sin( self ):
        sin = lambda x : math.sin(x)
        for num_points in range( 1, 7 ):
            self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = sin, num_points = num_points ), second = 0.0, delta = 1e-12 )

    def test_integrate_cos( self ):
        cos = lambda x : math.cos(x)
        self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = cos, num_points = 6 ), second = 2*math.sin(1), delta = 1e-4 )