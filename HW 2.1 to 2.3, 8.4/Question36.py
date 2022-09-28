import unittest
import math
import numpy
import numpy as np
import sympy as sp

def evalLegendreBasis1D( degree, variate ):
    degree += 1
    p = sp.Symbol('p')
    n = sp.Symbol('n')
    e = sp.Symbol('e')
    P = [p] * degree
    expr1 = e**n
    expr2 = sp.integrate(expr1 * (p), (e, -1, 1)) / sp.integrate((p) * (p), (e, -1, 1)) * p
    for i in range(degree):
        #P[i] = expr1.subs([(e,variate),(n,i)])
        P[i] = expr1.subs(n,i)
        for j in range(i-1):
            #P[i] -= expr2.subs([(e,variate),(p,P[j]),(n,i)])
            P[i] -= expr2.subs(p,P[i-1])
            print("Run " + str(i) + " " + str(j) + ": expr1 = " + str(expr1.subs(n,i)) + ": p = " + str(P[i-1].subs(n,i)))
        P[i] = P[i].subs([(n,i)])
    return P

print(evalLegendreBasis1D(degree = 4, variate = 1))
#evalLegendreBasis1D(degree = 1, variate = 1)

class Test_evalLegendreBasis1D( unittest.TestCase ):
    def test_basisAtBounds( self ):
        for p in range( 0, 2 ):
            if ( p % 2 == 0 ):
                self.assertAlmostEqual( first = evalLegendreBasis1D( degree = p, variate = -1 ), second = +1.0, delta = 1e-12 )
            else:
                self.assertAlmostEqual( first = evalLegendreBasis1D( degree = p, variate = -1 ), second = -1.0, delta = 1e-12 )
            self.assertAlmostEqual( first = evalLegendreBasis1D( degree = p, variate = +1 ), second = 1.0, delta = 1e-12 )

    def test_constant( self ):
        for x in numpy.linspace( -1, 1, 100 ):
            self.assertAlmostEqual( first = evalLegendreBasis1D( degree = 0, variate = x ), second = 1.0, delta = 1e-12 )

    def test_linear( self ):
        for x in numpy.linspace( -1, 1, 100 ):
            self.assertAlmostEqual( first = evalLegendreBasis1D( degree = 1, variate = x ), second = x, delta = 1e-12 )

    def test_quadratic_at_roots( self ):
        self.assertAlmostEqual( first = evalLegendreBasis1D( degree = 2, variate = -1.0 / math.sqrt(3.0) ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evalLegendreBasis1D( degree = 2, variate = +1.0 / math.sqrt(3.0) ), second = 0.0, delta = 1e-12 )

    def test_cubic_at_roots( self ):
        self.assertAlmostEqual( first = evalLegendreBasis1D( degree = 3, variate = -math.sqrt( 3 / 5 ) ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evalLegendreBasis1D( degree = 3, variate = 0 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evalLegendreBasis1D( degree = 3, variate = +math.sqrt( 3 / 5 ) ), second = 0.0, delta = 1e-12 )

