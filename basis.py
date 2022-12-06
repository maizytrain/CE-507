import unittest
import numpy
import sympy
import math

def evaluateBernsteinBasis1D(variate, degree, basis_idx):
    variate = (1 + variate) * .5
    zeta = sympy.symbols('zeta')
    expr = createBernsteinBasis1D(degree, basis_idx)
    val = expr.subs(zeta, variate)
    return val

def createBernsteinBasis1D(degree, basis_idx):
    zeta = sympy.symbols('zeta')
    scalar = math.factorial(degree) / (math.factorial(basis_idx) * math.factorial(degree - basis_idx))
    expr = scalar * zeta ** basis_idx * (1 - zeta) ** (degree - basis_idx)
    return expr

# print(evaluateBernsteinBasis1D(0,2,0))
# print(createBernsteinBasis1D(2,0))

def affine_mapping_1D(new_domain, old_domain, loc):
    mod = (old_domain[1]-old_domain[0])/(new_domain[1]-new_domain[0])
    val = (loc - new_domain[0]) * mod + old_domain[0]
    return val

def evalBernsteinBasisDeriv(degree, basis_idx, deriv, domain, variate):
    variate = affine_mapping_1D(domain, [0, 1], variate)
    zeta = sympy.symbols('zeta')
    expr = createBernsteinBasis1D(degree, basis_idx)
    for i in range(deriv):
        expr = sympy.diff(expr)
    val = expr.subs(zeta, variate)
    return val

def createBernsteinBasis1DDeriv(degree, basis_idx, deriv):
    zeta = sympy.symbols('zeta')
    scalar = math.factorial(degree) / (math.factorial(basis_idx) * math.factorial(degree - basis_idx))
    expr = scalar * zeta ** basis_idx * (1 - zeta) ** (degree - basis_idx)
    # for i in range(deriv):
    expr = sympy.diff(expr, zeta, deriv)
    return expr


def evalSplineBasisDeriv1D( extraction_operator, basis_idx, deriv, domain, variate ):
    degree = len(extraction_operator) - 1
    d = numpy.zeros(degree + 1)
    for i in range(degree + 1):
        d[i] = evalBernsteinBasisDeriv(degree, i, deriv, domain, variate) / math.pow((domain[1]-domain[0]), deriv)
    val = numpy.matmul(extraction_operator, d)
    return val[basis_idx]

#C = numpy.eye( 2 )
#print(evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 1, domain = [ -1, 1 ], variate = -1.0))



# x = [ -1.0 / math.sqrt(3.0) , 1.0 / math.sqrt(3.0) ]
# x = [ affine_mapping_1D( [-1, 1], [0, 1], xi ) for xi in x ]
# print(evalBernsteinBasisDeriv( degree = 2, basis_idx = 0, deriv = 1, domain = [0, 1], variate = x[0] ))




class Test_evaluateBernsteinBasis1D( unittest.TestCase ):
    def test_linearBernstein( self ):
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = -1, degree = 1, basis_idx = 0 ), second = 1.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = -1, degree = 1, basis_idx = 1 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = +1, degree = 1, basis_idx = 0 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = +1, degree = 1, basis_idx = 1 ), second = 1.0, delta = 1e-12 )

    def test_quadraticBernstein( self ):
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = -1, degree = 2, basis_idx = 0 ), second = 1.00, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = -1, degree = 2, basis_idx = 1 ), second = 0.00, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = -1, degree = 2, basis_idx = 2 ), second = 0.00, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate =  0, degree = 2, basis_idx = 0 ), second = 0.25, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate =  0, degree = 2, basis_idx = 1 ), second = 0.50, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate =  0, degree = 2, basis_idx = 2 ), second = 0.25, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = +1, degree = 2, basis_idx = 0 ), second = 0.00, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = +1, degree = 2, basis_idx = 1 ), second = 0.00, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = +1, degree = 2, basis_idx = 2 ), second = 1.00, delta = 1e-12 )
