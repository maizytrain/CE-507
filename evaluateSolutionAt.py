import unittest
import numpy
import sympy
from generate1DMesh import generateMesh1D

def evalLagrangeBasis1D(min, max, degree, basis_idx ):
    zeta = sympy.symbols('zeta')
    values = numpy.linspace(min, max, degree+1)
    expr = 1

    for i in range(degree+1):
        if (i != basis_idx):
            expr *= (zeta - values[i])/(values[basis_idx] - values[i])
    
    return expr

def solveLagrangeBasis1D(min, max, degree, func, variate):
    zeta = sympy.symbols('zeta')
    values = numpy.linspace(min, max, degree+1)
    expr = 0

    for i in range(degree+1):
        expr += evalLagrangeBasis1D(min, max, degree, i) * func(values[i])

    return expr.subs(zeta, variate)

#func = lambda x : x
#print (solveLagrangeBasis1D(-1, 1, 2, func, .78))


def evaluateSolutionAt(x, coeff, node_coords, connect, eval_basis, degree):
    zeta = sympy.symbols('zeta')
    expr = 0

    for i in range(len(node_coords)-degree):
        if (x != node_coords[len(node_coords)-1]):
            if (x >= node_coords[i] and x <= node_coords[i+degree]):
                for j in range(degree+1):
                    expr += evalLagrangeBasis1D(node_coords[i], node_coords[i+degree], degree, j) * coeff[(i+j)]
                return expr.subs(zeta, x)
        else:
            return coeff[len(coeff)-1]
    
    

#node_coords, connect = generateMesh1D( -1, 1, 2, 1 )
#coeff = numpy.array( [-1.0, 0, 1.0 ] )
#evaluateSolutionAt( x = -1, coeff = coeff, node_coords = node_coords, connect = connect, eval_basis = evalLagrangeBasis1D )
#print(evaluateSolutionAt( x = -1, coeff = coeff, node_coords = node_coords, connect = connect, eval_basis = evalLagrangeBasis1D ))

class Test_evaluateSolutionAt( unittest.TestCase ):
    def test_single_linear_element( self ):
        node_coords, connect = generateMesh1D( -1, 1, 1, 1 )
        coeff = numpy.array( [-1.0, 1.0 ] )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = -1.0, coeff = coeff, node_coords = node_coords, connect = connect, eval_basis = evalLagrangeBasis1D ), second = -1.0 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x =  0.0, coeff = coeff, node_coords = node_coords, connect = connect, eval_basis = evalLagrangeBasis1D ), second =  0.0 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = +1.0, coeff = coeff, node_coords = node_coords, connect = connect, eval_basis = evalLagrangeBasis1D ), second = +1.0 )

    def test_two_linear_elements( self ):
        node_coords, connect = generateMesh1D( -1, 1, 2, 1 )
        coeff = numpy.array( [ 1.0, 0.0, 1.0 ] )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = -1.0, coeff = coeff, node_coords = node_coords, connect = connect, eval_basis = evalLagrangeBasis1D ), second = +1.0 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x =  0.0, coeff = coeff, node_coords = node_coords, connect = connect, eval_basis = evalLagrangeBasis1D ), second =  0.0 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = +1.0, coeff = coeff, node_coords = node_coords, connect = connect, eval_basis = evalLagrangeBasis1D ), second = +1.0 )

    def test_single_quadratic_element( self ):
        node_coords, connect = generateMesh1D( -1, 1, 1, 2 )
        coeff = numpy.array( [+1.0, 0.0, 1.0 ] )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = -1.0, coeff = coeff, node_coords = node_coords, connect = connect, eval_basis = evalLagrangeBasis1D ), second = +1.0 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x =  0.0, coeff = coeff, node_coords = node_coords, connect = connect, eval_basis = evalLagrangeBasis1D ), second =  0.0 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = +1.0, coeff = coeff, node_coords = node_coords, connect = connect, eval_basis = evalLagrangeBasis1D ), second = +1.0 )

    def test_two_quadratic_elements( self ):
        node_coords, connect = generateMesh1D( -2, 2, 2, 2 )
        coeff = numpy.array( [ 1.0, 0.25, 0.5, 0.25, 1.0 ] )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = -2.0, coeff = coeff, node_coords = node_coords, connect = connect, eval_basis = evalLagrangeBasis1D ), second = +1.00 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = -1.0, coeff = coeff, node_coords = node_coords, connect = connect, eval_basis = evalLagrangeBasis1D ), second = +0.25 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x =  0.0, coeff = coeff, node_coords = node_coords, connect = connect, eval_basis = evalLagrangeBasis1D ), second = +0.50 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = +1.0, coeff = coeff, node_coords = node_coords, connect = connect, eval_basis = evalLagrangeBasis1D ), second = +0.25 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = +2.0, coeff = coeff, node_coords = node_coords, connect = connect, eval_basis = evalLagrangeBasis1D ), second = +1.00 )