import unittest
import numpy
import sympy
from generate1DMesh import generateMesh1D

def evalLagrangeBasis1D( variate, degree, basis_idx):
    zeta = sympy.symbols('zeta')
    values = numpy.linspace(-1,1, degree+1)
    expr = 1

    for i in range(degree+1):
        if (i != basis_idx):
            expr *= (zeta - values[i])/(values[basis_idx] - values[i])
    
    return expr.subs(zeta, variate)

def evaluateSolutionAt(x, coeff, node_coords, connect, eval_basis):
    zeta = sympy.symbols('zeta')
    
    for i in range(len(coeff)):
        if (x == node_coords[i]):
            return coeff[i]

    for i in range(len(coeff) - 1):
        if (x > node_coords[i] and x < node_coords[i+1]):
            return (node_coords[i] + eval_basis(x, 1, 1) * (node_coords[i+1] - node_coords[i]))
            

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