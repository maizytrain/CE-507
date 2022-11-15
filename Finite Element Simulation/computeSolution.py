import unittest
import numpy
from generate1DMesh import generateMesh1D


def computeSolution( target_fun, domain, num_elems, degree):
    node_coords, ien_array = generateMesh1D( xmin = domain[0], xmax = domain[1], num_elems = num_elems, degree = degree)
    test_solution = [None] * (num_elems * degree + 1)
    for i in range(num_elems * degree + 1):
        test_solution[i] = (target_fun)(node_coords[i])

    return test_solution, node_coords, ien_array







class Test_computeSolution( unittest.TestCase ):
    def test_single_linear_element_poly( self ):
        test_solution, node_coords, connect = computeSolution( target_fun = lambda x : x, domain = [-1.0, 1.0 ], num_elems = 1, degree = 1 )
        gold_solution = numpy.array( [ -1.0, 1.0 ] )
        self.assertTrue( numpy.allclose( test_solution, gold_solution ) )

    def test_single_quad_element_poly( self ):
        test_solution, node_coords, connect = computeSolution( target_fun = lambda x : x**2, domain = [-1.0, 1.0 ], num_elems = 1, degree = 2 )
        gold_solution = numpy.array( [ 1.0, 0.0, 1.0 ] )
        self.assertTrue( numpy.allclose( test_solution, gold_solution ) )

    def test_two_linear_element_poly( self ):
        test_solution, node_coords, connect = computeSolution( target_fun = lambda x : x**2, domain = [-1.0, 1.0 ], num_elems = 2, degree = 1 )
        gold_solution = numpy.array( [ 1.0, 0.0, 1.0 ] )
        self.assertTrue( numpy.allclose( test_solution, gold_solution ) )

    def test_four_quad_element_poly( self ):
        test_solution, node_coords, connect = computeSolution( target_fun = lambda x : x**2, domain = [-1.0, 1.0 ], num_elems = 4, degree = 1 )
        gold_solution = numpy.array( [ 1.0, 0.25, 0.0, 0.25, 1.0 ] )
        self.assertTrue( numpy.allclose( test_solution, gold_solution ) )