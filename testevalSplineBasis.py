import unittest
import numpy
import basis


class Test_evalSplineBasisDeriv1D( unittest.TestCase ):
       def test_C0_linear_0th_deriv_at_nodes( self ):
              C = numpy.eye( 2 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 0, domain = [ 0, 1 ], variate = 0.0 ), second = 1.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 0, domain = [ 0, 1 ], variate = 0.0 ), second = 0.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 0, domain = [ 0, 1 ], variate = 1.0 ), second = 0.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 0, domain = [ 0, 1 ], variate = 1.0 ), second = 1.0 )

       def test_C0_linear_1st_deriv_at_nodes( self ):
              C = numpy.eye( 2 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 1, domain = [ 0, 1 ], variate = 0.0 ), second = -1.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 1, domain = [ 0, 1 ], variate = 0.0 ), second = +1.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 1, domain = [ 0, 1 ], variate = 1.0 ), second = -1.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 1, domain = [ 0, 1 ], variate = 1.0 ), second = +1.0 )

       def test_C1_quadratic_0th_deriv_at_nodes( self ):
              C = numpy.array( [ [ 1.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.5 ], [ 0.0, 0.0, 0.5 ] ] )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 0, domain = [ 0, 1 ], variate = 0.0 ), second = 1.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 0, domain = [ 0, 1 ], variate = 0.0 ), second = 0.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 0, domain = [ 0, 1 ], variate = 0.0 ), second = 0.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 0, domain = [ 0, 1 ], variate = 0.5 ), second = 0.25 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 0, domain = [ 0, 1 ], variate = 0.5 ), second = 0.625 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 0, domain = [ 0, 1 ], variate = 0.5 ), second = 0.125 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 0, domain = [ 0, 1 ], variate = 1.0 ), second = 0.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 0, domain = [ 0, 1 ], variate = 1.0 ), second = 0.5 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 0, domain = [ 0, 1 ], variate = 1.0 ), second = 0.5 )

       def test_C1_quadratic_1st_deriv_at_nodes( self ):
              C = numpy.array( [ [ 1.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.5 ], [ 0.0, 0.0, 0.5 ] ] )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 1, domain = [ 0, 1 ], variate = 0.0 ), second = -2.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 1, domain = [ 0, 1 ], variate = 0.0 ), second = +2.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 1, domain = [ 0, 1 ], variate = 0.0 ), second = +0.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 1, domain = [ 0, 1 ], variate = 0.5 ), second = -1.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 1, domain = [ 0, 1 ], variate = 0.5 ), second = +0.5 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 1, domain = [ 0, 1 ], variate = 0.5 ), second = +0.5 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 1, domain = [ 0, 1 ], variate = 1.0 ), second = +0.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 1, domain = [ 0, 1 ], variate = 1.0 ), second = -1.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 1, domain = [ 0, 1 ], variate = 1.0 ), second = +1.0 )

       def test_C1_quadratic_2nd_deriv_at_nodes( self ):
              C = numpy.array( [ [ 1.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.5 ], [ 0.0, 0.0, 0.5 ] ] )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 2, domain = [ 0, 1 ], variate = 0.0 ), second = +2.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 2, domain = [ 0, 1 ], variate = 0.0 ), second = -3.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 2, domain = [ 0, 1 ], variate = 0.0 ), second = +1.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 2, domain = [ 0, 1 ], variate = 0.5 ), second = +2.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 2, domain = [ 0, 1 ], variate = 0.5 ), second = -3.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 2, domain = [ 0, 1 ], variate = 0.5 ), second = +1.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 2, domain = [ 0, 1 ], variate = 1.0 ), second = +2.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 2, domain = [ 0, 1 ], variate = 1.0 ), second = -3.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 2, domain = [ 0, 1 ], variate = 1.0 ), second = +1.0 )

       def test_biunit_C0_linear_0th_deriv_at_nodes( self ):
              C = numpy.eye( 2 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 0, domain = [ -1, 1 ], variate = -1.0 ), second = 1.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 0, domain = [ -1, 1 ], variate = -1.0 ), second = 0.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 0, domain = [ -1, 1 ], variate = +1.0 ), second = 0.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 0, domain = [ -1, 1 ], variate = +1.0 ), second = 1.0 )

       def test_biunit_C0_linear_1st_deriv_at_nodes( self ):
              C = numpy.eye( 2 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 1, domain = [ -1, 1 ], variate = -1.0 ), second = -0.5 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 1, domain = [ -1, 1 ], variate = -1.0 ), second = +0.5 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 1, domain = [ -1, 1 ], variate = +1.0 ), second = -0.5 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 1, domain = [ -1, 1 ], variate = +1.0 ), second = +0.5 )

       def test_biunit_C1_quadratic_0th_deriv_at_nodes( self ):
              C = numpy.array( [ [ 1.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.5 ], [ 0.0, 0.0, 0.5 ] ] )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 0, domain = [ -1, 1 ], variate = -1.0 ), second = 1.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 0, domain = [ -1, 1 ], variate = -1.0 ), second = 0.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 0, domain = [ -1, 1 ], variate = -1.0 ), second = 0.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 0, domain = [ -1, 1 ], variate = +0.0 ), second = 0.25 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 0, domain = [ -1, 1 ], variate = +0.0 ), second = 0.625 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 0, domain = [ -1, 1 ], variate = +0.0 ), second = 0.125 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 0, domain = [ -1, 1 ], variate = +1.0 ), second = 0.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 0, domain = [ -1, 1 ], variate = +1.0 ), second = 0.5 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 0, domain = [ -1, 1 ], variate = +1.0 ), second = 0.5 )

       def test_biunit_C1_quadratic_1st_deriv_at_nodes( self ):
              C = numpy.array( [ [ 1.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.5 ], [ 0.0, 0.0, 0.5 ] ] )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 1, domain = [ -1, 1 ], variate = -1.0 ), second = -1.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 1, domain = [ -1, 1 ], variate = -1.0 ), second = +1.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 1, domain = [ -1, 1 ], variate = -1.0 ), second = +0.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 1, domain = [ -1, 1 ], variate = +0.0 ), second = -0.5 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 1, domain = [ -1, 1 ], variate = +0.0 ), second = +0.25 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 1, domain = [ -1, 1 ], variate = +0.0 ), second = +0.25 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 1, domain = [ -1, 1 ], variate = +1.0 ), second = +0.0 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 1, domain = [ -1, 1 ], variate = +1.0 ), second = -0.5 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 1, domain = [ -1, 1 ], variate = +1.0 ), second = +0.5 )

       def test_biunit_C1_quadratic_2nd_deriv_at_nodes( self ):
              C = numpy.array( [ [ 1.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.5 ], [ 0.0, 0.0, 0.5 ] ] )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 2, domain = [ -1, 1 ], variate = -1.0 ), second = +0.50 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 2, domain = [ -1, 1 ], variate = -1.0 ), second = -0.75 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 2, domain = [ -1, 1 ], variate = -1.0 ), second = +0.25 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 2, domain = [ -1, 1 ], variate = +0.0 ), second = +0.50 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 2, domain = [ -1, 1 ], variate = +0.0 ), second = -0.75 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 2, domain = [ -1, 1 ], variate = +0.0 ), second = +0.25 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 0, deriv = 2, domain = [ -1, 1 ], variate = +1.0 ), second = +0.50 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 1, deriv = 2, domain = [ -1, 1 ], variate = +1.0 ), second = -0.75 )
              self.assertAlmostEqual( first = basis.evalSplineBasisDeriv1D( extraction_operator = C, basis_idx = 2, deriv = 2, domain = [ -1, 1 ], variate = +1.0 ), second = +0.25 )