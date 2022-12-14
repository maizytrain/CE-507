import unittest
import basis
import math




class Test_evalBernsteinBasisDeriv( unittest.TestCase ):
       def test_constant_at_nodes( self ):
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 0, basis_idx = 0, deriv = 0, domain = [0, 1], variate = 0.0 ), second = 1.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 0, basis_idx = 0, deriv = 0, domain = [0, 1], variate = 1.0 ), second = 1.0, delta = 1e-12 )

       def test_constant_1st_deriv_at_nodes( self ):
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 0, basis_idx = 0, deriv = 1, domain = [0, 1], variate = 0.0 ), second = 0.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 0, basis_idx = 0, deriv = 1, domain = [0, 1], variate = 1.0 ), second = 0.0, delta = 1e-12 )

       def test_constant_2nd_deriv_at_nodes( self ):
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 0, basis_idx = 0, deriv = 2, domain = [0, 1], variate = 0.0 ), second = 0.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 0, basis_idx = 0, deriv = 2, domain = [0, 1], variate = 1.0 ), second = 0.0, delta = 1e-12 )

       def test_linear_at_nodes( self ):
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 0, deriv = 0, domain = [0, 1], variate = 0.0 ), second = 1.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 0, deriv = 0, domain = [0, 1], variate = 1.0 ), second = 0.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 1, deriv = 0, domain = [0, 1], variate = 0.0 ), second = 0.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 1, deriv = 0, domain = [0, 1], variate = 1.0 ), second = 1.0, delta = 1e-12 )

       def test_linear_at_gauss_pts( self ):
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 0, deriv = 0, domain = [0, 1], variate =  0.5 ), second = 0.5, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 1, deriv = 0, domain = [0, 1], variate =  0.5 ), second = 0.5, delta = 1e-12 )

       def test_quadratic_at_nodes( self ):
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 0, deriv = 0, domain = [0, 1], variate = 0.0 ), second = 1.00, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 0, deriv = 0, domain = [0, 1], variate = 0.5 ), second = 0.25, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 0, deriv = 0, domain = [0, 1], variate = 1.0 ), second = 0.00, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 1, deriv = 0, domain = [0, 1], variate = 0.0 ), second = 0.00, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 1, deriv = 0, domain = [0, 1], variate = 0.5 ), second = 0.50, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 1, deriv = 0, domain = [0, 1], variate = 1.0 ), second = 0.00, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 2, deriv = 0, domain = [0, 1], variate = 0.0 ), second = 0.00, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 2, deriv = 0, domain = [0, 1], variate = 0.5 ), second = 0.25, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 2, deriv = 0, domain = [0, 1], variate = 1.0 ), second = 1.00, delta = 1e-12 )

       def test_quadratic_at_gauss_pts( self ):
              x = [ -1.0 / math.sqrt(3.0) , 1.0 / math.sqrt(3.0) ]
              x = [ basis.affine_mapping_1D( [-1, 1], [0, 1], xi ) for xi in x ]
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 0, deriv = 0, domain = [0, 1], variate = x[0] ), second = 0.62200846792814620, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 0, deriv = 0, domain = [0, 1], variate = x[1] ), second = 0.04465819873852045, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 1, deriv = 0, domain = [0, 1], variate = x[0] ), second = 0.33333333333333333, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 1, deriv = 0, domain = [0, 1], variate = x[1] ), second = 0.33333333333333333, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 2, deriv = 0, domain = [0, 1], variate = x[0] ), second = 0.04465819873852045, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 2, deriv = 0, domain = [0, 1], variate = x[1] ), second = 0.62200846792814620, delta = 1e-12 )

       def test_linear_1st_deriv_at_nodes( self ):
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 0, deriv = 1, domain = [0, 1], variate = 0.0 ), second = -1.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 0, deriv = 1, domain = [0, 1], variate = 1.0 ), second = -1.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 1, deriv = 1, domain = [0, 1], variate = 0.0 ), second = +1.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 1, deriv = 1, domain = [0, 1], variate = 1.0 ), second = +1.0, delta = 1e-12 )

       def test_linear_1st_deriv_at_gauss_pts( self ):
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 0, deriv = 1, domain = [0, 1], variate = 0.5 ), second = -1.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 1, deriv = 1, domain = [0, 1], variate = 0.5 ), second = +1.0, delta = 1e-12 )

       def test_linear_2nd_deriv_at_nodes( self ):
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 0, deriv = 2, domain = [0, 1], variate = 0.0 ), second = 0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 0, deriv = 2, domain = [0, 1], variate = 1.0 ), second = 0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 1, deriv = 2, domain = [0, 1], variate = 0.0 ), second = 0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 1, deriv = 2, domain = [0, 1], variate = 1.0 ), second = 0, delta = 1e-12 )

       def test_linear_2nd_deriv_at_gauss_pts( self ):
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 0, deriv = 2, domain = [0, 1], variate = 0.5 ), second = 0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 1, basis_idx = 1, deriv = 2, domain = [0, 1], variate = 0.5 ), second = 0, delta = 1e-12 )

       def test_quadratic_1st_deriv_at_nodes( self ):
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 0, deriv = 1, domain = [0, 1], variate = 0.0 ), second = -2.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 0, deriv = 1, domain = [0, 1], variate = 0.5 ), second = -1.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 0, deriv = 1, domain = [0, 1], variate = 1.0 ), second =  0.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 1, deriv = 1, domain = [0, 1], variate = 0.0 ), second = +2.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 1, deriv = 1, domain = [0, 1], variate = 0.5 ), second =  0.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 1, deriv = 1, domain = [0, 1], variate = 1.0 ), second = -2.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 2, deriv = 1, domain = [0, 1], variate = 0.0 ), second =  0.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 2, deriv = 1, domain = [0, 1], variate = 0.5 ), second =  1.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 2, deriv = 1, domain = [0, 1], variate = 1.0 ), second = +2.0, delta = 1e-12 )

       def test_quadratic_1st_deriv_at_gauss_pts( self ):
              x = [ -1.0 / math.sqrt(3.0) , 1.0 / math.sqrt(3.0) ]
              x = [ basis.affine_mapping_1D( [-1, 1], [0, 1], xi ) for xi in x ]
              print(x)
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 0, deriv = 1, domain = [0, 1], variate = x[0] ), second = -1.0 - 1/( math.sqrt(3) ), delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 0, deriv = 1, domain = [0, 1], variate = x[1] ), second = -1.0 + 1/( math.sqrt(3) ), delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 1, deriv = 1, domain = [0, 1], variate = x[0] ), second = +2.0 / math.sqrt(3), delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 1, deriv = 1, domain = [0, 1], variate = x[1] ), second = -2.0 / math.sqrt(3), delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 2, deriv = 1, domain = [0, 1], variate = x[0] ), second = +1.0 - 1/( math.sqrt(3) ), delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 2, deriv = 1, domain = [0, 1], variate = x[1] ), second = +1.0 + 1/( math.sqrt(3) ), delta = 1e-12 )

       def test_quadratic_2nd_deriv_at_nodes( self ):
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 0, deriv = 1, domain = [0, 1], variate = 0.0 ), second = -2.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 0, deriv = 1, domain = [0, 1], variate = 0.5 ), second = -1.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 0, deriv = 1, domain = [0, 1], variate = 1.0 ), second =  0.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 1, deriv = 1, domain = [0, 1], variate = 0.0 ), second = +2.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 1, deriv = 1, domain = [0, 1], variate = 0.5 ), second =  0.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 1, deriv = 1, domain = [0, 1], variate = 1.0 ), second = -2.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 2, deriv = 1, domain = [0, 1], variate = 0.0 ), second =  0.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 2, deriv = 1, domain = [0, 1], variate = 0.5 ), second =  1.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 2, deriv = 1, domain = [0, 1], variate = 1.0 ), second = +2.0, delta = 1e-12 )

       def test_quadratic_2nd_deriv_at_gauss_pts( self ):
              x = [ -1.0 / math.sqrt(3.0) , 1.0 / math.sqrt(3.0) ]
              x = [ basis.affine_mapping_1D( [-1, 1], [0, 1], xi ) for xi in x ]
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 0, deriv = 2, domain = [0, 1], variate = x[0] ), second = +2.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 0, deriv = 2, domain = [0, 1], variate = x[1] ), second = +2.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 1, deriv = 2, domain = [0, 1], variate = x[0] ), second = -4.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 1, deriv = 2, domain = [0, 1], variate = x[1] ), second = -4.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 2, deriv = 2, domain = [0, 1], variate = x[0] ), second = +2.0, delta = 1e-12 )
              self.assertAlmostEqual( first = basis.evalBernsteinBasisDeriv( degree = 2, basis_idx = 2, deriv = 2, domain = [0, 1], variate = x[1] ), second = +2.0, delta = 1e-12 )