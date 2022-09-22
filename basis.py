from string import digits
import unittest
import numpy
import sympy


def evalLegendreBasis1D( degree, variate):
       if degree == 0:
              val = 1.0
       elif degree == 1:
              val = variate
       else:
              i = degree - 1
              term_1 = i * evalLegendreBasis1D( degree = i - 1, variate = variate)
              term_2 = (2 * i + 1) * variate * evalLegendreBasis1D( degree = i, variate = variate)
              val = (term_2 - term_1) / (i + 1)
       return val

class Test_evalLegendreBasis1D( unittest.TestCase ):
       def test_linear( self ):
              for x in numpy.linspace ( -1, 1, 100):
                     self.assertAlmostEqual( first = evalLegendreBasis1D( degree = 1, variate = x), second = x, digits = 7)
              
       
