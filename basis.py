from string import digits
import unittest
import numpy
import sympy
import matplotlib.pyplot as plt

#Keeping this here literally just so I can see how to do some of this stuff
#def evalLegendreBasis1D( degree, variate):
#       if degree == 0:
#              val = 1.0
#       elif degree == 1:
#              val = variate
#       else:
#              i = degree - 1
#              term_1 = i * evalLegendreBasis1D( degree = i - 1, variate = variate)
#              term_2 = (2 * i + 1) * variate * evalLegendreBasis1D( degree = i, variate = variate)
#              val = (term_2 - term_1) / (i + 1)
#       return val

#class Test_evalLegendreBasis1D( unittest.TestCase ):
#       def test_linear( self ):
#              for x in numpy.linspace ( -1, 1, 100):
#                     self.assertAlmostEqual( first = evalLegendreBasis1D( degree = 1, variate = x), second = x, digits = 7)
              
       
#Question 32
def taylorExpansion( fun, a, order ):
       x = list( fun.atoms( sympy.Symbol ) )[0]
       t = 0
       for i in range( 0, order + 1 ):
              df = sympy.diff( fun, x, i )
              term = ( df.subs( x, a ) / sympy.factorial( i ) ) * ( x - a )**i
              t += term
       return t

class Test_evaltaylorExpansion( unittest.TestCase ):
       def test_linear( self ):
              j = sympy.Symbol('j')
              expr=sympy.sin(j)
              for i in numpy.linspace ( -1, 1, 3):
                     #print(taylorExpansion( fun = expr, a = i, order = 10).subs(j,i))
                     self.assertAlmostEqual(taylorExpansion( fun = expr, a = i, order = 10).subs(j,i), expr.subs(j,i), places = 4)

#Question 32 (2? I guess?)

j = sympy.Symbol('j')
expr=sympy.sin(j*numpy.pi)
x = numpy.linspace(-1,1,100)
x_actual = numpy.linspace(-1,1,100)
x_0=numpy.linspace(-1,1,100)
x_1=numpy.linspace(-1,1,100)
x_3=numpy.linspace(-1,1,100)
x_5=numpy.linspace(-1,1,100)
x_7=numpy.linspace(-1,1,100)
#replacements = [(j, x) for i in range(x)]
for i in range(len(x)):
       x_actual[i] = expr.subs(j,x[i])
       x_0[i] = taylorExpansion( fun = expr, a = 0, order = 0).subs(j,x[i])
       x_1[i] = taylorExpansion( fun = expr, a = 0, order = 1).subs(j,x[i])
       x_3[i] = taylorExpansion( fun = expr, a = 0, order = 3).subs(j,x[i])
       x_5[i] = taylorExpansion( fun = expr, a = 0, order = 5).subs(j,x[i])
       x_7[i] = taylorExpansion( fun = expr, a = 0, order = 7).subs(j,x[i])
plt.plot(x,x_actual, label = "actual")
plt.plot(x,x_0, label = "0th Order")
plt.plot(x,x_1, label = "1st Order")
plt.plot(x,x_3, label = "3rd Order")
plt.plot(x,x_5, label = "5th Order")
plt.plot(x,x_7, label = "7th Order")
plt.legend()
plt.show()

j = sympy.Symbol('j')
expr=sympy.exp(j)
x = numpy.linspace(-1,1,100)
x_actual = numpy.linspace(-1,1,100)
x_0=numpy.linspace(-1,1,100)
x_1=numpy.linspace(-1,1,100)
x_2=numpy.linspace(-1,1,100)
x_3=numpy.linspace(-1,1,100)
x_4=numpy.linspace(-1,1,100)
#replacements = [(j, x) for i in range(x)]
for i in range(len(x)):
       x_actual[i] = expr.subs(j,x[i])
       x_0[i] = taylorExpansion( fun = expr, a = 0, order = 0).subs(j,x[i])
       x_1[i] = taylorExpansion( fun = expr, a = 0, order = 1).subs(j,x[i])
       x_2[i] = taylorExpansion( fun = expr, a = 0, order = 2).subs(j,x[i])
       x_3[i] = taylorExpansion( fun = expr, a = 0, order = 3).subs(j,x[i])
       x_4[i] = taylorExpansion( fun = expr, a = 0, order = 4).subs(j,x[i])
plt.plot(x,x_actual, label = "actual")
plt.plot(x,x_0, label = "0th Order")
plt.plot(x,x_1, label = "1st Order")
plt.plot(x,x_2, label = "2rd Order")
plt.plot(x,x_3, label = "3th Order")
plt.plot(x,x_4, label = "4th Order")
plt.legend()
plt.show()

j = sympy.Symbol('j')
expr=sympy.erfc(j)
x = numpy.linspace(-2,2,100)
x_actual = numpy.linspace(-2,2,100)
x_0=numpy.linspace(-2,2,100)
x_1=numpy.linspace(-2,2,100)
x_3=numpy.linspace(-2,2,100)
x_5=numpy.linspace(-2,2,100)
x_7=numpy.linspace(-2,2,100)
x_9=numpy.linspace(-2,2,100)
x_11=numpy.linspace(-2,2,100)
#replacements = [(j, x) for i in range(x)]
for i in range(len(x)):
       x_actual[i] = expr.subs(j,x[i])
       x_0[i] = taylorExpansion( fun = expr, a = 0, order = 0).subs(j,x[i])
       x_1[i] = taylorExpansion( fun = expr, a = 0, order = 1).subs(j,x[i])
       x_3[i] = taylorExpansion( fun = expr, a = 0, order = 3).subs(j,x[i])
       x_5[i] = taylorExpansion( fun = expr, a = 0, order = 5).subs(j,x[i])
       x_7[i] = taylorExpansion( fun = expr, a = 0, order = 7).subs(j,x[i])
       x_9[i] = taylorExpansion( fun = expr, a = 0, order = 9).subs(j,x[i])
       x_11[i] = taylorExpansion( fun = expr, a = 0, order = 11).subs(j,x[i])
plt.plot(x,x_actual, label = "actual")
plt.plot(x,x_0, label = "0th Order")
plt.plot(x,x_1, label = "1st Order")
plt.plot(x,x_3, label = "3rd Order")
plt.plot(x,x_5, label = "5th Order")
plt.plot(x,x_7, label = "7th Order")
plt.plot(x,x_9, label = "9th Order")
plt.plot(x,x_11, label = "11th Order")
plt.legend()
plt.show()
                     