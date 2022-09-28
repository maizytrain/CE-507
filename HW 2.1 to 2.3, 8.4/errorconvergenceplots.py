from string import digits
import unittest
import numpy
import sympy
import matplotlib.pyplot as plt

def taylorExpansion( fun, a, order ):
       x = list( fun.atoms( sympy.Symbol ) )[0]
       t = 0
       for i in range( 0, order + 1 ):
              df = sympy.diff( fun, x, i )
              term = ( df.subs( x, a ) / sympy.factorial( i ) ) * ( x - a )**i
              t += term
       return t

#Question 32 (otra vez lol)
j = sympy.Symbol('j')
#expr = sympy.integrate(sympy.sin(numpy.pi*j) - taylorExpansion( fun = expr, a = 0, order = 1), j)
x = numpy.linspace(-1,1,100)
values = numpy.linspace(-1,1,100)
orders = numpy.linspace(0,10,11)
maxError = numpy.linspace(0,1,11)
x_actual = numpy.linspace(-1,1,100)
x_0 = numpy.linspace(-1,1,100)
function = sympy.erfc(j)
for k in range(11):
    print("running" + str(k))
    expr = (function - taylorExpansion( fun = function, a = 0, order = k))
    print(expr)
    #exprpos = sympy.integrate(expr, j)
    #exprneg = sympy.integrate(-expr, j)
    for i in range(len(x)):
        #expr = sympy.integrate((sympy.sin(numpy.pi*j) - taylorExpansion( fun = sympy.sin(numpy.pi*j), a = 0, order = k)), (j,-1,x[i]))
        if (expr.subs(j,x[i]) >= 0):
            #print(str(expr.subs(j,x[i])) + "POSITIVE")
            values[i] = expr.subs(j,x[i])
        else:
            #print(str(expr.subs(j,x[i])) + "NEGATIVE")
            values[i] = -expr.subs(j,x[i])
        x_actual[i] = function.subs(j,x[i])
        x_0[i] = taylorExpansion( fun = function, a = 0, order = k).subs(j,x[i])
    maxError[k] = max(abs(values))
    plt.plot(x,values, label = "Error Order " + str(k))
    plt.plot(x,x_actual, label = "Actual")
    plt.plot(x,x_0, label = "Order " + str(k))
    plt.legend()
    plt.show()

fig = plt.figure()
ax = fig.add_subplot()
plt.plot(orders,maxError)
ax.set_yscale('log')
print(maxError)
plt.show()