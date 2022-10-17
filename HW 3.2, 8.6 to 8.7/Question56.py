import numpy
import numpy as np
import sympy
import sympy as sp
import matplotlib.pyplot as plt
from Question52 import generateMesh1D
#from CE507 import computeNewtonCotesQuadrature

def createLagrangeBasis1D( degree, basis_idx):
    zeta = sympy.symbols('zeta')
    values = numpy.linspace(-1,1, degree+1)
    expr = 1

    for i in range(degree+1):
        if (i != basis_idx):
            expr *= (zeta - values[i])/(values[basis_idx] - values[i])

    #print(values)
    #print(expr)
    
    return expr

def computeNewtonCotesQuadrature( fun, num_points, degree, lowerBound, upperBound ):
    calcMat = np.zeros([degree+1,degree+1])
    for i in range(degree+1):
        expr = sp.Poly(sp.expand(sp.simplify(createLagrangeBasis1D( degree = degree, basis_idx = i))))
        coeffs = expr.all_coeffs()
        for j in range(degree+1):
            calcMat[i,j] = coeffs[j]

    calcMat = np.transpose(calcMat)

    xBounds = np.linspace(lowerBound,upperBound,num_points+1)
    finals = 0
    zeta = sp.symbols('zeta')

    xValues = np.linspace(lowerBound,upperBound,50*num_points)
    yActual = np.zeros(50*num_points)
    yTest = np.zeros(50*num_points)
    for i in range(num_points):
        eqExpr = calcSingleNewtonCotes(fun = fun, degree = degree, lowerBound = xBounds[i], upperBound = xBounds[i+1], calcMat = calcMat)
        #finals += sp.integrate(eqExpr, (zeta, xBounds[i], xBounds[i+1]))
        for j in range(50):
            yTest[j+50*i] = eqExpr.subs(zeta, xValues[j+50*i])

    for i in range(50*num_points):
        yActual[i] = fun.subs(x, xValues[i])

    plt.plot(xValues,yActual, label = "Actual Values")
    plt.plot(xValues,yTest, label = "Approximated Values")
    plt.legend()
    plt.show()

    

    return yActual


def calcSingleNewtonCotes( fun, degree, lowerBound, upperBound, calcMat ):
    xValues = np.linspace(lowerBound,upperBound,degree+1)
    fValues = np.zeros(degree+1)
    x = sp.symbols('x')
    for i in range(degree+1):
        fValues[i] = fun.subs(x,xValues[i])
    
    iMat = np.zeros([degree+1, degree+1])
    for i in range(degree + 1):
        for j in range(degree + 1):
            iMat[i, j] = xValues[i]**(degree - j)

    polyValues = np.linalg.solve(iMat,fValues)
    lValues = np.linalg.solve(calcMat,polyValues)

    zeta = sp.symbols('zeta')
    lExpr = zeta*0
    for i in range(degree+1):
        lExpr += lValues[i] * createLagrangeBasis1D( degree = degree, basis_idx = i)

    sExpr = (sp.simplify(lExpr))

    return sExpr

#xs = np.linspace(-1,1,3)
#for i in range(3):
#    createLagrangeBasis1D(3, 2)


zeta = sp.symbols('zeta')
x = sp.symbols('x')

print(computeNewtonCotesQuadrature(fun = sp.erfc(x), num_points = 3, degree = 3, lowerBound = -1, upperBound = 1))
