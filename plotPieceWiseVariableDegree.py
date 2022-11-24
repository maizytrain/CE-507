import numpy
import sympy
from evaluateSolutionAt import solveLagrangeBasis1D
from evaluateSolutionAt import evalLagrangeBasis1D
from generate1DMesh import generateMesh1D
from plotPieceWise import plotPieceWise
from matplotlib import pyplot as plt

def plotPieceWiseVariableDegree(func, min, max, num_elems, degreeList, resolution):
    xList = numpy.linspace(min, max, num_elems+1)

    for l in range(len(degreeList)):
        node_coords, connect = generateMesh1D(xList[l], xList[l+1], 1, degreeList[l])
        zeros = numpy.zeros(len(node_coords))
        plt.scatter(node_coords, zeros, color="red")

        for i in range(len(connect)):
            for j in range(len(connect[i])):
                ys = numpy.zeros(len(connect[i]))
                xs = numpy.zeros(len(connect[i]))
                for k in range(len(xs)):
                    index = int(connect[i][k])
                    xs[k] = node_coords[index]
                plt.plot(xs, ys, color="red")

        zeta = sympy.symbols('zeta')
        for i in range(1):
            start = int(connect[i][0])
            end = int(connect[i][len(connect[i])-1])
            sumexpr = 0
            for j in range(degreeList[l] + 1):
                current = int(connect[i][j])
                expr = evalLagrangeBasis1D(node_coords[start], node_coords[end], degreeList[l], j) * func(node_coords[current])
                sumexpr += expr
                xs = numpy.linspace(node_coords[start], node_coords[end], resolution)
                ys = numpy.zeros(len(xs))
                for k in range(len(xs)):
                    ys[k] = expr.subs(zeta, xs[k])
                plt.plot(xs,ys)
            xs = numpy.linspace(node_coords[start], node_coords[end], resolution)
            ys = numpy.zeros(len(xs))
            for j in range(len(xs)):
                ys[j] = sumexpr.subs(zeta, xs[j])
            plt.plot(xs,ys,color="blue")


    plt.show()

degreeList = [1,2,3,4]
plotPieceWiseVariableDegree(lambda x : x, -1, 1, 4, degreeList, 20)
plotPieceWiseVariableDegree(lambda x : x**2, -1, 1, 4, degreeList, 10)
plotPieceWiseVariableDegree(lambda x : numpy.sin(numpy.pi * x), -1, 1, 4, degreeList, 10)
plotPieceWiseVariableDegree(lambda x : numpy.exp(x), -1, 1, 4, degreeList, 10)
plotPieceWiseVariableDegree(lambda x : sympy.erfc(x), -1, 1, 4, degreeList, 10)