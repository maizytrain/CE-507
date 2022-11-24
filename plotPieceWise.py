import numpy
import sympy
from evaluateSolutionAt import solveLagrangeBasis1D
from evaluateSolutionAt import evalLagrangeBasis1D
from generate1DMesh import generateMesh1D
from matplotlib import pyplot as plt

def plotPieceWise(func, min, max, num_elems, degree, resolution):
    node_coords, connect = generateMesh1D(min, max, num_elems, degree)
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
    for i in range(num_elems):
        start = int(connect[i][0])
        end = int(connect[i][len(connect[i])-1])
        sumexpr = 0
        for j in range(degree + 1):
            current = int(connect[i][j])
            expr = evalLagrangeBasis1D(node_coords[start], node_coords[end], degree, j) * func(node_coords[current])
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

# plotPieceWise(lambda x : x, -1, 1, 5, 2, 10)
# plotPieceWise(lambda x : x**2, -1, 1, 5, 2, 10)
# plotPieceWise(lambda x : numpy.sin(numpy.pi * x), -1, 1, 5, 2, 10)
# plotPieceWise(lambda x : numpy.exp(x), -1, 1, 5, 2, 10)
# plotPieceWise(lambda x : sympy.erfc(x), -1, 1, 5, 2, 10)