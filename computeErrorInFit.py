from evaluateSolutionAt import evaluateSolutionAt
from generate1DMesh import generateMesh1D
from evaluateSolutionAt import evalLagrangeBasis1D
import scipy
from scipy import integrate
import numpy
from matplotlib import pyplot as plt

def computeFitError( target_fun, coeff, node_coords, connect, eval_basis, degree ):
    num_elems = connect.shape[0]
    domain = [ min( node_coords ), max( node_coords ) ]
    abs_err_fun = lambda x : abs( target_fun( x ) - evaluateSolutionAt( x, coeff, node_coords, connect, eval_basis, degree ) )
    #print(evaluateSolutionAt( .5, coeff, node_coords, connect, eval_basis, degree ))
    fit_error, residual = integrate.quad( func = abs_err_fun, a = domain[0], b = domain[1], epsrel = 1e-12, limit = num_elems * 100 )
    return fit_error, residual


def evaluateErrorInFunction(target_fun, max2power, min, max, degree):
    num_elems = 2 ** numpy.array( range( 0, max2power ) )
    errors = numpy.zeros(len(num_elems))
    residuals = numpy.zeros(len(num_elems))
    for i in range(len(num_elems)):
        node_coords, connect = generateMesh1D(min, max, num_elems[i], degree)
        coeff = numpy.linspace(min, max, len(node_coords))
        for j in range(len(node_coords)):
            coeff[j] = target_fun(coeff[j])
        errors[i], residuals[i] = computeFitError( target_fun, coeff, node_coords, connect, evalLagrangeBasis1D, degree)

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.plot(num_elems, errors)
    plt.show()

evaluateErrorInFunction(lambda x: x**2, 5, 0, 1, 1)
evaluateErrorInFunction(lambda x: x**3, 5, 0, 1, 2)
evaluateErrorInFunction(lambda x: numpy.sin(numpy.pi * x), 5, -1, 1, 1)
evaluateErrorInFunction(lambda x: numpy.sin(numpy.pi * x), 5, 0, 1, 2)