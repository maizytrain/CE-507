from evaluateSolutionAt import evaluateSolutionAt
from generate1DMesh import generateMesh1D
from evaluateSolutionAt import evalLagrangeBasis1D
import scipy
import numpy
from matplotlib import pyplot as plt

def computeFitError( target_fun, coeff, node_coords, connect, eval_basis ):
    num_elems = connect.shape[0]
    domain = [ min( node_coords ), max( node_coords ) ]
    abs_err_fun = lambda x : abs( target_fun( x ) - evaluateSolutionAt( x, coeff, node_coords, connect, eval_basis ) )
    fit_error, residual = scipy.integrate.quad( abs_err_fun, domain[0], domain[1], epsrel = 1e-12, limit = num_elems * 100 )
    return fit_error, residual


def evaluateErrorInFunction(target_fun, max2power, min, max, coeff):
    num_elems = 2 ** numpy.array( range( 0, max2power ) )
    errors = numpy.zeros(len(num_elems))
    residuals = numpy.zeros(len(num_elems))
    for i in range(len(num_elems)):
        node_coords, connect = generateMesh1D(min, max, num_elems[i], 4)
        errors[i], residuals[i] = computeFitError( target_fun, coeff, node_coords, connect, evalLagrangeBasis1D)

    plt.plot(num_elems, errors)
    plt.show()
    

coeff = numpy.array( [ 0.0, 1.0 ] )
#evaluateErrorInFunction(lambda x: x**2, 5, 0, 1, coeff)
node_coords, connect = generateMesh1D(0, 1, 16, 4)
computeFitError(lambda x: x**2, coeff, node_coords, connect, evalLagrangeBasis1D )