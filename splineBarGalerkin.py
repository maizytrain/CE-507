import numpy
import scipy
from scipy import integrate
import matplotlib
import matplotlib.pyplot as plt
import sympy
import readBEXTJSON as bext
import stiffnessMatrix
import basis
import math
from sympy.integrals import quadrature



## MAIN CODE
def computeSolution( problem, uspline_bext ):
    K = assembleStiffnessMatrix(problem, uspline_bext)
    F = assembleForceVector(problem, uspline_bext)
    K, F = applyDisplacement(problem, K, F, uspline_bext)
    coeff = numpy.matmul(numpy.linalg.inv(K), F)
    coeff = assembleSolution(coeff, problem, uspline_bext)
    return coeff

def assembleSolution( coeff, problem, uspline_bext ):
    # Your code goes here
    return coeff

def applyDisplacement( problem, stiffness_matrix, force_vector, uspline_bext ):
    # Your code goes here
    idx = bext.getElementIdContainingPoint(uspline_bext, problem["displacement"]["position"])
    print("idx = ", idx)
    return stiffness_matrix, force_vector

def applyTraction( problem, force_vector, uspline_bext ):
    zeta = sympy.symbols('zeta')
    x = sympy.symbols('x')
    elem_id = bext.getElementIdContainingPoint(uspline_bext, problem["traction"]["position"])
    elem_degree = bext.getElementDegree(uspline_bext, elem_id)
    elem_domain = bext.getElementDomain(uspline_bext, elem_id)
    elem_ext_oper = bext.getElementExtractionOperator(uspline_bext, elem_id)
    elem_bernstein_basis = [None] * (elem_degree + 1)
    change_expr = basis.affine_mapping_1D(elem_domain, [0, 1], x)

    for n in range( 0, elem_degree + 1 ):
        elem_bernstein_basis[n] = basis.createBernsteinBasis1DDeriv(elem_degree, n, 1)
        elem_bernstein_basis[n] = elem_bernstein_basis[n].subs(zeta, change_expr)
    elem_spline_basis = numpy.matmul(elem_ext_oper, elem_bernstein_basis)

    for i in range(elem_degree + 1):
        idx_a = bext.getElementNodeIds(uspline_bext, elem_id)[i]
        # print(elem_spline_basis[i] * problem["traction"]["value"])
        # print(idx_a)
        force_vector[idx_a] += elem_spline_basis[i] * problem["traction"]["value"]


    return force_vector

def evaluateConstitutiveModel( problem ):
    # Your code goes here
    return

def assembleStiffnessMatrix( problem, uspline_bext ):
    stiffness_matrix = stiffnessMatrix.assembleStiffnessMatrix(problem, uspline_bext)
    return stiffness_matrix

def assembleForceVector( problem, uspline_bext ):
    zeta = sympy.symbols('zeta')
    x = sympy.symbols('x')
    num_elems = bext.getNumElems(uspline_bext)
    num_nodes = bext.getNumNodes(uspline_bext)
    force_vector = numpy.zeros(num_nodes)

    for i in range(num_elems):
        elem_id = bext.elemIdFromElemIdx(uspline_bext, i)
        elem_degree = bext.getElementDegree(uspline_bext, elem_id)
        elem_domain = bext.getElementDomain(uspline_bext, elem_id)
        elem_ext_oper = bext.getElementExtractionOperator(uspline_bext, elem_id)
        elem_bernstein_basis = [None] * (elem_degree + 1)
        change_expr = basis.affine_mapping_1D(elem_domain, [0, 1], x)

        for n in range( 0, elem_degree + 1 ):
            elem_bernstein_basis[n] = basis.createBernsteinBasis1DDeriv(elem_degree, n, 1)
            elem_bernstein_basis[n] = elem_bernstein_basis[n].subs(zeta, change_expr)
        elem_spline_basis = numpy.matmul(elem_ext_oper, elem_bernstein_basis)

        for j in range(elem_degree + 1):
            idx_a = bext.getElementNodeIds(uspline_bext, elem_id)[j]
            #print(elem_bernstein_basis)
            force_vector[idx_a] += sympy.integrate(elem_spline_basis[j] * problem["body_force"], (x, elem_domain[0], elem_domain[1]))
            #print(K[idx_a, idx_b])

    force_vector = applyTraction(problem, force_vector, uspline_bext)

    return force_vector

problem = { "elastic_modulus": 100,
                       "area": 0.01,
                       "length": 1.0,
                       "traction": { "value": 1e-3, "position": 1.0 },
                       "displacement": { "value": 0.0, "position": 0.0 },
                       "body_force": 1e-3 }
uspline_bext = bext.readBEXT( "test_simple.json" )

print(computeSolution( problem = problem, uspline_bext = uspline_bext ))
#print(assembleForceVector(problem, uspline_bext))

## UTILITY CODE
def evaluateSolutionAt( x, coeff, uspline_bext ):
    elem_id = bext.getElementIdContainingPoint( uspline_bext, x )
    elem_nodes = bext.getElementNodeIds( uspline_bext, elem_id )
    elem_domain = bext.getElementDomain( uspline_bext, elem_id )
    elem_degree = bext.getElementDegree( uspline_bext, elem_id )
    elem_extraction_operator = bext.getElementExtractionOperator( uspline_bext, elem_id )
    sol = 0.0
    for n in range( 0, len( elem_nodes ) ):
        curr_node = elem_nodes[n]
        sol += coeff[curr_node] * basis.evalSplineBasis1D( extraction_operator = elem_extraction_operator, basis_idx = n, domain = elem_domain, variate = x )
    return sol

def computeElementFitError( problem, coeff, uspline_bext, elem_id ):
    domain = bext.getDomain( uspline_bext )
    elem_domain = bext.getElementDomain( uspline_bext, elem_id )
    elem_degree = bext.getElementDegree( uspline_bext, elem_id )
    num_qp = int( numpy.ceil( ( 2*(elem_degree - 1) + 1 ) / 2.0 ) + 1 )
    abs_err_fun = lambda x : abs( evaluateExactSolutionAt( problem, basis.affine_mapping_1D( [-1, 1], elem_domain, x ) ) - evaluateSolutionAt( basis.affine_mapping_1D( [-1, 1], elem_domain, x ), coeff, uspline_bext ) )
    abs_error = quadrature.quad( abs_err_fun, elem_domain, num_qp )
    return abs_error

def computeFitError( problem, coeff, uspline_bext ):
    num_elems = bext.getNumElems( uspline_bext )
    abs_error = 0.0
    for elem_idx in range( 0, num_elems ):
        elem_id = bext.elemIdFromElemIdx( uspline_bext, elem_idx )
        abs_error += computeElementFitError( problem, coeff, uspline_bext, elem_id )
    domain = bext.getDomain( uspline_bext )
    target_fun_norm, _ = scipy.integrate.quad( lambda x: abs( evaluateExactSolutionAt( problem, x ) ), domain[0], domain[1], epsrel = 1e-12, limit = num_elems * 100 )
    rel_error = abs_error / target_fun_norm
    return abs_error, rel_error

def plotCompareGoldTestSolution( gold_coeff, test_coeff, uspline_bext ):
    domain = bext.getDomain( uspline_bext )
    x = numpy.linspace( domain[0], domain[1], 1000 )
    yg = numpy.zeros( 1000 )
    yt = numpy.zeros( 1000 )
    for i in range(0, len(x) ):
        yg[i] = evaluateSolutionAt( x[i], test_coeff, uspline_bext )
        yt[i] = evaluateSolutionAt( x[i], gold_coeff, uspline_bext )
    plt.plot( x, yg )
    plt.plot( x, yt )
    plt.show()

def plotCompareFunToExactSolution( problem, test_coeff, uspline_bext ):
    domain = bext.getDomain( uspline_bext )
    x = numpy.linspace( domain[0], domain[1], 1000 )
    ya = numpy.zeros( 1000 )
    ye = numpy.zeros( 1000 )
    for i in range(0, len(x) ):
        ya[i] = evaluateSolutionAt( x[i], test_coeff, uspline_bext )
        ye[i] = evaluateExactSolutionAt( problem, x[i] )
    plt.plot( x, ya )
    plt.plot( x, ye )
    plt.show()

def computeConvergenceRate( num_entities, qoi ):
    def func( x, a, b, c ):
        return a * numpy.power( x, b ) + c
    fit = scipy.optimize.curve_fit(func, num_entities, qoi, method='trf', bounds = ([-numpy.inf, -numpy.inf, -numpy.inf ], [numpy.inf, 0.0, numpy.inf]) )
    a,b,c = fit[0]
    return b

def plotSolution( sol_coeff, uspline_bext ):
    domain = bext.getDomain( uspline_bext )
    x = numpy.linspace( domain[0], domain[1], 1000 )
    y = numpy.zeros( 1000 )
    for i in range(0, len(x) ):
        y[i] = evaluateSolutionAt( x[i], sol_coeff, uspline_bext )
    plt.plot( x, y )
    plt.plot( bext.getSplineNodes( uspline_bext )[:,0], sol_coeff, color = "k", marker = "o", markerfacecolor = "k" )
    plt.show()

def evaluateExactSolutionAt( problem, x ):
    term_1 = problem[ "traction" ][ "value" ] / evaluateConstitutiveModel( problem ) * x
    term_2 = problem[ "displacement" ][ "value" ]
    term_3 =  ( ( problem[ "length" ]**2.0 * problem[ "body_force" ] / 2 ) / evaluateConstitutiveModel( problem ) ) - ( ( ( problem[ "length" ] - x )**2.0 * problem[ "body_force" ] / 2 ) / evaluateConstitutiveModel( problem ) )
    sol = term_1 + term_2 + term_3
    return sol

def plotExactSolution( problem ):
    domain = [0, problem[ "length" ] ]
    x = numpy.linspace( domain[0], domain[1], 1000 )
    y = numpy.zeros( 1000 )
    for i in range(0, len(x) ):
        y[i] = evaluateExactSolutionAt( problem, x[i] )
    plt.plot( x, y )
    plt.show()