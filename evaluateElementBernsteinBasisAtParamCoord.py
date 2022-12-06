import numpy
import matplotlib
from matplotlib import pyplot as plt
import readBEXTJSON as read
import basis

def evaluateElementBernsteinBasisAtParamCoord( uspline, elem_id, param_coord ):
    elem_degree = read.getElementDegree(uspline, elem_id)
    elem_bernstein_basis = numpy.zeros( elem_degree + 1 )
    for n in range( 0, elem_degree + 1 ):
        elem_bernstein_basis[n] = basis.evaluateBernsteinBasis1D(param_coord, elem_degree, n)
    return elem_bernstein_basis

def evaluateElementSplineBasisAtParamCoord( uspline, elem_id, param_coord ):
    elem_ext_operator = read.getElementExtractionOperator(uspline, elem_id) # Get the extraction operator of the element
    #print(elem_ext_operator)
    elem_bernstein_basis = evaluateElementBernsteinBasisAtParamCoord( uspline, elem_id, param_coord )
    elem_spline_basis = numpy.matmul(elem_ext_operator, elem_bernstein_basis)
    #print("param_coord =" + str(param_coord) + ", basis = " + str(elem_spline_basis))
    return elem_spline_basis 

def plotUsplineBasis( uspline, color_by ):
    num_pts = 100
    xi = numpy.linspace( 0, 1, num_pts )
    num_elems = read.getNumElems(uspline)
    for elem_idx in range( 0, num_elems ):
        elem_id = read.elemIdFromElemIdx(uspline, elem_idx)
        elem_domain = read.getElementDomain(uspline, elem_id)
        elem_node_ids = read.getSplineNodes(uspline)
        elem_degree = read.getElementDegree(uspline, elem_id)
        x = numpy.linspace( elem_domain[0], elem_domain[1], num_pts )
        y = numpy.zeros( shape = ( elem_degree + 1, num_pts ) )
        mysum = numpy.zeros( num_pts )
        variate = numpy.linspace(-1, 1, num_pts)
        for i in range( 0, num_pts ):
            y[:,i] = evaluateElementSplineBasisAtParamCoord(uspline, elem_id, variate[i])# Evaluate the current element’s spline basis at the current coordinate
            mysum[i] = sum(y[:,i])
        # Do plotting for the current element
        for n in range( 0, elem_degree + 1 ):
            if color_by == "element":
                color = getLineColor( elem_idx )
            elif color_by == "node":
                color = getLineColor( elem_node_ids[n] )
            matplotlib.pyplot.plot( x, y[n,:], color = color )
            plt.plot(x, mysum, color = "blue")
    plt.show()

def getLineColor( idx ):
    colors = list( matplotlib.colors.TABLEAU_COLORS.keys() )
    num_colors = len( colors )
    mod_idx = idx % num_colors
    return matplotlib.colors.TABLEAU_COLORS[ colors[ mod_idx ] ]
    
#uspline = read.readBEXT( "quadratic_bspline.json" )
#plotUsplineBasis( uspline, "element" )
#plotUsplineBasis( uspline, "node" )