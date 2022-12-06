import unittest
import readBEXTJSON as bext
import sympy
import numpy
import evaluateElementBernsteinBasisAtParamCoord
import basis
import math

def assembleStiffnessMatrix(problem, uspline_bext):
    zeta = sympy.symbols('zeta')
    x = sympy.symbols('x')
    num_elems = bext.getNumElems(uspline_bext)
    num_nodes = bext.getNumNodes(uspline_bext)
    K = numpy.zeros([num_nodes, num_nodes])
    for i in range(num_elems):
        elem_id = bext.elemIdFromElemIdx(uspline_bext, i)
        elem_degree = bext.getElementDegree(uspline_bext, elem_id)
        elem_domain = bext.getElementDomain(uspline_bext, elem_id)
        elem_ext_oper = bext.getElementExtractionOperator(uspline_bext, elem_id)
        elem_bernstein_basis_deriv = [None] * (elem_degree + 1)
        change_expr = basis.affine_mapping_1D(elem_domain, [0, 1], x)
        for n in range( 0, elem_degree + 1 ):
            elem_bernstein_basis_deriv[n] = basis.createBernsteinBasis1DDeriv(elem_degree, n, 1)
            elem_bernstein_basis_deriv[n] = elem_bernstein_basis_deriv[n].subs(zeta, change_expr) / math.pow((elem_domain[1]-elem_domain[0]), 1)
        # print(elem_bernstein_basis)
        # sympy.plot(elem_bernstein_basis[0], elem_bernstein_basis[1], elem_bernstein_basis[2], (x, elem_domain[0], elem_domain[1]))
        elem_spline_basis_deriv = numpy.matmul(elem_ext_oper, elem_bernstein_basis_deriv)
        # print(elem_spline_basis)
        # sympy.plot(elem_spline_basis[0], elem_spline_basis[1], elem_spline_basis[2], (x, elem_domain[0], elem_domain[1]))
        for j in range(elem_degree + 1):
            idx_a = bext.getElementNodeIds(uspline_bext, elem_id)[j]
            #print(elem_bernstein_basis)
            for k in range(elem_degree + 1):
                idx_b = bext.getElementNodeIds(uspline_bext, elem_id)[k]
                print(idx_a, idx_b)
                K[idx_a, idx_b] += sympy.integrate(elem_spline_basis_deriv[j] * problem["elastic_modulus"] * problem["area"] * elem_spline_basis_deriv[k], (x, elem_domain[0], elem_domain[1]))
                #print(K[idx_a, idx_b])
    return K

problem = { "elastic_modulus": 100,
                     "area": 0.01,
                     "length": 1.0,
                     "traction": { "value": 1e-3, "position": 1.0 },
                     "displacement": { "value": 0.0, "position": 0.0 },
                     "body_force": 0.0 }
uspline_bext = bext.readBEXT( "test_one_linear_C0_element.json" )
print(assembleStiffnessMatrix( problem = problem, uspline_bext = uspline_bext ))



class test_assembleStressMatrix( unittest.TestCase ):
       def test_one_linear_C0_element( self ):
              problem = { "elastic_modulus": 100,
                     "area": 0.01,
                     "length": 1.0,
                     "traction": { "value": 1e-3, "position": 1.0 },
                     "displacement": { "value": 0.0, "position": 0.0 },
                     "body_force": 0.0 }
              #spline_space = { "domain": [0, problem[ "length" ]], "degree": [ 1 ], "continuity": [ -1, -1 ] }
              #uspline.make_uspline_mesh( spline_space, "test_one_linear_C0_element" )
              uspline_bext = bext.readBEXT( "test_one_linear_C0_element.json" )
              test_stiffness_matrix = assembleStiffnessMatrix( problem = problem, uspline_bext = uspline_bext )
              gold_stiffness_matrix = numpy.array( [ [ 1.0, -1.0 ], [ -1.0, 1.0 ] ] )
              self.assertTrue( numpy.allclose( test_stiffness_matrix, gold_stiffness_matrix ) )

       def test_two_linear_C0_element( self ):
              problem = { "elastic_modulus": 100,
                     "area": 0.01,
                     "length": 1.0,
                     "traction": { "value": 1e-3, "position": 1.0 },
                     "displacement": { "value": 0.0, "position": 0.0 },
                     "body_force": 0.0 }
              #spline_space = { "domain": [0, problem[ "length" ]], "degree": [ 1, 1 ], "continuity": [ -1, 0, -1 ] }
              #uspline.make_uspline_mesh( spline_space, "test_two_linear_C0_element" )
              uspline_bext = bext.readBEXT( "test_two_linear_C0_element.json" )
              test_stiffness_matrix = assembleStiffnessMatrix( problem = problem, uspline_bext = uspline_bext )
              gold_stiffness_matrix = numpy.array( [ [ 2.0, -2.0, 0.0 ], 
                                                    [ -2.0, 4.0, -2.0 ], 
                                                    [ 0.0, -2.0, 2.0 ] ] )
              self.assertTrue( numpy.allclose( test_stiffness_matrix, gold_stiffness_matrix ) )

       def test_one_quadratic_C0_element( self ):
              problem = { "elastic_modulus": 100,
                     "area": 0.01,
                     "length": 1.0,
                     "traction": { "value": 1e-3, "position": 1.0 },
                     "displacement": { "value": 0.0, "position": 0.0 },
                     "body_force": 0.0 }
              #spline_space = { "domain": [0, problem[ "length" ]], "degree": [ 2 ], "continuity": [ -1, -1 ] }
              #uspline.make_uspline_mesh( spline_space, "test_one_quadratic_C0_element" )
              uspline_bext = bext.readBEXT( "test_one_quadratic_C0_element.json" )
              test_stiffness_matrix = assembleStiffnessMatrix( problem = problem, uspline_bext = uspline_bext )
              gold_stiffness_matrix = numpy.array( [ [  4.0 / 3.0, -2.0 / 3.0, -2.0 / 3.0 ],
                                                 [ -2.0 / 3.0,  4.0 / 3.0, -2.0 / 3.0 ],
                                                 [ -2.0 / 3.0, -2.0 / 3.0,  4.0 / 3.0 ] ] )
              self.assertTrue( numpy.allclose( test_stiffness_matrix, gold_stiffness_matrix ) )

       def test_two_quadratic_C0_element( self ):
              problem = { "elastic_modulus": 100,
                     "area": 0.01,
                     "length": 1.0,
                     "traction": { "value": 1e-3, "position": 1.0 },
                     "displacement": { "value": 0.0, "position": 0.0 },
                     "body_force": 0.0 }
              #spline_space = { "domain": [0, problem[ "length" ]], "degree": [ 2, 2 ], "continuity": [ -1, 0, -1 ] }
              #uspline.make_uspline_mesh( spline_space, "test_two_quadratic_C0_element" )
              uspline_bext = bext.readBEXT( "test_two_quadratic_C0_element.json" )
              test_stiffness_matrix = assembleStiffnessMatrix( problem = problem, uspline_bext = uspline_bext )
              gold_stiffness_matrix = numpy.array( [ [  8.0 / 3.0, -4.0 / 3.0, -4.0 / 3.0,  0.0,        0.0 ],
                                                 [ -4.0 / 3.0,  8.0 / 3.0, -4.0 / 3.0,  0.0,        0.0 ],
                                                 [ -4.0 / 3.0, -4.0 / 3.0, 16.0 / 3.0, -4.0 / 3.0, -4.0 / 3.0 ],
                                                 [  0.0,        0.0,       -4.0 / 3.0,  8.0 / 3.0, -4.0 / 3.0 ],
                                                 [  0.0,        0.0,       -4.0 / 3.0, -4.0 / 3.0,  8.0 / 3.0 ] ] )
              self.assertTrue( numpy.allclose( test_stiffness_matrix, gold_stiffness_matrix ) )

       def test_two_quadratic_C1_element( self ):
              problem = { "elastic_modulus": 100,
                     "area": 0.01,
                     "length": 1.0,
                     "traction": { "value": 1e-3, "position": 1.0 },
                     "displacement": { "value": 0.0, "position": 0.0 },
                     "body_force": 0.0 }
              #spline_space = { "domain": [0, problem[ "length" ]], "degree": [ 2, 2 ], "continuity": [ -1, 1, -1 ] }
              #uspline.make_uspline_mesh( spline_space, "test_two_quadratic_C1_element" )
              uspline_bext = bext.readBEXT( "test_two_quadratic_C1_element.json" )
              test_stiffness_matrix = assembleStiffnessMatrix( problem = problem, uspline_bext = uspline_bext )
              gold_stiffness_matrix = numpy.array( [ [  8.0 / 3.0, -2.0,       -2.0/ 3.0,   0.0 ],
                                                 [ -2.0,        8.0 / 3.0,  0.0,       -2.0 / 3.0 ],
                                                 [ -2.0 / 3.0,  0.0,        8.0 / 3.0, -2.0 ],
                                                 [  0.0,       -2.0 / 3.0, -2.0,        8.0 / 3.0 ] ] )
              self.assertTrue( numpy.allclose( test_stiffness_matrix, gold_stiffness_matrix ) )