import unittest
import readBEXTJSON as read


class Test_assembleGramMatrix( unittest.TestCase ):
    def test_two_element_linear_bspline( self ):
        target_fun = lambda x: x**0
        spline_space = { "domain": [0, 2], "degree": [ 1, 1 ], "continuity": [ -1, 0, -1 ] }
        uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        uspline_bext = bext.readBEXT( "temp_uspline.json" )
        test_gram_matrix = splineApproxGalerkin.assembleGramMatrix( uspline_bext = uspline_bext )
        gold_gram_matrix = numpy.array( [ [ 1.0/3.0, 1.0/6.0, 0.0 ],
                                          [ 1.0/6.0, 2.0/3.0, 1.0/6.0 ],
                                          [ 0.0, 1.0/6.0, 1.0/3.0 ] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )

    def test_two_element_quadratic_bspline( self ):
        target_fun = lambda x: x**0
        spline_space = { "domain": [0, 2], "degree": [ 2, 2 ], "continuity": [ -1, 1, -1 ] }
        uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        uspline_bext = bext.readBEXT( "temp_uspline.json" )
        test_gram_matrix = splineApproxGalerkin.assembleGramMatrix( uspline_bext = uspline_bext )
        gold_gram_matrix = numpy.array( [ [ 1.0/5.0, 7.0/60.0, 1.0/60.0, 0.0 ],
                                          [ 7.0/60.0, 1.0/3.0, 1.0/5.0, 1.0/60.0],
                                          [ 1.0/60.0, 1.0/5.0, 1.0/3.0, 7.0/60.0 ],
                                          [ 0.0, 1.0/60.0, 7.0/60.0, 1.0/5.0] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )

    def test_two_element_cubic_bspline( self ):
        spline_space = { "domain": [0, 2], "degree": [ 3, 3 ], "continuity": [ -1, 2, -1 ] }
        uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        uspline_bext = bext.readBEXT( "temp_uspline.json" )
        test_gram_matrix = splineApproxGalerkin.assembleGramMatrix( uspline_bext = uspline_bext )
        gold_gram_matrix = numpy.array( [ [ 1.0/7.0, 7.0/80.0, 1.0/56.0, 1.0/560.0, 0.0 ],
                                          [ 7.0/80.0, 31.0/140.0, 39.0/280.0, 1.0/20.0, 1.0/560.0 ],
                                          [ 1.0/56.0, 39.0/280.0, 13.0/70.0, 39.0/280.0, 1.0/56.0 ],
                                          [ 1.0/560.0, 1.0/20.0, 39.0/280.0, 31.0/140.0, 7.0/80.0 ],
                                          [ 0.0, 1.0/560.0, 1.0/56.0, 7.0/80.0, 1.0/7.0 ] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )