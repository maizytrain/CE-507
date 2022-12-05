import sys
# add cubit libraries to your path
if __name__ == "CubitPythonInterpreter_2":
    # We are running within the Coreform Cubit application, cubit python module is already available
    pass
else:
    sys.path.append("/Applications/Coreform-Cubit-2022.11.app/Contents/MacOS")
import cubit
cubit.init([])

cubit.cmd( "reset" )
# Create geometry
cubit.cmd( "create curve location 0 0 0 location 1 0 0" )
# Create initial mesh topology
cubit.cmd( "curve 1 interval 2" )
cubit.cmd( "mesh curve 1" )
# Create U-spline fit to the geometry
cubit.cmd( "set uspline curve 1 degree 2 continuity 1" )
cubit.cmd( "build uspline curve 1 as 1" )
cubit.cmd( "fit uspline 1" )
# Export U-spline Bezier extraction data in JSON file format
cubit.cmd( "export uspline 1 json ’quadratic_bspline’" )