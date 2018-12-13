from thetis import *

mesh2d = Mesh('inputs/mesh.msh')


#I changed the parameters for it to be quicker (but I probably should change the tides parameters too for it to be prettier)
t_end = 3600. # 48 * total duration in seconds
t_export = 200.0    # export interval in seconds
Dt = 100             # timestep


# Omega for amplitude of seaward boundary - assumed to be a sinusoidal function with a given amplitude
amplitude = 4.0
period = 12.42 * 3600
omega = (2 * pi) / period

# Bathymetry and viscosity field
P1_2d = FunctionSpace(mesh2d, 'DG', 1)
bathymetry_2d = Function(P1_2d, name='Bathymetry')
viscosity_2d = Function(P1_2d, name='viscosity')

# Bathymetry combination of bell curve and sloping bathymetry
x, y = SpatialCoordinate(mesh2d)

depth_oce = 20.0
depth_riv = -10.0
mi, sigma = 0, 2000

bathymetry_2d.interpolate(-2 * 1e5 * (-1 / (sigma * sqrt(2 * pi)) * exp(-(y - 3000 - mi)**2 / (2 * sigma**2))) +
                          (depth_riv - depth_oce) * x/16e3)

# Viscosity sponge at the seaward boundary
viscosity_2d.interpolate(conditional(le(x, 2e3), 1e3 * (2e3+1 - x)/2e3, 1))

#manning field
#V = FunctionSpace(mesh2d, 'CG', 1)
#mu_manning = Function(V).interpolate(conditional(le(bathymetry_2d,30),0.022, 0.022))
#File('outputs/manning.pvd').write(mu_manning)
#chk = DumbCheckpoint("inputs/manning2D", mode=FILE_CREATE)
#chk.store(mu_manning, name="manning")#

# create solver
solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry_2d)
options = solver_obj.options
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.check_volume_conservation_2d = True
options.timestepper_type = 'CrankNicolson'
options.fields_to_export = ['uv_2d', 'elev_2d']
#options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
options.timestepper_options.implicitness_theta = 0.5
options.timestepper_options.use_semi_implicit_linearization = True
options.use_wetting_and_drying = True
options.wetting_and_drying_alpha = Constant(0.5)
options.manning_drag_coefficient = Constant(0.022)
options.horizontal_viscosity = viscosity_2d
options.timestep = Dt

# set initial condition for elevation, piecewise linear function
elev_init = Function(P1_2d)
elev_init.assign(0.0)

tidal_elev = Function(bathymetry_2d.function_space())

solver_obj.bnd_functions['shallow_water'] = {4: {'elev': tidal_elev},}
                                             # 5: {'flux': Constant(-1e3)},}


solver_obj.assign_initial_conditions(uv=as_vector((1e-5, 0.0)), elev=elev_init)
#uv, elev = solver_obj.timestepper.solution.split()

##than's code
def wd_bathymetry_displacement(solver,functionspace):
    """
    Returns wetting and drying bathymetry displacement as described in:
    Karna et al.,  2011.
    """
    H = solver.fields["bathymetry_2d"]+solver.fields["elev_2d"]
    disp = Function(functionspace).interpolate(0.5 * (sqrt(H ** 2 + solver.options.wetting_and_drying_alpha ** 2) - H))
    return disp

def compute_total_depth(solver,functionspace):
    """
    Returns effective depth by accounting for the wetting and drying algorithm
    """
    if hasattr(solver.options, 'use_wetting_and_drying') and solver.options.use_wetting_and_drying:
        return Function(functionspace).interpolate(solver.fields["bathymetry_2d"]+solver.fields["elev_2d"]+
                                              wd_bathymetry_displacement(solver, functionspace))

    else:
        return Function(functionspace).assign(solver.fields["bathymetry_2d"]+solver.fields["elev_2d"])


def compute_bed_shear_term(solver,functionspace):

    # Parameters
    g_grav = 9.807
    dens = 1025
    uv, elev = solver.timestepper.solution.split()

    C_D = 0
    if solver_obj.options.quadratic_drag_coefficient is not None:
        C_D = Function(functionspace).assign(solver_obj.options.quadratic_drag_coefficient)

    elif solver_obj.options.manning_drag_coefficient is not None:
        C_D = Function(functionspace).assign(g_grav* solver_obj.options.manning_drag_coefficient**2/
                                             compute_total_depth(solver,functionspace)**(1./3.))



    shear_stress_vector = Function(VectorFunctionSpace(solver.mesh2d, "DG",1)).\
        interpolate(conditional(le(elev+solver.fields["bathymetry_2d"],0), as_vector((0.0,0.0)),dens * C_D * sqrt(dot(uv,uv)) * uv))

    return shear_stress_vector





outfile = File('outputs/bed_shear.pvd')

#P1_2d = FunctionSpace(mesh2d, 'CG', 1)

shear = Function(VectorFunctionSpace(mesh2d,'DG',1))


outfile2 = File('outputs/averge_bed_shear_stress.pvd')
average_shear = Function(VectorFunctionSpace(mesh2d, 'DG', 1))

outfile3 = File('outputs/binned_bed_shear_stress.pvd')
binned_shear = Function(FunctionSpace(mesh2d, 'DG', 1))

#def export_average_bed_shear():
   # t1, t2 = 300, 1000  # define start and end time (t1 nd t2 respectively)
   # nbr_Dt = (t2-t1)/Dt  # define how many shear vector fields will be averaged
   # add = VectorFunctionSpace(mesh2d, 'DG', 1)  # create a vector function space for adding the fields together
   # average_shear = Function(VectorFunctionSpace(mesh2d, 'DG', 1))  # create a function on the vector function space for interpolation

  #  if solver_obj.simulation_time == 1000:  # at start time assign the first shear field
  #      add.interpolate(compute_bed_shear_term(solver_obj,P1_2d))
  #  outfile3.write(add)

   # elif solver_obj.simulation_time > t1:  # over the time range, add up the shear fields in the 'add' function space? should I rather interpolate it?
   #     if solver_obj.simulation_time<=t2:
   #         add.assign(add + compute_bed_shear_term(solver_obj,P1_2d))

  #  elif solver_obj.simulation_time >t2:  # when we re pas the end time, interpolate the average on the function
   #     average_shear.interpolate(add/nbr_Dt)

  #  outfile2.write(average_shear)  # export the field of average bed shear stress

add = Function(VectorFunctionSpace(mesh2d, 'DG', 1))

add_velocity = Function(VectorFunctionSpace(mesh2d, 'DG', 1))
average_velocity = Function(VectorFunctionSpace(mesh2d, 'DG', 1))
# define the timestep over which the bed shear stress will be averaged
t_start=1000
t_end=3000
uv, p = solver_obj.timestepper.solution.split()



def update_forcings(t_new,):
    tidal_elev.assign(Constant(tanh((t_new)/(4*3600.)) * amplitude * sin(omega * t_new)))



    shear.interpolate(compute_bed_shear_term(solver_obj, P1_2d))
    if t_new % t_export == 0:
        outfile.write(shear)

    #calculat average bed shear stress field over a defined timestep
    if t_new == t_start:
        add.assign(shear)
        add_velocity.assign(uv)

    if (t_start - Dt) < t_new < (t_end + Dt):
        add.assign(add + shear)
        add_velocity.assign(add_velocity + uv)

    if t_new == t_end:
        #average_shear = Function(VectorFunctionSpace(mesh2d, 'DG', 1)).interpolate(add/((t_end-t_start)/Dt))
        average_shear.interpolate(add / ((t_end - t_start) / Dt))
        average_velocity.interpolate(add_velocity / ((t_end - t_start) / Dt))
        File('outputs/average_velocity.pvd').write(average_velocity)
        outfile2.write(average_shear)


        binned_shear.interpolate(conditional(le(sqrt(average_shear[0]**2 + average_shear[1]**2 ),0.015), 1,
                                             conditional(le(sqrt(average_shear[0]**2 + average_shear[1]**2 ),0.04),2,3)))
        outfile3.write(binned_shear)


# Solve the system
solver_obj.iterate(update_forcings=update_forcings,)

###################################### try to output vlue at coordintes ################################################

#txt fle ith x, y and numbers
"""
import numpy as np
x, y= SpatialCoordinate(mesh2d)
x_vector, y_vector = interpolate(x, Function(P1_2d)).dat.data, interpolate(y, Function(P1_2d)).dat.data

#print (uv.dat.data[5,])

for i in range(int(t_end/t_export)+1):# going through the exported data
      #print('Reading h5 files. Time ',i,i*t_export)
    solver_obj.load_state(i)
    print(uv.dat.data[5,0], 'and',uv.dat.data[5,1] )


import csv
with open('data.csv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(binned_shear.dat.data,x_vector,y_vector))

locations = list(zip(x_vector, y_vector))
numbers = list(binned_shear.dat.data)

import pyproj
import shapely.geometry
import fiona
import fiona.crs

UTM_ZONE30 = pyproj.Proj(
        proj='utm',
        zone=30,
        datum='WGS84',
        units='m',
        errcheck=True)
LL_WGS84 = pyproj.Proj(proj='latlong', datum='WGS84', errcheck=True)

schema = {'geometry': 'Point', 'properties': {'numbers': 'str'}}
crs = fiona.crs.from_string(UTM_ZONE30.srs)
with fiona.collection("benned_shear.shp", "w", "ESRI Shapefile", schema, crs=crs) as output:
    for xy, numbers in zip(locations, numbers):
        point = shapely.geometry.Point(xy[0], xy[1])
        output.write({'properties': {'numbers': numbers}, 'geometry': shapely.geometry.mapping(point)})
"""