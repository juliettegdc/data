
###################################################
# First, grab the mesh.
m = V.ufl_domain()

# Now make the VectorFunctionSpace corresponding to V.
W = VectorFunctionSpace(m, V.ufl_element())

# Next, interpolate the coordinates onto the nodes of W.
X = interpolate(m.coordinates, W)

# Make an output function.
f = Function(V)

# Use the external data function to interpolate the values of f.
f.dat.data[:] = mydata(X.dat.data_ro)

####################################################
# Expression:
f = interpolate(Expression("sin(x[0]*pi)"), V)

# UFL equivalent:
x = SpatialCoordinate(V.mesh())
f = interpolate(sin(x[0] * math.pi), V)

# Expression with a Constant parameter:
f = interpolate(Expression('sin(x[0]*t)', t=t), V)

# UFL equivalent:
x = SpatialCoordinate(V.mesh())
f = interpolate(sin(x[0] * t), V)

##################################################

# create a vectorfunction space for uv
W = VectorFunctionSpace(mesh2d, "DG", 1)
uv = Function(W, name='vel_2d')
print(uv.dat.data.shape)

# create arrays
u_data_set = np.empty((int(t_end/t_export)+1,uv.dat.data.shape[0]))
v_data_set = np.empty((int(t_end/t_export)+1,uv.dat.data.shape[0]))
elev_data_set = np.empty((int(t_end/t_export)+1,uv.dat.data.shape[0]))


count = 0
for i in range(0,int(t_end/t_export)+1): # going through the exported data
      #print('Reading h5 files. Time ',i,i*t_export)
    solver_obj.load_state(i)
    u_data_set[count, :] = solver_obj.fields.uv_2d.dat.data[:, 0]  # extract u componnt at time-step i
    v_data_set[count, :] = solver_obj.fields.uv_2d.dat.data[:, 1]  #extract v component at time-step i
    elev_data_set[count, :] = solver_obj.fields.elev_2d.dat.data[:]  # extract elevation at timestep i


    #print(u_data_set[1], v_data_set[1], u[1])
    count += 1

# extract x and y to use in a function
DG_2d = FunctionSpace(mesh2d, 'DG', 1)
x = SpatialCoordinate(mesh2d)
x_vector, y_vector = interpolate(x[0], Function(DG_2d)).dat.data, interpolate(x[1], Function(DG_2d)).dat.data


#try to see velo field

#u, v = interpolate(solver_obj.fields.uv_2d[0], Function(DG_2d)).dat.data, interpolate(solver_obj.fields.uv_2d[1], Function(DG_2d)).dat.data

#print(u[:], u_data_set[1], v_data_set[1])



##############################################################################################



  options.manning_drag_coefficient = manning
  #options.quadratic_drag_coefficient = Constant(0.003)
  options.element_family = "dg-dg"
  options.timestepper_type = 'CrankNicolson'
  options.timestepper_options.use_semi_implicit_linearization = False
  options.timestepper_options.solver_parameters = {
        'snes_type': 'newtonls',
        'snes_rtol': 1e-2,
        'ksp_rtol': 1e-3,
        'snes_monitor': True,
        'ksp_converged_reason': True,
        'ksp_type': 'gmres',
        'pc_type': 'fieldsplit',
    }
#      'snes_type': 'newtonls',
#      'snes_rtol': 1e-3,
#    #     'snes_linesearch_monitor': True,
#      'snes_linesearch_type': 'bt',
#      'snes_monitor': True,
#      'snes_max_it': 20,
#    #     'snes_view': True,
#      'ksp_type': 'preonly',
#      'snes_converged_reason': False,
#      'ksp_converged_reason': False,
#      'pc_type': 'lu',
#      'pc_factor_mat_solver_package': 'mumps',
#    }

  print_output("#### Sovler options ####")
  print_output(solverObj.options)
  print_output("#### Timestepper options ####")
  print_output(solverObj.options.timestepper_options)


  # boundary conditions
  tidal_elev = Function(bathymetry2d.function_space())
  solverObj.bnd_functions['shallow_water'] = {
          2000: {'elev': tidal_elev},
          1000: {'un': 0.0},
  }

  solverObj.assign_initial_conditions(uv=Constant((1.0e-10,0.0)))

  if solverObj.comm.rank == 0:
    det_file = open(os.path.join(output_dir,'detectors.txt'), 'wb')

  uv, p = solverObj.timestepper.solution.split()
  valid_xy = detectors.get_valid_detectors(p)
  print_output("Number of valid detectors {}".format(len(valid_xy)))

  def update_forcings(t):
    with timed_stage('writing detectors'):
      uvxy = array(uv(valid_xy)).T
      pxy = array(p(valid_xy))[newaxis,:]
      if solverObj.comm.rank == 0:
        savetxt(det_file, uvxy)
        savetxt(det_file, pxy)
    with timed_stage('update forcings'):
      print_output("Updating tidal field at t={}".format(t))
      tidal_forcing.set_tidal_field(tidal_elev, t)
      print_output("Done updating tidal field")

  speed_max_pvd = File(os.path.join(output_dir, 'Speed_max.pvd'))
  speed_max = Function(V)
  speed = sqrt(uv[0]**2 + uv[1]**2)

  def export_func():
    if solverObj.simulation_time > 86400.0 :
      speed_max.interpolate(Max(speed, speed_max))
    speed_max_pvd.write(speed_max)


  #solverObj.load_state(245, outputdir='outputs_2013_cd25',
  #   iteration=245*36., t=245*3600.)

solverObj.iterate(update_forcings=update_forcings, export_func=export_func)
speed_max_pvd.close()
if solverObj.comm.rank == 0:
    det_file.close()


###########################################################################################################

uv, elev = solver_obj.timestepper.solution.split()

def wd_bathymetry_displacement(solver,functionspace):
    """
    Returns wetting and drying bathymetry displacement as described in:
    Karna et al.,  2011.
    """

    H = solver.fields["bathymetry_2d"]+solver.fields["elev_2d"]
    disp = Function(functionspace).assign(0.5 * (sqrt(H ** 2 + solver.options.wetting_and_drying_alpha ** 2) - H))
    return disp


def compute_total_depth(solver,functionspace):
    """
    Returns effective depth by accounting for the wetting and drying algorithm
    """

    if hasattr(solver.options, 'use_wetting_and_drying') and solver.options.use_wetting_and_drying:
        return Function(functionspace).assign(solver.fields["bathymetry_2d"]+solver.fields["elev_2d"]+
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

    shear_stress_vector = Function(VectorFunctionSpace(solver.mesh2d, "DG",1)).interpolate(dens * C_D * sqrt(dot(uv,uv)) * uv)
    return shear_stress_vector

outfile = File('outputs/bed_shear.pvd')
shear = Function(VectorFunctionSpace(mesh2d,'DG',1))

def update_forcings(t_new,):
    tidal_elev.assign(Constant(tanh((t_new)/(4*3600.)) * amplitude * sin(omega * t_new)))

    if t_new % t_export == 0:
        shear.interpolate(compute_bed_shear_term(solver_obj,P1_2d))
        outfile.write(shear)


#########################################################################################################
##################### attempt 1 at averaging the bed shear stress over a time step ######################

outfile2 = File('output/average_bed_shear.pvd')

def export_average_bed_shear():
    t1, t2 = 300, 1000  # define start and end time (t1 nd t2 respectively)
    nbr_Dt = (t2-t1)/Dt  # define how many shear vector fields will be averaged
    add = VectorFunctionSpace(mesh2d, 'DG', 1)  # create a vector function space for adding the fields together
    average_shear = Function(VectorFunctionSpace(mesh2d, 'DG', 1))  # create a function on the vector function space for interpolation

    if solver_obj.simulation_time == t1:  # at start time assign the first shear field
        add.assign(shear)

    elif solver_obj.simulation_time > t1:  # over the time range, add up the shear fields in the 'add' function space? should I rather interpolate it?
        if solver_obj.simulation_time<=t2:
            add.assign(add + shear)

    elif solver_obj.simulation_time >t2:  # when we re pas the end time, interpolate the average on the function
        average_shear.interpolate(add/nbr_Dt)

    outfile2.write(average_shear)  # export the field of average bed shear stress



def update_forcings(t_new,):
    tidal_elev.assign(Constant(tanh((t_new)/(4*3600.)) * amplitude * sin(omega * t_new)))

    shear.interpolate(compute_bed_shear_term(solver_obj, P1_2d))
    if t_new % t_export == 0:
        outfile.write(shear)



# Solve the system
solver_obj.iterate(update_forcings=update_forcings,export_func=export_average_bed_shear)



######################### error coming up :Writing different set of functions#####################################


##################################################################################################################
##################### attempt 2 at exporting average bed shear stress over a time step ###########################

# Solve the system
solver_obj.iterate(update_forcings=update_forcings,)

outfile2 = File('output/average_bed_shear.pvd')
#average_shear = Function(VectorFunctionSpace(mesh2d,'DG',1))

t1, t2 = 300, 1000  # define start and end time (t1 nd t2 respectively)
nbr_Dt = (t2-t1)/Dt  # define how many shear vector fields will be averaged
add = VectorFunctionSpace(mesh2d, 'DG', 1)  # create a vector function space for adding the fields together
average_shear = Function(VectorFunctionSpace(mesh2d, 'DG', 1))  # create a function on the vector function space for interpolation

if solver_obj.simulation_time == t1:  # at start time assign the first shear field
    add.assign(compute_bed_shear_term(solver_obj, P1_2d))

elif solver_obj.simulation_time > t1:  # over the time range, add up the shear fields in the 'add' function space? should I rather interpolate it?
    while solver_obj.simulation_time <= t2:
            add.assign(add + compute_bed_shear_term(solver_obj, P1_2d))

elif solver_obj.simulation_time > t2:  # when we re pas the end time, interpolate the average on the function
    average_shear.interpolate(add/nbr_Dt)

    outfile2.write(average_shear)

######################################################################################################################

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

P1_2d = FunctionSpace(mesh2d, 'CG', 1)

shear = Function(VectorFunctionSpace(mesh2d,'DG',1))


#########################################################################################################
#atemps at extracting long and lat from mesh / convrting x, y to long lqt coordinqtes

import utm

x_0, y_0, utm_zone, zone_letter = utm.from_latlon(lat, 0)


#########################################################################################################

from thetis import *
import pyproj
import numpy as np
import inputs.input_file_paths

tidegauge_file = inputs.input_file_paths.tidegauge_file

UTM_ZONE30 = pyproj.Proj(
        proj='utm',
        zone=30,
        datum='WGS84',
        units='m',
        errcheck=True)
LL_WGS84 = pyproj.Proj(proj='latlong', datum='WGS84', errcheck=True)


def get_detectors(mesh2d):
    gauge_names = np.loadtxt(tidegauge_file, skiprows=1, usecols=(0,), dtype=str, delimiter=',')
    gauge_xy = np.loadtxt(tidegauge_file, skiprows=1, usecols=(3,4), delimiter=',')
    ind = np.argsort(gauge_names)
    gauge_names = list(gauge_names[ind])
    gauge_xy = list(gauge_xy[ind])

    #make names unique
    unique_names = []
    last_name = ''; ctr = 0
    for name in gauge_names:
        if name==last_name:
            unique_names.append(name + '_' + str(ctr))
            ctr += 1
        else:
            unique_names.append(name)
            ctr = 0
        last_name = name

    # add custom detectors examples:

    #gauge_latlon.append([58.659477538216585, -3.141400563744444])
    #unique_names.append('Exeter2013')
    #xy = [pyproj.transform(LL_WGS84, crs, lon, lat) for lat, lon in gauge_latlon]
    return select_and_move_detectors(mesh2d, gauge_xy, unique_names, maximum_distance=10e3)

if __name__ == "__main__":
    mesh2d = Mesh(inputs.input_file_paths.mesh_file)

    locations, names = get_detectors(mesh2d)
    if mesh2d.comm.rank == 0: # only processor 0
        print_output("Found detectors: {}".format(names))
        # write out shape-file
        import shapely.geometry
        import fiona
        import fiona.crs

        schema = {'geometry': 'Point', 'properties': {'name': 'str'}}
        crs = fiona.crs.from_string(UTM_ZONE30.srs)
        with fiona.collection("data/detectors.shp", "w", "ESRI Shapefile", schema, crs=crs) as output:
            for xy, name in zip(locations, names):
                point = shapely.geometry.Point(xy[0], xy[1])
                output.write({'properties': {'name': name}, 'geometry': shapely.geometry.mapping(point)})


mesh2d = Mesh(inputs.input_file_paths.mesh_file)
locations, numbers = xy.dat.data, binned_shear.dat.data

import shapely.geometry
import fiona
import fiona.crs

schema = {'geometry': 'Point', 'properties': {'numbers': 'str'}}
crs = fiona.crs.from_string(UTM_ZONE30.srs)
with fiona.collection("benned_shear.shp", "w", "ESRI Shapefile", schema, crs=crs) as output:
    for xy, numbers in zip(locations, numbers):
        point = shapely.geometry.Point(xy[0], xy[1])
        output.write({'properties': {'name': numbers}, 'geometry': shapely.geometry.mapping(point)})





x,y =SpatialCoordinate(mesh)

x.dat.data

################################################################################################################
from qgis.core import *
import processing

layer1 = processing.getObject('MyPointsLayer')
layer2 = processing.getObject('MyPolygonsLayer')

index = QgsSpatialIndex() # Spatial index
for ft in layer1.getFeatures():
    index.insertFeature(ft)

selection = [] # This list stores the features which contains at least one point
for feat in layer2.getFeatures():
    inGeom = feat.geometry()
    idsList = index.intersects(inGeom.boundingBox())
    if idsList:
        selection.append(feat)

# Select all the polygon features which contains at least one point
layer2.setSelectedFeatures([k.id() for k in selection])


###################################################################################################################

def out_shapefile_csv(csv_file_name, field_to_export, shape_file_name):

    x, y = SpatialCoordinate(mesh2d)
    x_vector, y_vector = interpolate(x, Function(P1_2d)).dat.data, interpolate(y, Function(P1_2d)).dat.data

    import csv
    with open(csv_file_name, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(zip(x_vector, y_vector, field_to_export.dat.data))

    locations = list(zip(x_vector, y_vector))
    numbers = list(field_to_export.dat.data)

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
    with fiona.collection(shape_file_name, "w", "ESRI Shapefile", schema, crs=crs) as output:
        for xy, numbers in zip(locations, numbers):
            point = shapely.geometry.Point(xy[0], xy[1])
            output.write({'properties': {'numbers': numbers}, 'geometry': shapely.geometry.mapping(point)})


################################################################################################################
import numpy as np
from matplotlib import pyplot as plt
import datetime
import uptide
from firedrake import *
from thetis import *
from netCDF4 import Dataset as NetCDFFile
from scipy.interpolate import interp2d
import utm
import sys
import pyproj
import uptide.tidal_netcdf
import h5py
import matplotlib.pyplot as plt

output_dir = 'comparison'
create_directory(output_dir)
file_location1 = 'modified/tidal_model/output/hdf5/'  # location of the Elevation2d output files
file_location2 = '1995/output/hdf5/'
constituent_order = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', 'M4', 'MS4', 'MN4']  # LEAVE THIS
constituents = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', 'M4', 'MS4',
                'MN4']  # CHOSEN CONSTITUENTS MUST BE IN THE ORDER ABOVE

# t_end =  2724300 #2377 hdf5 files
t_end = 450000  # 2377 hdf5 files
start_file = 192
t_export = 900

t_n = int(t_end / t_export + 1)
thetis_times = t_export * np.arange(t_n) + t_export

mesh1 = Mesh('modified/tidal_model/runfiles/mesh.msh')
mesh2 = Mesh('1995/runfiles/mesh.msh')

xyvector1 = mesh1.coordinates.dat.data  # UTM51, nodes
xyvector2 = mesh2.coordinates.dat.data

P1DG_1 = FunctionSpace(mesh1, "DG", 1)
P1DG_2 = FunctionSpace(mesh2, "DG", 1)

x = xyvector1[:, 0]
y = xyvector1[:, 1]

max_x = 344135.7783181509
min_x = 322639.8772734722
max_y = 474965.0681569836
min_y = 458630.5856613767

n = 400
gridxy = np.mgrid[min_x:max_x:400j, min_y:max_y:400j].T
gridxy = np.reshape(gridxy, (n * n, 2))

'''
n = 400 # sqrt(no. points)
gridxy = np.mgrid[min(x):max(x):400j, min(y):max(y):400j].T # Ny x Nx x 2
gridxy = np.reshape(gridxy, (n*n, 2)) # N x 2
'''

gridx, gridy = gridxy.T

elev1 = Function(P1DG_1, name='elev_2d')
elev1_data_set = np.empty((t_n, gridxy.shape[0]))

elev2 = Function(P1DG_2, name='elev_2d')
elev2_data_set = np.empty((t_n, gridxy.shape[0]))

# for i,xy in enumerate(gridxy):
#    for i in range(len(gridxy)):
#        print (xy)

# Reading in elevation data
for i in range(start_file, int(t_end / t_export) + 1):
    print('Reading h5 files. Time  ', i, i * t_export)
    checkpoint_file1 = DumbCheckpoint(file_location1 + '/Elevation2d_{:05}'.format(i), mode=FILE_READ)
    checkpoint_file1.load(elev1)
    checkpoint_file1.close()
    f1 = elev1.at(gridxy, dont_raise=True)  # interpolate elev1 at gridxy points
    f1 = [np.array([np.nan]) if x == None else x for x in f1]
    f1 = np.array(f1)
    elev1_data_set[i, :] = f1[:]

    checkpoint_file2 = DumbCheckpoint(file_location2 + '/Elevation2d_{:05}'.format(i), mode=FILE_READ)
    checkpoint_file2.load(elev2)
    checkpoint_file2.close()
    f2 = elev2.at(gridxy, dont_raise=True)  # interpolate elev2 at gridxy points
    f2 = [np.array([np.nan]) if x == None else x for x in f2]
    f2 = np.array(f2)
    elev2_data_set[i, :] = f2[:]

    # empty_array[:] = f1[:] # work around an error
    # f1 = empty_array[:]

    # f1 = np.array(f1)
    # f1 = f1.reshape(n,n)

    # plt.imshow(f1, origin='lower', cmap='coolwarm')
    # plt.colorbar()
    # plt.show()

detector_amp1 = []
detector_phase1 = []

detector_amp2 = []
detector_phase2 = []

detector_ampd = []
detector_phased = []

print(elev1_data_set.shape)
print('\nHarmonic analysis.\n')
for i in range(elev1_data_set.shape[1]):  # loop over each point
    # Mesh1
    thetis_elev1 = elev1_data_set[:, i] * 10  # take elevation of point at all timesteps

    tide = uptide.Tides(constituents)
    tide.set_initial_time(datetime.datetime(2018, 1, 1, 0, 0))

    thetis_elev1 = thetis_elev1 - thetis_elev1.mean()
    thetis_amplitudes1, thetis_phases1 = uptide.analysis.harmonic_analysis(tide, thetis_elev1[int(start_file):],
                                                                           thetis_times[int(start_file):])

    detector_amp1.append(thetis_amplitudes1)
    detector_phase1.append(thetis_phases1)

    # Mesh2
    thetis_elev2 = elev2_data_set[:, i] * 10

    thetis_elev2 = thetis_elev2 - thetis_elev2.mean()
    thetis_amplitudes2, thetis_phases2 = uptide.analysis.harmonic_analysis(tide, thetis_elev2[int(start_file):],
                                                                           thetis_times[int(start_file):])

    detector_amp2.append(thetis_amplitudes2)
    detector_phase2.append(thetis_phases2)

    # Difference
    amp_diff = thetis_amplitudes2 - thetis_amplitudes1
    phase_diff = thetis_phases2 - thetis_phases1
    detector_ampd.append(amp_diff)
    detector_phased.append(phase_diff)

M2_amp_1 = np.array(detector_amp1)[:, constituents.index('M2')]
M2_amp_2 = np.array(detector_amp2)[:, constituents.index('M2')]
M2_amp_d = np.array(detector_ampd)[:, constituents.index('M2')]  # shape = N = Nx . Ny

M2_amp_1 = M2_amp_1.reshape(n, n)
M2_amp_2 = M2_amp_2.reshape(n, n)
M2_amp_d = M2_amp_d.reshape(n, n)  # shape = Nx x Ny

# M2_amp_sihwa = M2_amp_d[:int(n/2),:int(n/2)]

if COMM_WORLD.rank == 0:
    figd = plt.imshow(M2_amp_d, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
                      cmap='coolwarm')
    plt.colorbar()
    plt.savefig('M2_amp_diff')
    plt.show()

    fig1 = plt.imshow(M2_amp_1, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
                      cmap='coolwarm')
    plt.colorbar()
    plt.savefig('M2_amp_1')
    plt.show()

    fig2 = plt.imshow(M2_amp_1, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
                      cmap='coolwarm')
    plt.colorbar()
    plt.savefig('M2_amp_2')
    plt.show()

#############################################################################################################
# def treatlines():
#     boundaries = qmesh.vector.Shapes()
#     boundaries.fromFile('lagoon_Swansea.shp')
#     loopShapes = qmesh.vector.identifyLoops(boundaries,
#                                             isGlobal=False, defaultPhysID=1000,
#                                             fixOpenLoops=True)
#     loopShapes.writeFile('Extended_BCM_lagoon_Swansea_1.shp')



def mesh(name):
    '''Todo: add docstring '''
    # Reading in the shapefile describing the domain boundaries, and creating a gmsh file.

    boundaries = qmesh.vector.Shapes()
    boundaries.fromFile('severn_cardiff.shp')

    loopShapes = qmesh.vector.identifyLoops(boundaries,
                                            isGlobal=False, defaultPhysID=1000, fixOpenLoops=True)
   polygonShapes = qmesh.vector.identifyPolygons(loopShapes, smallestNotMeshedArea=300,
                                                  meshedAreaPhysID=1)


    inner = qmesh.vector.Shapes()
    inner.fromFile('inner_swan.shp')
    inner_loops = qmesh.vector.identifyLoops(inner,
                                             fixOpenLoops=True, extraPointsPerVertex=10)

    inner_polygon = qmesh.vector.identifyPolygons(inner_loops, meshedAreaPhysID = 3)



    inner2 = qmesh.vector.Shapes()
    inner2.fromFile('inner_2.shp')
    inner2_loops = qmesh.vector.identifyLoops(inner2,
                                             fixOpenLoops=True, extraPointsPerVertex=10)
    inner2_polygon = qmesh.vector.identifyPolygons(inner2_loops, meshedAreaPhysID = 3)


    outer = qmesh.vector.Shapes()
    outer.fromFile('outer_swan.shp')
    outer_loops = qmesh.vector.identifyLoops(outer,
                                             fixOpenLoops=True, extraPointsPerVertex=10)
    outer_polygon = qmesh.vector.identifyPolygons(outer_loops, meshedAreaPhysID = 4)


    outer2 = qmesh.vector.Shapes()
    outer2.fromFile('outer_2.shp')
    outer2_loops = qmesh.vector.identifyLoops(outer2,
                                             fixOpenLoops=True, extraPointsPerVertex=10)
    outer2_polygon = qmesh.vector.identifyPolygons(outer2_loops, meshedAreaPhysID = 4)


    # Create raster for mesh gradation towards full-resolution shorelines.


    GSHHS_fine_boundaries = qmesh.vector.Shapes()
    GSHHS_fine_boundaries.fromFile('grad_isl.shp')
    grad_0 = qmesh.raster.meshMetricTools.gradationToShapes()
    grad_0.setShapes(GSHHS_fine_boundaries)
    grad_0.setRasterBounds(-8.0, -2.0, 50.0, 53.0)
    grad_0.setRasterResolution(1300, 1300)
    grad_0.setGradationParameters(100.0, 2000.0, 0.5, 0.001)
    grad_0.calculateLinearGradation()
    grad_0.writeNetCDF('grad_isl.nc')


    GSHHS_fine_boundaries = qmesh.vector.Shapes()
    GSHHS_fine_boundaries.fromFile('grad_100.shp')
    grad_1 = qmesh.raster.meshMetricTools.gradationToShapes()
    grad_1.setShapes(GSHHS_fine_boundaries)
    grad_1.setRasterBounds(-8.0, -2.0, 50.0, 53.0)
    grad_1.setRasterResolution(1300, 1300)
    grad_1.setGradationParameters(150.0, 8000.0, 1)
    grad_1.calculateLinearGradation()
    grad_1.writeNetCDF('grad_100.nc')



    GSHHS_fine_boundaries = qmesh.vector.Shapes()
    GSHHS_fine_boundaries.fromFile('grad_200.shp')
    grad_2 = qmesh.raster.meshMetricTools.gradationToShapes()
    grad_2.setShapes(GSHHS_fine_boundaries)
    grad_2.setRasterBounds(-8.0, -2.0, 50.0, 53.0)
    grad_2.setRasterResolution(1000, 1000)
    grad_2.setGradationParameters(300.0, 8000.0, 1)
    grad_2.calculateLinearGradation()
    grad_2.writeNetCDF('grad_200.nc')



    GSHHS_coarser_boundaries = qmesh.vector.Shapes()
    GSHHS_coarser_boundaries.fromFile('grad_500.shp')
    grad_3 = qmesh.raster.meshMetricTools.gradationToShapes()
    grad_3.setShapes(GSHHS_coarser_boundaries)
    grad_3.setRasterBounds(-8.0, -2.0, 50.0, 53.0)
    grad_3.setRasterResolution(1300, 1300)
    grad_3.setGradationParameters(1000.0, 8000.0, 1.0)
    grad_3.calculateLinearGradation()
    grad_3.writeNetCDF('grad_500.nc')



    GSHHS_coarser_boundaries = qmesh.vector.Shapes()
    GSHHS_coarser_boundaries.fromFile('cardiff_ref_structures.shp')
    grad_4 = qmesh.raster.meshMetricTools.gradationToShapes()
    grad_4.setShapes(GSHHS_coarser_boundaries)
    grad_4.setRasterBounds(-8.0, -2.0, 50.0, 53.0)
    grad_4.setRasterResolution(1300, 1300)
    grad_4.setGradationParameters(60., 10000.0, 1.0, 0.001)
    grad_4.calculateLinearGradation()



    GSHHS_coarser_boundaries = qmesh.vector.Shapes()
    GSHHS_coarser_boundaries.fromFile('cardiff_ref_impoundment.shp')
    grad_5 = qmesh.raster.meshMetricTools.gradationToShapes()
    grad_5.setShapes(GSHHS_coarser_boundaries)
    grad_5.setRasterBounds(-8.0, -2.0, 50.0, 53.0)
    grad_5.setRasterResolution(1300, 1300)
    grad_5.setGradationParameters(150., 10000.0, 1.0, 0.001)
    grad_5.calculateLinearGradation()



    # domainLines, domainPolygons = qmesh.vector.insertRegions(loopShapes, polygonShapes, inner_loops, inner_polygon)
    # domainLines, domainPolygons = qmesh.vector.insertRegions(domainLines, domainPolygons, outer_loops, outer_polygon)
    domainLines, domainPolygons = qmesh.vector.insertRegions(loopShapes, polygonShapes, inner2_loops, inner2_polygon)
    domainLines, domainPolygons = qmesh.vector.insertRegions(domainLines, domainPolygons, outer2_loops, outer2_polygon)


    # Calculate overall mesh-metric raster
    meshMetricRaster = qmesh.raster.meshMetricTools.minimumRaster([grad_1, grad_2, grad_3, grad_0, grad_4, grad_5]) #grad_0
    meshMetricRaster.writeNetCDF('meshMetric.nc')
    # Create domain object and write gmsh files.
    domain = qmesh.mesh.Domain()
    domain.setGeometry(domainLines, domainPolygons)
    # domain.setGeometry(loopShapes, polygonShapes)
    domain.setMeshMetricField(meshMetricRaster)
    domain.setTargetCoordRefSystem('EPSG:32630', fldFillValue=1000.0)
    # Meshing


    domain.gmsh(geoFilename= name + '.geo', \
                fldFilename= name + '.fld', \
                mshFilename= name + '.msh', \
                )

#

def convertMesh(name):
    mesh = qmesh.mesh.Mesh()
    mesh.readGmsh( name + '.msh', 'EPSG:32630')
    mesh.writeShapefile(name)





if __name__ == '__main__':
    import qmesh
    import os

    print os.getcwd()
    os.chdir("cardiff_2018/")
    name = "cardiff_lagoon_2018"
    # Initialising qgis API
    qmesh.initialise()

#    treatlines()
    mesh(name)
    convertMesh(name)