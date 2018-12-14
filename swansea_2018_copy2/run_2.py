# 2D shallow water equations
# ================================================================================
from thetis import *
import tidal_forcing
import detectors
import gauges_deterctors
import pickle
import inputs.input_file_paths
#from tools import thetis_support_scripts
import tools.thetis_support_scripts as thetis_support_scripts
from time import gmtime, strftime
import sys

from modules import input_barrages, tools
from modules.lagoon_operation import *

inputdir = 'inputs'
outputdir = inputs.input_file_paths.paraview_output_folder
datadir = 'data'

with timed_stage('reading mesh'):
    mesh2d = Mesh(inputs.input_file_paths.mesh_file)
    #mesh2d = Mesh("inputs/severn_mesh_2.msh")


print_output('Loaded mesh '+ mesh2d.name)
print_output('Exporting to '+ outputdir)

# simulation
#identifier = int(sys.argv[1])
identifier = 0
operation = None

print_output('Simulation identifier : '+str (identifier))


t_end = 30* 24 *3600                               # Simulation duration in sec
t_start = 0 #2530000 #identifier * 30 * 24 * 3600           # Simulation start time relative to tidal_forcing
Dt = 100.                                       # Crank Nicolson timestep
t_export = 1000.0                               # Export time if necessary
wd_alpha = 0.5                                  # Wetting and drying

lat_coriolis = 51                               # Coriolis calculation parameters

P1_2d = FunctionSpace(mesh2d, 'CG', 1)          # Functionspace

bathymetry_2d, h_viscosity, elev_init, uv_init, mu_manning= thetis_support_scripts.initialise_fields(mesh2d, inputdir,
                                                                                                 outputdir, identifier,manning= True)
coriolis_2d = thetis_support_scripts.coriolis(mesh2d, lat_coriolis)
#lagoon_input = input_barrages.input_barrage("inputs/LagoonSpecs.dat")
#lagoon_status = pickle.load(open(inputdir + "/barrage_status_"+str(identifier)+".p", "rb"))
"""
if operation != None: 
   adaptive_control = pickle.load(open(inputdir + "/" + operation +".p","rb"))
   holding_timing = pickle.load(open(inputdir + "/holding_" + operation + ".p","rb"))
else: adaptive_control = None

# Set initial mode time to start of simulation.
Barrages = len(lagoon_status)
for i in range(Barrages):   # initial
    lagoon_status[i]["E"] = 0
"""
with timed_stage('initialisation'):
  # --- create solver ---
  solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry_2d)
  options = solver_obj.options
  options.cfl_2d = 1.0
  options.use_nonlinear_equations = True
  options.simulation_export_time = t_export
  options.simulation_end_time = t_end
  options.coriolis_frequency = coriolis_2d
  options.output_directory = outputdir
  options.check_volume_conservation_2d = True
  options.fields_to_export = ['elev_2d', 'uv_2d']
  options.fields_to_export_hdf5 = ['elev_2d', 'uv_2d']
  options.element_family = "dg-dg"
  options.timestepper_type = 'CrankNicolson'
  options.timestepper_options.implicitness_theta = 1.0
  options.timestepper_options.use_semi_implicit_linearization = True
  options.use_wetting_and_drying = True
  options.wetting_and_drying_alpha = Constant(wd_alpha)
  options.manning_drag_coefficient = mu_manning
  options.horizontal_viscosity = h_viscosity
  options.use_grad_div_viscosity_term = True
  options.use_grad_depth_viscosity_term = False
  options.timestep = Dt  # override dt for CrankNicolson (semi-implicit)
  options.timestepper_options.solver_parameters = {
      'snes_type': 'newtonls',
      'snes_rtol': 1e-3,
      'snes_linesearch_monitor': False,
      'snes_linesearch_type': 'bt',
      'snes_monitor': False,
      'snes_max_it': 20,
      'snes_view': False,
      'ksp_type': 'preonly',
      'snes_converged_reason': True,
      'pc_type': 'lu',
      'pc_factor_mat_solver_package': 'mumps',
  }

# boundary conditions
tidal_elev = Function(bathymetry_2d.function_space())
#lagoon_hydraulic_structures = {"sl_i": Constant(0.), "sl_o": Constant(0.), "tb_i": Constant(0.), "tb_o": Constant(0.)}

river1 = Constant(-10.)
river2 = Constant(-20.)
river3 = Constant(-5.)
Gloucester = Constant(-50.)
# 1 was 11 and 5 was 1
solver_obj.bnd_functions['shallow_water'] = { 1: {'flux': river1}, 2: {'flux': river2},3: {'flux': river3},
                                              4: {'flux': Gloucester},
                                              5: {'elev': tidal_elev},
                                              #6: {'flux': lagoon_hydraulic_structures["tb_i"]},
                                              #7: {'flux': lagoon_hydraulic_structures["sl_i"]},
                                             # 8: {'flux': lagoon_hydraulic_structures["tb_o"]},
                                              #9: {'flux': lagoon_hydraulic_structures["sl_o"]},
                                              }

# Detectors
f_detector = open(datadir + "/detectors_" + str(identifier) + ".dat", "w")

"""
f1 = []
for i in range(Barrages):
    if identifier == 0:
        f1.append(open(datadir+"/Lagoon"+str(i)+"_"+str(operation)+".dat", "w"))
    else:
        f1.append(open(datadir+"/Lagoon"+str(i)+"_"+str(operation)+".dat", "a"))
"""
# Simulation preliminaries
solver_obj.assign_initial_conditions(elev= elev_init, uv=uv_init)
#solver_obj.load_state(2530)
det_xy, det_names = detectors.get_detectors(mesh2d)
cb = DetectorsCallback(solver_obj, det_xy, ['elev_2d', 'uv_2d'], name='detectors', detector_names=det_names)
solver_obj.add_callback(cb, 'timestep')

#gauge loc
det_gauge_xy, det_gauge_names = gauges_deterctors.get_detectors(mesh2d)
cb_gauges = DetectorsCallback(solver_obj, det_gauge_xy, ['elev_2d'], name='gauges_detectors', detector_names=det_gauge_names)
solver_obj.add_callback(cb_gauges, 'timestep')
#cb_lagoon = tools.LagoonCallback(solver_obj, {"inner": dx(3), "outer": dx(4)}, lagoon_input,
#                           thetis_boundaries=lagoon_hydraulic_structures, time=t_start,
##                           number_timesteps=3, name="lagoon_"+str(identifier)+"_"+str(operation),
#                           status=lagoon_status[0], adaptive_control=adaptive_control, export_to_hdf5=True)
#solver_obj.add_callback(cb_lagoon, 'timestep')

uv, elev = solver_obj.timestepper.solution.split()

def intermediate_steps(t,):

    # Exporting to data file - useful for quick sampling etc.
    f_detector.write(str(t))
    for item in elev(inputs.input_file_paths.elevation_detectors):
        f_detector.write(" %s" % item)
    f_detector.write("\n")
    f_detector.flush()

    # Output files for results in text file.
   # i = 0
   # for item in cb_lagoon.output:
   #     f1[i].write(" {:8.3f} ".format(item))
  #  f1[i].write("\n")
 #   f1[i].flush()

    # Exporting to data file - useful for quick sampling etc.
    if inputs.input_file_paths.spatial_harmonics_distribution == True \
            and t % inputs.input_file_paths.elevation_output_interval == 0:

        print_output("Exporting elevation field for harmonic analysis")
        elev_CG = Function(P1_2d, name='elev_CG').project(elev)
        checkpoint_file = checkpointing.DumbCheckpoint(outputdir + '/elev_' + str(t))
        checkpoint_file.store(elev_CG)
        checkpoint_file.close()

    # Export final state that can be picked up later - like load state but including lagoon state
    if t == t_end: thetis_support_scripts.export_final_state(inputdir, identifier, uv, elev)# lagoon = [cb_lagoon.status])

####################### than's code #########################


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


#######################################################################DG
    shear_stress_vector = Function(VectorFunctionSpace(solver.mesh2d, "CG",1)).\
        interpolate(conditional(le(elev+solver.fields["bathymetry_2d"],0), as_vector((0.0,0.0)),dens * C_D * sqrt(dot(uv,uv)) * uv))

    return shear_stress_vector


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


outfile = File(outputdir + '/bed_shear.pvd')
shear = Function(VectorFunctionSpace(mesh2d,'CG',1)) #DG

bed_shear_max = Function(FunctionSpace(mesh2d,'CG',1)) #DG

# for adding up the sher fields over time lapse for averaging
#add = Function(VectorFunctionSpace(mesh2d, 'CG', 1)) #DG
add = Function(FunctionSpace(mesh2d, 'CG', 1)) #DG

#for average velocity
# add_velocity = Function(VectorFunctionSpace(mesh2d, 'DG', 1)) #DG
add_velocity = Function(FunctionSpace(mesh2d, 'CG', 1)) #DG
average_velocity = Function(FunctionSpace(mesh2d, 'CG', 1)) #DG
max_velocity = Function(FunctionSpace(mesh2d, 'CG', 1)) #DG

# define the time-step over which the bed shear stress will be averaged
t_s= t_start
t_e= t_end

def update_forcings(t):

    with timed_stage('lagoon_flux'):
        print_output(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        intermediate_steps(float(t + t_start))

    with timed_stage('update forcings'):
        print_output("Updating tidal field at t={}".format(t_start + t))
        tidal_forcing.set_tidal_field(tidal_elev, t + int(t_start))
        print_output("Done updating tidal field")

        shear.interpolate(compute_bed_shear_term(solver_obj, P1_2d))
        bed_shear_max.interpolate(Max(sqrt(shear[0] ** 2 + shear[1] ** 2), bed_shear_max))
        max_velocity.interpolate(Max(sqrt(uv[0] ** 2 + uv[1] ** 2), max_velocity))

    if t % t_export == 0:
        outfile.write(shear)
        print_output("Exporting bed shear")
        shear_CG = Function(VectorFunctionSpace(mesh2d, 'CG', 1), name='shear').project(shear)  # DG
        checkpoint_file = DumbCheckpoint(outputdir + '/shear_' + str(t))
        checkpoint_file.store(shear_CG)
        checkpoint_file.close()

    # calculates average bed shear stress field over a given time
    if t == t_s:
        add.interpolate(sqrt(shear[0] ** 2 + shear[1] ** 2))
        #add_velocity.assign(uv)
        add_velocity.interpolate(sqrt(uv[0] ** 2 + uv[1] ** 2))

    if (t_s - Dt) < t < (t_e + Dt):
        add.interpolate(add + sqrt(shear[0] ** 2 + shear[1] ** 2))
        #add_velocity.assign(add_velocity + uv)
        add_velocity.interpolate(add_velocity + sqrt(uv[0] ** 2 + uv[1] ** 2))

    if t == t_e:

        #average bed shear
        #average_shear = Function(FunctionSpace(mesh2d, 'CG', 1)).interpolate(sqrt(add[0]**2 +add[1]**2) / ((t_e - t_s) / Dt)) #DG
        average_shear = Function(FunctionSpace(mesh2d, 'CG', 1)).interpolate(add / ((t_e - t_s) / Dt))  # DG
        File(outputdir + '/average_bed_shear_stress.pvd').write(average_shear)

        print_output("Exporting average_shear")
        average_shear_CG = Function(FunctionSpace(mesh2d, 'CG', 1), name='average_shear').project(average_shear) #DG
        checkpoint_file = DumbCheckpoint(outputdir + '/average_shear_' + str(t))
        checkpoint_file.store(average_shear_CG)
        checkpoint_file.close()

        # max bed shear stress
        File(outputdir + '/bed_shear_max.pvd').write(bed_shear_max)
        print_output("Exporting bed_shear_max")
        bed_shear_max_CG = Function(FunctionSpace(mesh2d, 'CG', 1), name='bed_shear_max').project(bed_shear_max)  # DG
        checkpoint_file = checkpointing.DumbCheckpoint(outputdir + '/bed_shear_max_' + str(t))
        checkpoint_file.store(bed_shear_max_CG)
        checkpoint_file.close()

        # average velocity
        # average_velocity.interpolate(sqrt(add_velocity[0]**2 + add_velocity[1]**2) / ((t_end - t_start) / Dt))
        average_velocity.interpolate(add_velocity / ((t_end - t_start) / Dt))
        File(outputdir + '/average_velocity.pvd').write(average_velocity)

        print_output("Exporting average_velocity")
        velocity_CG = Function(FunctionSpace(mesh2d, 'CG', 1), name='average_velocity').project(average_velocity)  # DG
        checkpoint_file = checkpointing.DumbCheckpoint(outputdir + '/average_velocity_' + str(t))
        checkpoint_file.store(velocity_CG)
        checkpoint_file.close()

        #max velocity
        File(outputdir + '/max_velocity.pvd').write(max_velocity)
        print_output("Exporting max_velocity")
        max_velocity_CG = Function(FunctionSpace(mesh2d, 'CG', 1), name='max_velocity').project(max_velocity)  # DG
        checkpoint_file = checkpointing.DumbCheckpoint(outputdir + '/max_velocity_' + str(t))
        checkpoint_file.store(max_velocity_CG)
        checkpoint_file.close()



        # binned_shear = Function(FunctionSpace(mesh2d, 'DG', 1)).interpolate(
        #     conditional(le(sqrt(average_shear[0] ** 2 + average_shear[1] ** 2), 0.25), 1.,
        #                 conditional(le(sqrt(average_shear[0] ** 2 + average_shear[1] ** 2), 0.5), 2.,
        #                             conditional(le(sqrt(average_shear[0] ** 2 + average_shear[1] ** 2),2.0), 3,
        #                                         conditional(le(sqrt(average_shear[0] ** 2 + average_shear[1] ** 2), 8.0), 4,
        #                                                     conditional(le(sqrt(average_shear[0] ** 2 + average_shear[1] ** 2),16.0),5,6))))))
        # File('outputs/binned_bed_shear_stress.pvd').write(binned_shear)
        #
        # out_shapefile_csv('outputs/data_binned_shear.csv', binned_shear, "outputs/binned_shear.shp")


        # binned_max_shear = Function(FunctionSpace(mesh2d, 'DG', 1)).interpolate(
        #     conditional(le(bed_shear_max, 0.25), 1.,
        #                 conditional(le(bed_shear_max, 0.5), 2.,
        #                             conditional(le(bed_shear_max, 2.0), 3,
        #                                         conditional(le(bed_shear_max, 8.0), 4,
        #                                                     conditional(le(bed_shear_max, 16.0), 5, 6))))))
        #
        # File('outputs/binned_max_bed_shear_stress.pvd').write(binned_max_shear)
        #
        # out_shapefile_csv('outputs/data_binned_max_shear.csv', binned_max_shear, "outputs/binned_max_shear.shp")

# Simulation iteration
solver_obj.iterate(update_forcings=update_forcings)
