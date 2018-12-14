# 2D shallow water equations in a closed channel with (currently) a single turbine
# ================================================================================
from thetis import *
import tidal_forcing
import inputs.input_file_paths
import tools.thetis_support_scripts
from modules import input_barrages
from modules.tools import LagoonCallback
import pickle

inputdir = 'inputs'
outputdir = inputs.input_file_paths.paraview_output_folder
datadir = 'data'

with timed_stage('reading mesh'):
    mesh2d = Mesh(inputs.input_file_paths.mesh_file)

print_output('Loaded mesh '+mesh2d.name)
print_output('Exporting to '+outputdir)

# simulation
#identifier = -1
#print_output('Simulation identifier : '+str (identifier))

t_end = 30 * 24 * 3600                          # Simulation duration in sec
t_start =  0#identifier * 30 * 24 * 3600           # Simulation start time relative to tidal_forcing
Dt = 100.                                       # Crank Nicolson timestep
t_export = 1000.0                               # Export time if necessary
wd_alpha = 0.5                                  # Wetting and drying

lat_coriolis = 51                               # Coriolis calculation parameters

CG_2d = FunctionSpace(mesh2d, 'CG', 1)
P1_2d = FunctionSpace(mesh2d, 'CG', 1)

# simulation
# identifier = int(sys.argv[1])
# ctrl_ID = int(sys.argv[2])
#
#
#
# if ctrl_ID == 0 : operation = None
# elif ctrl_ID == 1 : operation = 'cardiff_ebb'
# elif ctrl_ID == 2 : operation = 'cardiff_two-way'

identifier = 0
operation = None

bathymetry_2d, h_viscosity,  elev_init, uv_init, mu_manning = tools.thetis_support_scripts.initialise_fields(mesh2d, inputdir, outputdir, identifier, manning=True)
coriolis_2d = tools.thetis_support_scripts.coriolis(mesh2d, lat_coriolis)

# lagoon_input = input_barrages.input_barrage("inputs/LagoonSpecs.dat")
# lagoon_swansea_input = input_barrages.input_barrage("/data/swansea_2018_copy/inputs/LagoonSpecs_original.dat")
# lagoon_cardiff_input = input_barrages.input_barrage("inputs/LagoonSpecs.dat")

lagoon_swansea_input = input_barrages.input_barrage("../swansea_2018_copy/inputs/LagoonSpecs.dat")
lagoon_cardiff_input = input_barrages.input_barrage("inputs/LagoonSpecs.dat") # inputs/LagoonSpecs.dat

lagoon_status = pickle.load(open(inputdir + "/barrage_status_"+str(identifier)+".p", "rb")) #should add one?

if operation != None:
   adaptive_control = pickle.load(open(inputdir + "/" + operation +".p","rb"))
   holding_timing = pickle.load(open(inputdir + "/holding_" + operation + ".p","rb"))
else: adaptive_control = None

Barrages = len(lagoon_status)
for i in range(Barrages):   # initial
    lagoon_status[i]["E"] = 0


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
  options.fields_to_export = ['uv_2d','elev_2d']
  options.fields_to_export_hdf5 = []
  options.element_family = "dg-dg"
  options.timestepper_type = 'CrankNicolson'
  options.timestepper_options.implicitness_theta = 1.0
  options.timestepper_options.use_semi_implicit_linearization = True
  options.use_wetting_and_drying = True
  options.wetting_and_drying_alpha = Constant(wd_alpha)
  options.manning_drag_coefficient = Constant(0.022)
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
lagoon_hydraulic_structures_swansea = {"sl_i": Constant(0.), "sl_o": Constant(0.), "tb_i": Constant(0.), "tb_o": Constant(0.)}
lagoon_hydraulic_structures_cardiff = {"sl_i": Constant(0.), "sl_o": Constant(0.), "tb_i": Constant(0.), "tb_o": Constant(0.)}
river1 = Constant(-10.)
river2 = Constant(-20.)
river3 = Constant(-5.)
Gloucester = Constant(-50.)

solver_obj.bnd_functions['shallow_water'] = { 1: {'flux': river1}, 2: {'flux': river2}, 3: {'flux': river3},
                                              4: {'flux': Gloucester},
                                              5: {'elev': tidal_elev},

                                              6: {'flux': lagoon_hydraulic_structures_swansea["tb_i"]},
                                              7: {'flux': lagoon_hydraulic_structures_swansea["sl_i"]},
                                              8: {'flux': lagoon_hydraulic_structures_swansea["tb_o"]},
                                              9: {'flux': lagoon_hydraulic_structures_swansea["sl_o"]},

                                              10: {'flux': lagoon_hydraulic_structures_cardiff["tb_i"]},
                                              11: {'flux': lagoon_hydraulic_structures_cardiff["sl_i"]},
                                              12: {'flux': lagoon_hydraulic_structures_cardiff["tb_o"]},
                                              13: {'flux': lagoon_hydraulic_structures_cardiff["sl_o"]},
                                              }

elev_init = Function(CG_2d)
elev_init.assign(0.0)

f2 = open(datadir+"/detectors_"+str(identifier)+".dat", "w")
f1 = open(datadir+"/Lagoon_swansea_"+str(identifier)+".dat", "w")
f3 = open(datadir+"/Lagoon_cardiff_"+str(identifier)+".dat", "w")

solver_obj.assign_initial_conditions(uv=as_vector((1e-3, 0.0)), elev = elev_init)

cb_lagoon1 = LagoonCallback(solver_obj, {"inner": dx(5), "outer": dx(6)}, lagoon_swansea_input,
                           thetis_boundaries=lagoon_hydraulic_structures_swansea, time=t_start,
                           number_timesteps=5, name="lagoon_swansea_"+str(identifier), adaptive_control=None, export_to_hdf5=True)

cb_lagoon2 = LagoonCallback(solver_obj, {"inner": dx(3), "outer": dx(4)}, lagoon_cardiff_input,
                           thetis_boundaries=lagoon_hydraulic_structures_cardiff, time=t_start,
                           number_timesteps=5, name="lagoon_cardiff_"+str(identifier), adaptive_control=None, export_to_hdf5=True)

solver_obj.add_callback(cb_lagoon1, 'timestep')
solver_obj.add_callback(cb_lagoon2, 'timestep')

uv, elev = solver_obj.timestepper.solution.split()

def intermediate_steps(t):

    # Exporting to data file - useful for quick sampling etc.
    f2.write(str(t))
    for item in elev(inputs.input_file_paths.elevation_detectors):
        f2.write(" %s" % item)
    f2.write("\n")
    f2.flush()

    for item in cb_lagoon1.output:
        f1.write(" {:8.3f} ".format(item))
    f1.write("\n")
    f1.flush()

    for item in cb_lagoon2.output:
        f3.write(" {:8.3f} ".format(item))
    f3.write("\n")
    f3.flush()

    # Exporting to data file - useful for quick sampling etc.
    if inputs.input_file_paths.spatial_harmonics_distribution == True \
            and t % inputs.input_file_paths.elevation_output_interval== 0:

        print_output("Exporting elevation field for harmonic analysis")
        elev_CG = Function(CG_2d, name='elev_CG').project(elev)
        checkpoint_file = checkpointing.DumbCheckpoint(outputdir + '/elev_' + str(t))
        checkpoint_file.store(elev_CG)
        checkpoint_file.close()

    # Export final state that can be picked up later - like load state but including lagoon state
    if t == t_end: tools.thetis_support_scripts.export_final_state(inputdir, identifier, uv, elev, lagoon = [cb_lagoon1.status, cb_lagoon2.status])


# def update_forcings(t):
#
#     intermediate_steps(float(t + t_start))
#
#     print_output("Updating tidal field at t={}".format(t_start + t))
#     tidal_forcing_ramp.set_tidal_field(tidal_elev, t + int(t_start), t_start)

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


# DGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
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
# add = Function(VectorFunctionSpace(mesh2d, 'CG', 1)) #DG
add = Function(FunctionSpace(mesh2d, 'CG', 1)) #DG

#for average velocity
# add_velocity = Function(VectorFunctionSpace(mesh2d, 'DG', 1)) #DG
add_velocity = Function(FunctionSpace(mesh2d, 'CG', 1)) #DG
average_velocity = Function(FunctionSpace(mesh2d, 'CG', 1)) #DG
max_velocity = Function(FunctionSpace(mesh2d, 'CG', 1))

# define the time-step over which the bed shear stress will be averaged
t_s= t_start
t_e= t_end

def update_forcings(t):

    """
    with timed_stage('lagoon_flux'):
      print_output(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
      intermediate_steps(float(t + t_start))
    """
    with timed_stage('update forcings'):
        print_output("Updating tidal field at t={}".format(t_start + t))
        tidal_forcing.set_tidal_field(tidal_elev, t + int(t_start))
        print_output("Done updating tidal field")

    shear.interpolate(compute_bed_shear_term(solver_obj, P1_2d))
    bed_shear_max.interpolate(Max(sqrt(shear[0] ** 2 + shear[1] ** 2), bed_shear_max))
    max_velocity.interpolate(Max(sqrt(uv[0]**2 + uv[1]**2), max_velocity))

    if t % t_export == 0:
        outfile.write(shear)
        shear_CG = Function(VectorFunctionSpace(mesh2d, 'CG', 1), name='shear').project(shear)  # DG
        checkpoint_file = DumbCheckpoint(outputdir + '/shear_' + str(t))
        checkpoint_file.store(shear_CG)
        checkpoint_file.close()

    # calculates average bed shear stress field over a given time
    if t == t_s:
        #add.assign(abs(shear))
        add.interpolate(sqrt(shear[0] ** 2 + shear[1] ** 2))
        # add_velocity.assign(abs(uv))
        add_velocity.interpolate(sqrt(uv[0]**2 + uv[1]**2))

    if (t_s - Dt) < t < (t_e + Dt):
        # add.assign(add + abs(shear))
        add.interpolate(add + sqrt(shear[0] ** 2 + shear[1] ** 2))
        # add_velocity.assign(add_velocity + abs(uv))
        add_velocity.interpolate(add_velocity + sqrt(uv[0]**2 + uv[1]**2))


    if t == t_e:

        #average bed shear
        # average_shear = Function(FunctionSpace(mesh2d, 'CG', 1)).interpolate(sqrt(add[0]**2 + add[1]**2) / ((t_e - t_s) / Dt)) #DG
        average_shear = Function(FunctionSpace(mesh2d, 'CG', 1)).interpolate(add / ((t_e - t_s) / Dt))  # DG
        File(outputdir + '/average_bed_shear_stress.pvd').write(average_shear)

        print_output("Exporting average_shear")
        average_shear_CG = Function(FunctionSpace(mesh2d, 'CG', 1), name='average_shear').project(average_shear) #DG
        checkpoint_file = DumbCheckpoint(outputdir + '/average_shear_test_' + str(t))
        checkpoint_file.store(average_shear_CG)
        checkpoint_file.close()

        # max bed shear stress
        File(outputdir + '/bed_shear_max.pvd').write(bed_shear_max)
        print_output("Exporting bed_shear_max")
        bed_shear_max_CG = Function(FunctionSpace(mesh2d, 'CG', 1), name='bed_shear_max').project(bed_shear_max)  # DG
        checkpoint_file = checkpointing.DumbCheckpoint(outputdir + '/bed_shear_max_test_' + str(t))
        checkpoint_file.store(bed_shear_max_CG)
        checkpoint_file.close()

        #averge velocity
        #average_velocity.interpolate(sqrt(add_velocity[0]**2 + add_velocity[1]**2) / ((t_end - t_start) / Dt))
        average_velocity.interpolate(add_velocity / ((t_end - t_start) / Dt))
        File(outputdir + '/average_velocity.pvd').write(average_velocity)

        print_output("Exporting average_velocity")
        velocity_CG = Function(FunctionSpace(mesh2d, 'CG', 1), name='average_velocity').project(average_velocity) #DG
        checkpoint_file = checkpointing.DumbCheckpoint(outputdir + '/average_velocity_test_' + str(t))
        checkpoint_file.store(velocity_CG)
        checkpoint_file.close()

        # max velocity
        File(outputdir + '/max_velocity.pvd').write(max_velocity)
        print_output("Exporting max_velocity")
        max_velocity_CG = Function(FunctionSpace(mesh2d, 'CG', 1), name='max_velocity').project(max_velocity)  # DG
        checkpoint_file = checkpointing.DumbCheckpoint(outputdir + '/max_velocity_' + str(t))
        checkpoint_file.store(max_velocity_CG)
        checkpoint_file.close()

solver_obj.iterate(update_forcings=update_forcings)
