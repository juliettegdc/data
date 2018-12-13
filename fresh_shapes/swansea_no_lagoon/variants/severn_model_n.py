# 2D Bristol Channel and Severn Estuary simulation
# ================================================
#

from scipy.interpolate import interp1d
from thetis import *
from time import gmtime, strftime
import math

# Imports for producing bathymetry/viscosity files and tidal forcing
# import distance_from_boundary_severn
import tidal_forcing

# Imports for validation points
import detector_coordinates
mu_manning = 0.023                  # Bed friction
print(str(mu_manning)[3:])
outputdir = 'outputs'
datadir = 'data'

mesh2d = Mesh('severn_mesh_2.msh')
print_output('Loaded mesh '+mesh2d.name)
print_output('Exporting to '+outputdir)

# simulation
identifier = 0
print_output('Simulation identifier : '+str (identifier))

t_end = 30 *  24 * 3600  + identifier * 30 * 24 * 3600                 # Simulation duration
start_time = 0.0 + identifier * 30 * 24 * 3600                         # Simulation start time relative to tidal_forcing
Dt = 100.                           # Crank Nicolson timestep
t_export = 1000.0                   # Export time if necessary
wd_alpha = 0.5                      # Wetting and drying
mu_manning = 0.023                  # Bed friction
theta=51./(180*math.pi)             # Coriolis theta


P1_2d = FunctionSpace(mesh2d, 'CG', 1)          # Initialise Functionspaces
Input_2d = FunctionSpace(mesh2d, 'DG', 1)

# bathymetry
with timed_stage('reading bathymetry'):
  chk = DumbCheckpoint("bathymetry2D", mode=FILE_READ)
  bathymetry_2d = Function(P1_2d, name="bathymetry")
  chk.load(bathymetry_2d)
  File(outputdir+"/bath.pvd").write(bathymetry_2d)
  chk.close()
outfile = File(outputdir+"/depth.pvd")

# viscosity
with timed_stage('initialising viscosity'):
  chk = DumbCheckpoint("viscosity", mode=FILE_READ)
  h_viscosity = Function(P1_2d, name="viscosity")
  chk.load(h_viscosity)
  File(outputdir+"/viscosity.pvd").write(h_viscosity)
  chk.close()

# elevation
with timed_stage('initialising elevation field'):
  chk = DumbCheckpoint("elevation"+str(identifier), mode=FILE_READ)
  elev_init = Function(Input_2d, name="elevation")
  chk.load(elev_init)
  File(outputdir+"/elevation_imported.pvd").write(elev_init)
  chk.close()

# velocity
with timed_stage('initialising velocity field'):
  chk = DumbCheckpoint("velocity"+str(identifier), mode=FILE_READ)
  V = VectorFunctionSpace(mesh2d, 'DG', 1)
  uv_init = Function(V, name="velocity")
  chk.load(uv_init)
  File(outputdir+"/velocity_imported.pvd").write(uv_init)
  chk.close()

with timed_stage('initialisation'):
  solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry_2d)
  options = solver_obj.options
  options.cfl_2d = 1.0
  options.use_nonlinear_equations = True
  options.simulation_export_time = t_export
  options.simulation_end_time = t_end
  options.coriolis_frequency = Constant(sin(theta) * 2 * 2 * math.pi /86400.)
  options.output_directory = outputdir
  options.check_volume_conservation_2d = True
  options.fields_to_export = []
  options.fields_to_export_hdf5 = []
  options.element_family = "dg-dg"
  options.timestepper_type = 'CrankNicolson'
  options.timestepper_options.implicitness_theta = 1.0
  options.timestepper_options.use_semi_implicit_linearization = True
  options.use_wetting_and_drying = True
  options.wetting_and_drying_alpha = Constant(wd_alpha)
  options.manning_drag_coefficient = Constant(mu_manning)
  options.horizontal_viscosity = h_viscosity
  options.use_grad_div_viscosity_term = True
  options.use_grad_depth_viscosity_term = False
  options.timestep = Dt  # override dt for CrankNicolson (semi-implicit)
  options.timestepper_options.solver_parameters = {
      'snes_type': 'newtonls',
      'snes_rtol': 1e-3,
      'snes_linesearch_monitor': True,
      'snes_linesearch_type': 'bt',
      'snes_monitor': True,
      'snes_max_it': 20,
      'snes_view': True,
      'ksp_type': 'preonly',
      'snes_converged_reason': True,
      'pc_type': 'lu',
      'pc_factor_mat_solver_package': 'mumps',
  }
# Hydrodynamic Model boundary conditions
tidal_elev =  Function(bathymetry_2d.function_space())
river1 = Constant(-10.)
river2 = Constant(-20.)
river3 = Constant(-5.)
Gloucester = Constant(-50.)
solver_obj.bnd_functions['shallow_water'] = {2: {'flux': river1}, 3: {'flux': river2}, 4: {'flux': river3}, 6: {'flux': Gloucester},}

for i in (5,):
    solver_obj.bnd_functions['shallow_water'][i] = {'elev': tidal_elev}

# AA10012017 setting monitor points upstream and downstream (including more points in case it balances)
coordinates_down, coordinates_up = [], []

def split(start, end, segments):
    x_delta = (end[0] - start[0]) / float(segments)
    y_delta = (end[1] - start[1]) / float(segments)
    points = []
    for i in range(1, segments):
        points.append([start[0] + i * x_delta, start[1] + i * y_delta])
    return np.array([start] + points + [end])

# Swansea
coordinates_down.append(split((435318,5715559),(435702,5715162),25))
coordinates_up.append(split((435475,5715713),(435862,5715315),25))
# Cardiff
coordinates_down.append(split((490448,5698694),(491349,5698134),25))
coordinates_up.append(split((491845,5699612),(492635,5699052),25))
# Barrage
coordinates_down.append(split((489580,5690318),(494152,5686878),25))
coordinates_up.append(split((489712,5690788),(494475,5687128),25))

f2 = open("detectors.dat", "w")


def lagflux(t_new):
    uv, elev = solver_obj.timestepper.solution.split()
    # Output files for results in text file.
    f2.write(str(t_new))
    for item in elev(detector_coordinates.xy2):
        f2.write(" %s" % item)
    f2.write("\n")

    if t_new == t_end:
        print ("Checkpointing")
        chk = DumbCheckpoint("velocity"+str(identifier+1), mode=FILE_CREATE)
        chk.store(uv, name="velocity")
        File('outputs/velocityout.pvd').write(uv)
        chk.close()
        chk = DumbCheckpoint("elevation"+str(identifier+1), mode=FILE_CREATE)
        chk.store(elev, name="elevation")
        File('outputs/elevationout.pvd').write(elev)
        chk.close()

def update_forcings(t_new):
    with timed_stage('lagoon_flux'):
      print_output(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
      lagflux(float(t_new + start_time))

    with timed_stage('update forcings'):
      print_output("Updating tidal field at t={}".format(t_new + int(start_time)))
      tidal_forcing.set_tidal_field(tidal_elev, t_new + int(start_time), start_time)
      print_output("Done updating tidal field")

parameters['coffee'] = {}

solver_obj.assign_initial_conditions(elev=elev_init, uv=uv_init)
solver_obj.iterate(update_forcings=update_forcings)

f2.close()


