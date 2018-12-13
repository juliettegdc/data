# 2D shallow water equations in a closed channel with (currently) a single turbine
# ================================================================================
from thetis import *
import os.path
import tidal_forcing_ramp
import detectors
import inputs.input_file_paths
import tools.thetis_support_scripts
from time import gmtime, strftime

inputdir = 'inputs'
outputdir = inputs.input_file_paths.paraview_output_folder
datadir = 'data'

with timed_stage('reading mesh'):
  mesh2d = Mesh(inputs.input_file_paths.mesh_file)

print_output('Loaded mesh '+mesh2d.name)
print_output('Exporting to '+outputdir)

# simulation
identifier = -1
print_output('Simulation identifier : '+str (identifier))

ramptime = 2 * 24 * 3600
t_start = - ramptime             # Simulation start time relative to tidal_forcing
t_end = ramptime + t_start       # Simulation duration in sec
Dt = 100.                           # Crank Nicolson timestep
t_export = 20000.0                   # Export time if necessary
wd_alpha = 0.5                      # Wetting and drying
#mu_manning = 0.02                  # Bed friction

lat_coriolis = 51                   # Coriolis calculation parameters

CG_2d = FunctionSpace(mesh2d, 'CG', 1)

bathymetry_2d, h_viscosity, mu_manning = tools.thetis_support_scripts.initialise_fields(mesh2d, inputdir, outputdir, identifier, manning=True)
coriolis_2d = tools.thetis_support_scripts.coriolis(mesh2d, lat_coriolis)

with timed_stage('initialisation'):
  # --- create solver ---
  solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry_2d)
  options = solver_obj.options
  options.cfl_2d = 1.0
  options.use_nonlinear_equations = True
  options.simulation_export_time = t_export
  options.simulation_end_time = ramptime
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
  options.manning_drag_coefficient = Constant(0.02)
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

tidal_elev = Function(bathymetry_2d.function_space())
solver_obj.bnd_functions['shallow_water'] = {1: {'elev': tidal_elev},}

elev_init = Function(CG_2d)
elev_init.assign(0.0)

f2 = open(datadir+"/detectors_"+str(identifier)+".dat", "w")

solver_obj.assign_initial_conditions(uv=as_vector((1e-3, 0.0)), elev = elev_init)
# det_xy, det_names = detectors.get_detectors(mesh2d)
# # cb = DetectorsCallback(solver_obj, det_xy, ['elev_2d', 'uv_2d'], name='detectors', detector_names=det_names)
# solver_obj.add_callback(cb, 'timestep')
uv, elev = solver_obj.timestepper.solution.split()

def intermediate_steps(t):

    # Exporting to data file - useful for quick sampling etc.
    f2.write(str(t))
    for item in elev(inputs.input_file_paths.elevation_detectors):
        f2.write(" %s" % item)
    f2.write("\n")
    f2.flush()

    # Exporting to data file - useful for quick sampling etc.
    if inputs.input_file_paths.spatial_harmonics_distribution == True \
            and t % inputs.input_file_paths.elevation_output_interval== 0:

        print_output("Exporting elevation field for harmonic analysis")
        elev_CG = Function(CG_2d, name='elev_CG').project(elev)
        checkpoint_file = checkpointing.DumbCheckpoint(outputdir + '/elev_' + str(t))
        checkpoint_file.store(elev_CG)
        checkpoint_file.close()

    # Export final state that can be picked up later - like load state but including lagoon state
    if t == t_end: tools.thetis_support_scripts.export_final_state(inputdir, identifier, uv, elev, lagoon = None )

def update_forcings(t):

    intermediate_steps(float(t + t_start))

    print_output("Updating tidal field at t={}".format(t_start + t))
    tidal_forcing_ramp.set_tidal_field(tidal_elev, t + int(t_start), t_start)

solver_obj.iterate(update_forcings=update_forcings)
