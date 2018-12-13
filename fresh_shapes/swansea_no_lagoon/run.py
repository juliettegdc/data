# 2D shallow water equations
# ================================================================================
from thetis import *
from time import gmtime, strftime
import tidal_forcing
import detectors

import inputs.input_file_paths
import tools.thetis_support_scripts


inputdir = 'inputs'
outputdir = inputs.input_file_paths.paraview_output_folder
datadir = 'data'

with timed_stage('reading mesh'):
  mesh2d = Mesh(inputs.input_file_paths.mesh_file)

print_output('Loaded mesh '+ mesh2d.name)
print_output('Exporting to '+ outputdir)

# simulation
identifier = 0
print_output('Simulation identifier : '+str (identifier))

runtime = 15 * 24 * 3600
t_start = identifier * runtime
t_end = t_start + runtime
Dt = 100.                                       # Crank Nicolson timestep
t_export = 1000.0                               # Export time if necessary
wd_alpha = 0.5                                  # Wetting and drying

lat_coriolis = 51                               # Coriolis calculation parameters

P1_2d = FunctionSpace(mesh2d, 'CG', 1)          # Functionspace

bathymetry_2d, h_viscosity, elev_init, uv_init, mu_manning= tools.thetis_support_scripts.initialise_fields(mesh2d, inputdir,
                                                                                                 outputdir, identifier,manning= True)
coriolis_2d = tools.thetis_support_scripts.coriolis(mesh2d, lat_coriolis)

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
  options.fields_to_export = []
  options.fields_to_export_hdf5 = []
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

# Boundary conditions
tidal_elev = Function(bathymetry_2d.function_space())
solver_obj.bnd_functions['shallow_water'] = {1: {'elev': tidal_elev},}

# Detectors
f_detector = open(datadir + "/detectors_" + str(identifier) + ".dat", "w")

# Simulation preliminaries
solver_obj.assign_initial_conditions(elev= elev_init, uv=uv_init)
det_xy, det_names = detectors.get_detectors(mesh2d)
cb = DetectorsCallback(solver_obj, det_xy, ['elev_2d', 'uv_2d'], name='detectors', detector_names=det_names)
solver_obj.add_callback(cb, 'timestep')
uv, elev = solver_obj.timestepper.solution.split()

def intermediate_steps(t):

    # Exporting to data file - useful for quick sampling etc.
    f_detector.write(str(t))
    for item in elev(inputs.input_file_paths.elevation_detectors):
        f_detector.write(" %s" % item)
    f_detector.write("\n")
    f_detector.flush()

    # Exporting to data file - useful for quick sampling etc.
    if inputs.input_file_paths.spatial_harmonics_distribution == True \
            and t % inputs.input_file_paths.elevation_output_interval == 0:

        print_output("Exporting elevation field for harmonic analysis")
        elev_CG = Function(P1_2d, name='elev_CG').project(elev)
        checkpoint_file = checkpointing.DumbCheckpoint(outputdir + '/elev_' + str(t))
        checkpoint_file.store(elev_CG)
        checkpoint_file.close()

    # Export final state that can be picked up later - like load state but including lagoon state
    if t == t_end: tools.thetis_support_scripts.export_final_state(inputdir, identifier, uv, elev, lagoon = None )

def update_forcings(t):

    with timed_stage('lagoon_flux'):
      print_output(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
      intermediate_steps(float(t + t_start))

    with timed_stage('update forcings'):
      print_output("Updating tidal field at t={}".format(t_start + t))
      tidal_forcing.set_tidal_field(tidal_elev, t + int(t_start))
      print_output("Done updating tidal field")

# Simulation iteration
solver_obj.iterate(update_forcings=update_forcings)
