# 2D shallow water equations in a closed channel with (currently) a single turbine
# ================================================================================
from thetis import *
import tidal_forcing_ramp
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
identifier = -1
print_output('Simulation identifier : '+str (identifier))

t_end = 2000 #30 * 24 * 3600                          # Simulation duration in sec
t_start =  0#identifier * 30 * 24 * 3600           # Simulation start time relative to tidal_forcing
Dt = 100.                                       # Crank Nicolson timestep
t_export = 1000.0                               # Export time if necessary
wd_alpha = 0.5                                  # Wetting and drying

lat_coriolis = 51                               # Coriolis calculation parameters

CG_2d = FunctionSpace(mesh2d, 'CG', 1)


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

lagoon_input = input_barrages.input_barrage("inputs/LagoonSpecs.dat")
lagoon_status = pickle.load(open(inputdir + "/barrage_status_"+str(identifier)+".p", "rb"))

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
lagoon_hydraulic_structures = {"sl_i": Constant(0.), "sl_o": Constant(0.), "tb_i": Constant(0.), "tb_o": Constant(0.)}

river1 = Constant(-10.)
river2 = Constant(-20.)
river3 = Constant(-5.)
Gloucester = Constant(-50.)

solver_obj.bnd_functions['shallow_water'] = { 100: {'flux': river1}, 200: {'flux': river2}, 300: {'flux': river3},
                                              400: {'flux': Gloucester},
                                              1: {'elev': tidal_elev},
                                              10: {'flux': lagoon_hydraulic_structures["tb_i"]},
                                              11: {'flux': lagoon_hydraulic_structures["sl_i"]},
                                              12: {'flux': lagoon_hydraulic_structures["tb_o"]},
                                              13: {'flux': lagoon_hydraulic_structures["sl_o"]},
                                              }


elev_init = Function(CG_2d)
elev_init.assign(0.0)

f2 = open(datadir+"/detectors_"+str(identifier)+".dat", "w")
f1 = open(datadir+"/Lagoon_"+str(identifier)+".dat", "w")

solver_obj.assign_initial_conditions(uv=as_vector((1e-3, 0.0)), elev = elev_init)

cb_lagoon = LagoonCallback(solver_obj, {"inner": dx(3), "outer": dx(4)}, lagoon_input,
                           thetis_boundaries=lagoon_hydraulic_structures, time=t_start,
                           number_timesteps=5, name="lagoon_"+str(identifier), adaptive_control=None, export_to_hdf5=True)

solver_obj.add_callback(cb_lagoon, 'timestep')

uv, elev = solver_obj.timestepper.solution.split()

def intermediate_steps(t):

    # Exporting to data file - useful for quick sampling etc.
    f2.write(str(t))
    for item in elev(inputs.input_file_paths.elevation_detectors):
        f2.write(" %s" % item)
    f2.write("\n")
    f2.flush()

    for item in cb_lagoon.output:
        f1.write(" {:8.3f} ".format(item))
    f1.write("\n")
    f1.flush()

    # Exporting to data file - useful for quick sampling etc.
    if inputs.input_file_paths.spatial_harmonics_distribution == True \
            and t % inputs.input_file_paths.elevation_output_interval== 0:

        print_output("Exporting elevation field for harmonic analysis")
        elev_CG = Function(CG_2d, name='elev_CG').project(elev)
        checkpoint_file = checkpointing.DumbCheckpoint(outputdir + '/elev_' + str(t))
        checkpoint_file.store(elev_CG)
        checkpoint_file.close()

    # Export final state that can be picked up later - like load state but including lagoon state
    if t == t_end: tools.thetis_support_scripts.export_final_state(inputdir, identifier, uv, elev, lagoon = [cb_lagoon.status])


def update_forcings(t):

    intermediate_steps(float(t + t_start))

    print_output("Updating tidal field at t={}".format(t_start + t))
    tidal_forcing_ramp.set_tidal_field(tidal_elev, t + int(t_start), t_start)

solver_obj.iterate(update_forcings=update_forcings)
