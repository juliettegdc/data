from thetis import *

mesh2d = Mesh('inputs/severn_mesh_2.msh')


#I changed the parameters for it to be quicker (but I probably should change the tides parameters too for it to be prettier)
t_end = 30000. # 48 * total duration in seconds
t_export = 1000.0    # export interval in seconds
Dt = 200             # timestep


# Omega for amplitude of seaward boundary - assumed to be a sinusoidal function with a given amplitude
amplitude = 4.0
period = 12.42 * 3600
omega = (2 * pi) / period

# Bathymetry and viscosity field
P1_2d = FunctionSpace(mesh2d, 'DG', 1)
bathymetry_2d = Function(P1_2d, name='Bathymetry')
viscosity_2d = Function(P1_2d, name='viscosity')

#simple bathymetry
#depth = 20.0
#bathymetry_2d.assign(depth)

# Bathymetry combination of bell curve and sloping bathymetry
x, y = SpatialCoordinate(mesh2d)

depth_oce = 60.0
depth_riv = -10.0
mi, sigma = 0, 2000

bathymetry_2d.interpolate(-2 * 1e5 * (-1 / (sigma * sqrt(2 * pi)) * exp(-(y - 3000 - mi)**2 / (2 * sigma**2))) +
                          (depth_riv - depth_oce) * x/16e3)


# Viscosity sponge at the seaward boundary
"""
viscosity_2d.interpolate(conditional(le(x, 2e3), 1e3 * (2e3+1 - x)/2e3, 1))
"""
##########################################################################################################################
L = 1e3

#possibly useless in our case?
# Step 0 - Calculate lowest astronomical tide
"""
V = FunctionSpace(mesh2d, 'DG', 1)
lat = Function(V)
tidal_amplitude.get_lowest_astronomical_tide(lat)
"""
V = FunctionSpace(mesh2d, 'CG', 1) #keep it CG for som reason it doesn't want it DG

# Step 1 - Calculate distance for viscosity
print ("Calculate distance for viscosity")

# 5 is the boundary id of the Water Level boundary
# 3 is the inlet of the basin - I'm letting the distance function there start at 1e5
# so we don't accidentaly apply the higher viscosity there as well
bcs = [DirichletBC(V, 0.0, 1)]

v = TestFunction(V)
u = Function(V)

solver_parameters = {
    'snes_type': 'ksponly',
    'ksp_rtol': 1e-4,
    'ksp_type': 'preonly',
    'pc_type': 'lu',
    'pc_factor_mat_solver_packages': 'mumps',
    'ksp_monitor_true_residual': True
    }

# # Before we solve the Eikonal equation, let's solve a Laplace equation to
# # generate an initial guess
F = L**2*(inner(grad(u), grad(v))) * dx - v * dx
solve(F == 0, u, bcs, solver_parameters=solver_parameters)

solver_parameters = {
    'snes_type': 'newtonls',
    'snes_monitor': True,
    'ksp_rtol': 1e-4,
        'ksp_type': 'preonly',
        'pc_type': 'lu',
        'pc_factor_mat_solver_packages': 'mumps',
        }
# from the overleaf doc: epss values set the accuracy (in meters)
#  of the final "distance to boundary function.
# to make more accurate add in extra iterations, eg 500., 250., etc.
#this may result in the solver not converging.
epss = [20000., 1000., 500., 400.]
for i, eps in enumerate(epss):
  print ("Solving Eikonal with eps == ", float(eps))
  F = inner(sqrt(inner(grad(u), grad(u))), v) * dx - v * dx + eps*inner(grad(u), grad(v)) * dx
  solve(F == 0, u, bcs, solver_parameters=solver_parameters)

File("outputs/dist.pvd").write(u)


chk = DumbCheckpoint("inputs/viscosity", mode=FILE_CREATE)
with timed_stage('initialising viscosity'):
    h_viscosity = Function(V, name="viscosity")
    h_viscosity.interpolate(Max(1., 1000 * (1. - u / 2e4)))
    chk.store(h_viscosity, name="viscosity")
    File('outputs//viscosity.pvd').write(h_viscosity)
#####################################################################################################################








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

solver_obj.bnd_functions['shallow_water'] = {1: {'elev': tidal_elev},}
                                             # 5: {'flux': Constant(-1e3)},}


solver_obj.assign_initial_conditions(uv=as_vector((1e-3, 0.0)), elev=elev_init)
#uv, elev = solver_obj.timestepper.solution.split()

##than's code
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

outfile2 = File('outputs/averge_bed_shear_stress.pvd')

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

# define the timestep over which the bed shear stress will be averaged
t_start=1000
t_end=3000

def update_forcings(t_new,):
    tidal_elev.assign(Constant(tanh((t_new)/(4*3600.)) * amplitude * sin(omega * t_new)))

    shear.interpolate(compute_bed_shear_term(solver_obj, P1_2d))
    if t_new % t_export == 0:
        outfile.write(shear)

    #calculat average bed shear stress field over a defined timestep
    if t_new == t_start:
        add.assign(shear)

    if (t_start - Dt) < t_new < (t_end + Dt):
        add.assign(add + shear)

    if t_new == t_end:
        average_shear = Function(VectorFunctionSpace(mesh2d, 'DG', 1)).interpolate(add/((t_end-t_start)/Dt))
        outfile2.write(average_shear)


# Solve the system
solver_obj.iterate(update_forcings=update_forcings,)


#def update_forcings(t_new,):
    #tidal_elev.assign(Constant(tanh((t_new)/(4*3600.)) * amplitude * sin(omega * t_new)))
    #if t_new % t_export == 0:
       # shear.interpolate(compute_bed_shear_term(solver_obj, P1_2d))
        #outfile.write(shear)


