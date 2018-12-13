from firedrake import *
from thetis import *
import bathymetry
import tidal_amplitude
import inputs.input_file_paths


# Solution of the Eikonal Function using Firedrake (initial script from Roan, modified by Stephan and Than)
# typical length scale
L = 1e3
inputdir = "inputs"
outputdir = inputs.input_file_paths.paraview_output_folder
mesh = Mesh(inputs.input_file_paths.mesh_file)

# Step 0 - Calculate lowest astronomical tide
V = FunctionSpace(mesh, 'CG', 1)
lat = Function(V)
tidal_amplitude.get_lowest_astronomical_tide(lat)

# Step 1 - Calculate distance for viscosity
print ("Calculate distance for viscosity")

# 5 is the boundary id of the Water Level boundary
# 3 is the inlet of the basin - I'm letting the distance function there start at 1e5 
# so we don't accidentaly apply the higher viscosity there as well
bcs = [DirichletBC(V, 0.0, 1)]   # make sure third digit matches PhysID of open boundaries

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

# Before we solve the Eikonal equation, let's solve a Laplace equation to
# generate an initial guess
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


# epss values set the accuracy (in meters) of the final "distance to boundary" function.
# To make more accurate add in extra iterations, eg, 500., 250., etc. This may result in the solver not converging.

epss = [1000., 500., 400.]

#calculates distance to ocean boundaries of domain in order to calculate velocity and Manning fields below and produces the dist.pvd file
for i, eps in enumerate(epss):
  print ("Solving Eikonal with eps == ", float(eps))
  F = inner(sqrt(inner(grad(u), grad(u))), v) * dx - v * dx + eps*inner(grad(u), grad(v)) * dx
  solve(F == 0, u, bcs, solver_parameters=solver_parameters)

File(outputdir + "/dist.pvd").write(u)

# calculates viscosity fields and produces viscosity.pvd file
chk = DumbCheckpoint(inputdir + "/viscosity", mode=FILE_CREATE)
with timed_stage('initialising viscosity'):
  h_viscosity = Function(V, name="viscosity")
  h_viscosity.interpolate(Max(1., 1000 * (1. - u / 2e4)))
  chk.store(h_viscosity, name="viscosity")
  File(outputdir + '/viscosity.pvd').write(h_viscosity)


# Testing conditional statements
def zbedf(bathy, distance):
    zbed = conditional(ge(bathy,25 *(1.- distance / 10000.)),bathy,25 *(1.- distance / 10000.))
    return zbed
bath = bathymetry.get_bathymetry(inputs.input_file_paths.bathymetry_file, mesh)

bath.assign(bath+lat)
bathymetry.smoothen_bathymetry(bath)
bath.interpolate(Max(zbedf(bath,u), -50.))                          # testing conditional statements


# print ("River boundary bathymetries")
#
# L = 1e3
# # 4 is the boundary id of the inlet
# # 3 is the inlet of the basin - I'm letting the distance function there start at 1e5
# # Set up hydraulic structures
# bc2 = [DirichletBC(V, 0.0, 3),DirichletBC(V, 0.0, 2), DirichletBC(V, 0.0, 4)]          # Double check that
#
# v2 = TestFunction(V)
# u2 = Function(V)
#
# solver_parameters = {
#     'snes_type': 'ksponly',
#     'ksp_rtol': 1e-7,
#     'ksp_type': 'preonly',
#     'pc_type': 'lu',
#     'pc_factor_mat_solver_packages': 'mumps',
#     'ksp_monitor_true_residual': True
#     }
# # # Before we solve the Eikonal equation, let's solve a Laplace equation to
# # # generate an initial guess
# F = L**2*(inner(grad(u2), grad(v2))) * dx - v2 * dx
# solve(F == 0, u2, bc2, solver_parameters=solver_parameters)
# solver_parameters = {
#     'snes_type': 'newtonls',
#     'snes_monitor': True,
#     'ksp_rtol': 1e-7,
#         'ksp_type': 'preonly',
#         'pc_type': 'lu',
#         'pc_factor_mat_solver_packages': 'mumps',
#         }
#
# epss = [1000., 500., 200., 100., 50.]
# for i, eps in enumerate(epss):
#   print ("Solving Eikonal with eps == ", float(eps))
#   F = inner(sqrt(inner(grad(u2), grad(u2))), v2) * dx - v2 * dx + eps*inner(grad(u2), grad(v2)) * dx
#   solve(F == 0, u2, bc2, solver_parameters=solver_parameters)
#
# def zbedf(bathy, distance):
#     zbed = conditional(le(abs(distance),75), 7.0,
#                         conditional(le(abs(distance),125),Max( 7.0 *(1.- distance / 100),bathy),bathy))
#     return zbed
#
# bath.interpolate(zbedf(bath,u2))                          # testing conditional statements

mu_manning = Function(V).interpolate(conditional(le(bath,30),0.023, 0.015))
File(outputdir + '/manning.pvd').write(mu_manning)
chk = DumbCheckpoint(inputdir + "/manning2D", mode=FILE_CREATE)
chk.store(mu_manning, name="manning")

File(outputdir + '/bath.pvd').write(bath)
chk = DumbCheckpoint(inputdir + "/bathymetry2D", mode=FILE_CREATE)
chk.store(bath, name="bathymetry")
