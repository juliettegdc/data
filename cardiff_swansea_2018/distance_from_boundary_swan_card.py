from thetis import *
import inputs.input_file_paths
import bathymetry
import tidal_amplitude

import os

# Solution of the Eikonal Function using Firedrake (initial script from Roan, modified by Stephan and Angeloudis)
# typical length scale
L = 1e3
inputdir = "inputs"
outputdir = inputs.input_file_paths.paraview_output_folder
mesh = Mesh(inputs.input_file_paths.mesh_file)


# Step 0 - Calculate lowest astronomical tide
V = FunctionSpace(mesh, 'CG', 1)
lat = Function(V)
tidal_amplitude.get_lowest_astronomical_tide(lat)
# File('lat.pvd').write(lat)

# Step 1 - Calculate distance for viscosity
print ("Calculate distance for viscosity")


# 5 is the boundary id of the Water Level boundary
# 3 is the inlet of the basin - I'm letting the distance function there start at 1e5
# so we don't accidentaly apply the higher viscosity there as well
bcs = [DirichletBC(V, 0.0, 5), DirichletBC(V, 1e5, 7), DirichletBC(V, 1e5, 10)]

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
    'snes_type': 'ksponly',
    'ksp_rtol': 1e-4,
    'ksp_type': 'preonly',
    'pc_type': 'lu',
    'pc_factor_mat_solver_packages': 'mumps',
    'ksp_monitor_true_residual': True #his line is different than other codes
    }

epss = [20000., 1000., 500., 400.]
for i, eps in enumerate(epss):
    print("Solving Eikonal with eps == ", float(eps))
    F = inner(sqrt(inner(grad(u), grad(u))), v) * dx - v * dx + eps * inner(grad(u), grad(v)) * dx
    solve(F == 0, u, bcs, solver_parameters=solver_parameters)

File(outputdir + "/dist.pvd").write(u)

chk = DumbCheckpoint(inputdir + "/viscosity", mode=FILE_CREATE)
with timed_stage('initialising viscosity'):
    h_viscosity = Function(V, name="viscosity")
    h_viscosity.interpolate(Max(1, 1000 * (1. - u / 2e4)))
    chk.store(h_viscosity, name="viscosity")
    File(outputdir + '/viscosity.pvd').write(h_viscosity)


# Testing conditional statements
def zbedf(bathy, distance):
    zbed = conditional(ge(bathy,50 *(1.- distance / 10000.)),bathy,50 *(1.- distance / 10000.))
    return zbed
bath = bathymetry.get_bathymetry(inputs.input_file_paths.bathymetry_file, mesh)

# # bath.assign(bath + lat)
# bath.interpolate(Max(zbedf(bath,u), -50.))                          # testing conditional statements
# # bath.interpolate(Max(Max(bath, 50*(1.-u/10000.)), -50.))          # this works fine too!
# bathymetry.smoothen_bathymetry(bath)
# bathymetry.smoothen_bathymetry(bath)

bath.assign(bath + lat)
bathymetry.smoothen_bathymetry(bath)
bath.interpolate(Max(zbedf(bath,u), -50.))                          # testing conditional statements


#File('outputs/bath.pvd').write(bath)


print("Swansea Lagoon boundary bathymetries")

L = 1e3
# 4 is the boundary id of the inlet
# 3 is the inlet of the basin - I'm letting the distance function there start at 1e5
# Set up hydraulic structures
bc2 = [DirichletBC(V, 0.0, 6), DirichletBC(V, 0.0, 7),          # Double check that
       DirichletBC(V, 0.0, 9), DirichletBC(V, 0.0, 8),
       DirichletBC(V, 1e5, 10) ]

v2 = TestFunction(V)
u2 = Function(V)

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

F = L**2*(inner(grad(u2), grad(v2))) * dx - v2 * dx
solve(F == 0, u2, bc2, solver_parameters=solver_parameters)
solver_parameters = {
    'snes_type': 'newtonls',
    'snes_monitor': True,
    'ksp_rtol': 1e-4,
        'ksp_type': 'preonly',
        'pc_type': 'lu',
        'pc_factor_mat_solver_packages': 'mumps',
        }

epss =[50000., 25000., 10000., 5000., 1000., 500., 250., 100., 75.]
for i, eps in enumerate(epss):
    print("Solving Eikonal with eps == ", float(eps))
    F = inner(sqrt(inner(grad(u2), grad(u2))), v2) * dx - v2 * dx + eps * inner(grad(u2), grad(v2)) * dx
    solve(F == 0, u2, bc2, solver_parameters=solver_parameters)

# File("outputs/dist.pvd").write(u2)

# Testing conditional statements   / bath has been introduced previously
def zbedf(bathy, distance):
    zbed = conditional(le(distance,75),25.0 - distance/75 * 5,
                       conditional(le(distance,200),20,
                                   conditional(le(distance,600),Max(20-distance/400*15,bathy),bathy)))
    return zbed
bath.interpolate(zbedf(bath,u2))                          # testing conditional statements



print("Cardiff Lagoon boundary bathymetries")

L = 1e3
# 4 is the boundary id of the inlet
# 3 is the inlet of the basin - I'm letting the distance function there start at 1e5
# Set up hydraulic structures
bc3 = [DirichletBC(V, 0.0, 10), DirichletBC(V, 0.0, 11),  # Double check that
       DirichletBC(V, 0.0, 12), DirichletBC(V, 0.0, 13),
       DirichletBC(V, 1e5, 7)  ]

v3 = TestFunction(V)
u3 = Function(V)

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
F = L**2*(inner(grad(u3), grad(v3))) * dx - v3 * dx
solve(F == 0, u3, bc3, solver_parameters=solver_parameters)
solver_parameters = {
    'snes_type': 'newtonls',
    'snes_monitor': True,
    'ksp_rtol': 1e-4,
        'ksp_type': 'preonly',
        'pc_type': 'lu',
        'pc_factor_mat_solver_packages': 'mumps',
        }

epss = [50000., 25000., 10000., 5000., 1000., 500., 250., 100., 75.]

for i, eps in enumerate(epss):
    print("Solving Eikonal with eps == ", float(eps))
    F = inner(sqrt(inner(grad(u3), grad(u3))), v3) * dx - v3 * dx + eps*inner(grad(u3), grad(v3)) * dx
    solve(F == 0, u3, bc3, solver_parameters=solver_parameters)


# File("outputs/dist.pvd").write(u2)

# Testing conditional statements   / bath has been introduced previously

def zbedf(bathy, distance):
    zbed = conditional(le(distance,75),30.0 - distance/75 * 10,
                       conditional(le(distance,200),20,
                                   conditional(le(distance,600),Max(20-distance/400*15,bathy),bathy)))
    return zbed
bath.interpolate(zbedf(bath,u3))                          # testing conditional statements





print("River boundary bathymetries")

L = 1e3

# 4 is the boundary id of the inlet
# 3 is the inlet of the basin - I'm letting the distance function there start at 1e5
# Set up hydraulic structures
# bc2 = [DirichletBC(V, 1e5, 5), DirichletBC(V, 0.0, 3), DirichletBC(V, 0.0, 6),
#        DirichletBC(V, 0.0, 4), DirichletBC(V, 0.0, 7)]

bc2 = [DirichletBC(V, 0.0, 1), DirichletBC(V, 0.0, 2),          # Double check that
        DirichletBC(V, 0.0, 3), DirichletBC(V, 0.0, 4),
        DirichletBC(V, 1e5, 7), DirichletBC(V, 1e5, 10) ]

v2 = TestFunction(V)
u2 = Function(V)

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

F = L**2*(inner(grad(u2), grad(v2))) * dx - v2 * dx
solve(F == 0, u2, bc2, solver_parameters=solver_parameters)
solver_parameters = {
    'snes_type': 'newtonls',
    'snes_monitor': True,
    'ksp_rtol': 1e-4,
        'ksp_type': 'preonly',
        'pc_type': 'lu',
        'pc_factor_mat_solver_packages': 'mumps',
        }



epss = [1000., 500., 200.,  ]
for i, eps in enumerate(epss):
    print("Solving Eikonal with eps == ", float(eps))
    F = inner(sqrt(inner(grad(u2), grad(u2))), v2) * dx - v2 * dx + eps*inner(grad(u2), grad(v2)) * dx
    solve(F == 0, u2, bc2, solver_parameters=solver_parameters)


def zbedf(bathy, distance):
   zbed = conditional(le(abs(distance),75),6.0,
                         conditional(le(abs(distance),250),Max(6.0-distance/175 *2,bathy),bathy))
   return zbed



bath.interpolate(zbedf(bath,u2))
#bath.assign(bath+lat)

mu_manning = Function(V).interpolate(conditional(le(bath,30),0.023, 0.023))
File(outputdir + '/manning.pvd').write(mu_manning)
chk = DumbCheckpoint(inputdir + "/manning2D", mode=FILE_CREATE)
chk.store(mu_manning, name="manning")

File(outputdir + '/bath.pvd').write(bath)
chk = DumbCheckpoint(inputdir + "/bathymetry2D", mode=FILE_CREATE)
chk.store(bath, name="bathymetry")


