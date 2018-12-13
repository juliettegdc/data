from thetis import*
from firedrake import*
import numpy as np
import matplotlib.pyplot as plt


#creating a mesh using numpy mgrid, we will be using this grid to compare the results of the simulations

#hardcode the coordinate dimensions, because we are using more than one processes
max_x = 545141.9572129188#5739115.327273086 #344135.7783181509
min_x = 359906.8240550338  #322639.8772734722
max_y = 5741909.44544 #474965.0681569836
min_y = 5621942.00438#458630.5856613767

n = 400 #ensure adequate mesh points to refine our solution
gridxy = np.mgrid[min_x:max_x:400j, min_y:max_y:400j].T 
gridxy = np.reshape(gridxy, (n*n, 2))
gridx, gridy = gridxy.T


###############################Simulation A###########################################

#Firt we define the meshes and the FunctionSpaces for the Simulation A

mesh_1 = Mesh('/data/swansea_2018_copy/inputs/swansea_2018_2.msh')
P1_2d_1 = FunctionSpace(mesh_1, 'DG', 1) #for the components of velocity and magnitude
V_1 = VectorFunctionSpace(mesh_1, 'DG', 1) #for the velocity


#Next we load the velocities of simulation A

# with timed_stage('loading mean velocity of simulation A'):
# 	checkpoint_file1 = DumbCheckpoint('/data/swansea_2018_copy/outputs/average_velocity_test_2592000.0', mode=FILE_READ)
# 	uv_init_1 = Function(P1_2d_1, name = 'average_velocity') #the name 'uv_2d' is specified because thetis saved the velocity in this field name during the simulation
# 	checkpoint_file1.load(uv_init_1)

with timed_stage('loading velocity of simulation A'):
	checkpoint_file1 = DumbCheckpoint('/data/swansea_2018_copy/outputs/hdf5/Velocity2d_02000', mode=FILE_READ)
	uv_init_1 = Function(V_1, name = 'uv_2d') #the name 'uv_2d' is specified because thetis saved the velocity in this field name during the simulation
	checkpoint_file1.load(uv_init_1)
'''If one wants to access the compents of the velocity as functions :

x_vel = Function(P1_2d).interpolate(uv_init[0])
y_velocity = Function(P1_2d).interpolate(uv_init[1])

'''

#Now we use the  the velocity magnitude solver built in thetis to solve for the magnitude.

with timed_stage('solving the magnitude'):
	mag_1 = Function(P1_2d_1, name = 'mag')
	ex_1 = VelocityMagnitudeSolver(mag_1, uv_init_1, 0,  min_val=0, solver_parameters= {'snes_type': 'newtonls'})
	ex_1.solve()

'''
If one wants to store the magnitudes as an hdf file which can then be loaded back for processing : 

chk = DumbCheckpoint('Magnitude_1', mode=FILE_CREATE)
chk.store(mag_1, name='mag')

'''

# Now we use the firedrake point evaluation function (at) to evaluate the values of the 
# magnitude function at our grid points, note that we hardcoded ou mesh coordinates to make our np grid  
# because the point evaluation requries the same set of points to be supplied on each process , check # 612 firedrake for details

f1 = mag_1.at(gridxy, dont_raise=True) # dont_raise = True is used to ensure that points outside the domain are ignored, other wise results in errors
f1 = [np.array([np.nan]) if x == None else x for x in f1] # This will replace all None 
f1 = np.array(f1) # Make our list to a numpy array
f1 = f1.astype(np.float64) # Unless the array is converted to float 64, an error will be produced becase matplotlib canot convert the data to float
#f1 = f1.reshape(n,n) #Reshape our array to match the mesh if you want to make a plot


###############################Simulation B###########################################

# Create Function Spaces

mesh_2 = Mesh('/data/swansea_2018_copy2/inputs/severn_mesh_2.msh')
P1_2d_2 = FunctionSpace(mesh_2, 'DG', 1) 
V_2 = VectorFunctionSpace(mesh_2, 'DG', 1) 


#Next we load the velocities of simulation B

# with timed_stage('loading mean velocity of simulation B'):
# 	checkpoint_file2 = DumbCheckpoint('/data/swansea_2018_copy2/outputs/average_velocity_test_2592000.0', mode=FILE_READ)
# 	uv_init_2 = Function(P1_2d_2, name = 'average_velocity') #the name 'uv_2d' is specified because thetis saved the velocity in this field name during the simulation
# 	checkpoint_file2.load(uv_init_2)

with timed_stage('loading velocity of simulation B'):
	checkpoint_file2 = DumbCheckpoint('/data/swansea_2018_copy2/outputs/hdf5/Velocity2d_02000', mode=FILE_READ)
	uv_init_2 = Function(V_2, name = 'uv_2d') #the name 'uv_2d' is specified because thetis saved the velocity in this field name during the simulation
	checkpoint_file2.load(uv_init_2)
# Now we use the the velocity magnitude solver built in thetis to solve for the magnitude.

with timed_stage('solving the magnitude'):
	mag_2 = Function(P1_2d_2, name = 'mag')
	ex_2 = VelocityMagnitudeSolver(mag_2, uv_init_2, 0,  min_val=0, solver_parameters= {'snes_type': 'newtonls'})
	ex_2.solve()


# Now we use the firedrake point evaluation function (at) to evaluate the values of the magnitude at grid points

f2 = mag_2.at(gridxy, dont_raise=True)
f2 = [np.array([np.nan]) if x == None else x for x in f2] 
f2 = np.array(f2) 
f2 = f2.astype(np.float64) #what is this?
#f2 = f2.reshape(n,n) 

#####################################Make the Diff##############################

#Now make the diff of our arrays f1 and f2

f3 = f1 - f2
f3 = f3.reshape(n,n)

############################Create the Plot#########################################


if COMM_WORLD.rank == 0: # Using more then one processes will fail
	plt.imshow(f3, extent = [min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower', cmap='coolwarm')
	plt.colorbar()
	plt.show()
