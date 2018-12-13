import numpy as np
from matplotlib import pyplot as plt
import datetime
import uptide
from firedrake import *
from thetis import *
from netCDF4 import Dataset as NetCDFFile
from scipy.interpolate import interp2d
import utm
import sys
import pyproj
import uptide.tidal_netcdf
import h5py
import matplotlib.pyplot as plt
import matplotlib.colors as mpl

# output_dir = 'comparison'
# create_directory(output_dir)
# # t_end =  2724300 #2377 hdf5 files
# t_end = 2000  # 2377 hdf5 files
# start_file = 0
# t_export = 500
#
# t_n = int(t_end / t_export + 1)
# thetis_times = t_export * np.arange(t_n) + t_export

######
file_location1 = '/data/swansea_2018_copy/outputs'
file_location2 = '/data/swansea_2018_copy2/outputs'


mesh1 = Mesh('/data/swansea_2018_copy/inputs/swansea_2018_4.msh')
mesh2 = Mesh('/data/swansea_2018_copy2/inputs/severn_mesh_2.msh')

xyvector1 = mesh1.coordinates.dat.data  # UTM51, nodes
xyvector2 = mesh2.coordinates.dat.data

P1DG_1 = FunctionSpace(mesh1, "CG", 1)
P1DG_2 = FunctionSpace(mesh2, "CG", 1)

x = xyvector2[:, 0]
y = xyvector2[:, 1]

max_x = 545141.9572129188#5739115.327273086 #344135.7783181509
min_x = 359906.8240550338  #322639.8772734722
max_y = 5741909.44544 #474965.0681569836
min_y = 5621942.00438#458630.5856613767

n = 400
gridxy = np.mgrid[min_x:max_x:400j, min_y:max_y:400j].T
gridxy = np.reshape(gridxy, (n * n, 2))


# n = 400 # sqrt(no. points)
# gridxy = np.mgrid[min(x):max(x):400j, min(y):max(y):400j].T # Ny x Nx x 2
# gridxy = np.reshape(gridxy, (n*n, 2)) # N x 2


gridx, gridy = gridxy.T

field1 = Function(P1DG_1, name='max_velocity')
field1_data_set = np.empty(gridxy.shape[0])

field2 = Function(P1DG_2, name='max_velocity')
field2_data_set = np.empty(gridxy.shape[0])

#for i in range(start_file, int(t_end / t_export) + 1):
#print('Reading h5 files. Time  ', i, i * t_export)

checkpoint_file1 = DumbCheckpoint(file_location1 + '/max_velocity_2592000.0', mode=FILE_READ)
checkpoint_file1.load(field1)
checkpoint_file1.close()
f1 = field1.at(gridxy, dont_raise=True)  # interpolate elev1 at gridxy points
f1 = [np.array([np.nan]) if x == None else x for x in f1]  # thi is a polem wth vector functions
f1 = np.array(f1)
field1_data_set[:] = f1[:]


checkpoint_file2 = DumbCheckpoint(file_location2 + '/max_velocity_2592000.0', mode=FILE_READ)
checkpoint_file2.load(field2)
checkpoint_file2.close()
f2 = field2.at(gridxy, dont_raise=True)  # interpolate elev2 at gridxy points
f2 = [np.array([np.nan]) if x == None else x for x in f2]
f2 = np.array(f2)
field2_data_set[:] = f2[:]


thetis_field1 = []
thetis_field2 = []
thetis_diff_field = []
max_diff = 0
for i in range(field2_data_set.shape[0]):
    thetis_field1.append(field1_data_set[i])
    thetis_field2.append(field2_data_set[i])
    diff_field = (field1_data_set[i] - field2_data_set[i])

    if abs(diff_field)>max_diff:
        max_diff=abs(diff_field)

    # thetis_diff_field.append(diff_field)
    thetis_diff_field.append((abs(diff_field) / max_diff) * 100)

print (max_diff)


field1_plot = np.array(thetis_field1)[:]
field1_plot = field1_plot.reshape(n, n)

field2_plot = np.array(thetis_field2)[:]
field2_plot = field2_plot.reshape(n, n)

fieldd_plot = np.array(thetis_diff_field)[:]
fieldd_plot = fieldd_plot.reshape(n, n)




fig1 = plt.imshow(field1_plot, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
                  cmap='coolwarm')
plt.colorbar()
plt.savefig('velocity_max_swansea')
plt.show()

fig2 = plt.imshow(field2_plot, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
                  cmap='coolwarm')
plt.colorbar()
plt.savefig('velocity_max_plot')
plt.show()

figd = plt.imshow(fieldd_plot, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
                  cmap='rainbow')
plt.colorbar()
plt.savefig('velocity_max_diff_swansea_plot')
plt.show()