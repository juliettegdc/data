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

n = 800
gridxy = np.mgrid[min_x:max_x:800j, min_y:max_y:800j].T
gridxy = np.reshape(gridxy, (n * n, 2))


# n = 400 # sqrt(no. points)
# gridxy = np.mgrid[min(x):max(x):400j, min(y):max(y):400j].T # Ny x Nx x 2
# gridxy = np.reshape(gridxy, (n*n, 2)) # N x 2


gridx, gridy = gridxy.T

field1 = Function(FunctionSpace(mesh1, 'CG', 1), name='average_shear')
field1_data_set = np.empty(gridxy.shape[0])

field2 = Function(FunctionSpace(mesh2, 'CG', 1), name='average_shear')
field2_data_set = np.empty(gridxy.shape[0])

#for i in range(start_file, int(t_end / t_export) + 1):
#print('Reading h5 files. Time  ', i, i * t_export)

checkpoint_file1 = DumbCheckpoint(file_location1 + '/average_shear_test_2592000.0', mode=FILE_READ)
checkpoint_file1.load(field1)
checkpoint_file1.close()
f1 = field1.at(gridxy, dont_raise=True)  # interpolate elev1 at gridxy points
f1 = [np.array([np.nan]) if x == None else x for x in f1]  # thi is a polem wth vector functions
f1 = np.array(f1)
field1_data_set[:] = f1[:]


checkpoint_file2 = DumbCheckpoint(file_location2 + '/average_shear_test_2592000.0', mode=FILE_READ)
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
    # diff_field = field1_data_set[i] - field2_data_set[i]
    # thetis_diff_field.append(diff_field)
    # if abs(diff_field)>max_diff:
    #     max_diff=abs(diff_field)

thetis_binned1 = np.empty(len(thetis_field1))
thetis_binned2 = np.empty(len(thetis_field1))


for i in range(len(thetis_field1)):
    if float(field1_data_set[i])<0.194:
        thetis_binned1[i] = 1
    elif 0.194 <=float(field1_data_set[i]) < 0.27:
        thetis_binned1[i] = 2
    elif 0.27<= float(field1_data_set[i])<1.26:
        thetis_binned1[i] = 3
    elif 1.26<=float(field1_data_set[i])<5.7:
        thetis_binned1[i] = 4
    elif 5.7<= float(field1_data_set[i])<12.2:
        thetis_binned1[i] = 5
    elif 12.2<= float(field1_data_set[i])<26.0:
        thetis_binned1[i] = 6
    elif float(field1_data_set[i])>=26.0:
        thetis_binned1[i] = 7
    else:
        thetis_binned1[i] = np.nan#field1_data_set[i]

for i in range(len(thetis_field2)):
    if float(field2_data_set[i])<0.194:
        thetis_binned2[i] = 1
    elif 0.194 <= float(field2_data_set[i]) < 0.27:
        thetis_binned2[i] = 2
    elif 0.27 <= float(field2_data_set[i])<1.26:
        thetis_binned2[i] = 3
    elif 1.26 <= float(field2_data_set[i]) < 5.7:
        thetis_binned2[i] = 4
       # print ('here')
    elif 5.7 <= float(field2_data_set[i])<12.2:
        thetis_binned2[i] = 5
    elif 12.2 <= float(field2_data_set[i])<26.0:
        thetis_binned2[i] = 6
    elif float(field2_data_set[i]) >= 26.0:
        thetis_binned2[i] = 7
    else:
        thetis_binned2[i] = np.nan #field2_data_set[i]

    # diff_field = abs(field1_data_set[i] - field2_data_set[i])
    # thetis_diff_field.append((diff_field/max_diff)*100)

thetis_binned_diff = np.empty(len(thetis_field1))

for i in range(len(thetis_binned2)):

    if thetis_binned2[i] >= 1:
        diff_field = thetis_binned2[i] - thetis_binned1[i]
        if diff_field == 0:
            thetis_binned_diff[i] = 0
        elif diff_field > 0:
            thetis_binned_diff[i] = -1  # means new field has smaller magnitudes -> decrease hence -
        elif diff_field < 0:
            thetis_binned_diff[i] = 1  # new field bigger magnitudes -> increase hence +
    else:
        thetis_binned_diff[i] = np.nan



    #thetis_binned_diff.append(diff_field)



field1_plot = np.array(thetis_field1)[:]
field1_plot = field1_plot.reshape(n, n)

fieldbinned1_plot = np.array(thetis_binned1)[:]
fieldbinned1_plot = fieldbinned1_plot.reshape(n, n)

field2_plot = np.array(thetis_field2)[:]
field2_plot = field2_plot.reshape(n, n)

fieldbinned2_plot = np.array(thetis_binned2)[:]
fieldbinned2_plot = fieldbinned2_plot.reshape(n, n)

# fieldd_plot = np.array(thetis_diff_field)[:]
# fieldd_plot = fieldd_plot.reshape(n, n)

fieldd_binned_plot = np.array(thetis_binned_diff)[:]
fieldd_binned_plot = fieldd_binned_plot.reshape(n, n)



#
cmap = mpl.ListedColormap(["midnightblue", "darkgreen", "palegreen", "gold", "darkorange", 'sienna', 'darkred'])
bounds=[1, 2, 3, 4, 5, 6, 7,8]

norm = mpl.BoundaryNorm(bounds, cmap.N)

fig11 = plt.imshow(fieldbinned1_plot, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
                   cmap = cmap, norm = norm)

plt.colorbar(fig11,cmap=cmap,
                norm=norm,boundaries=bounds,ticks=[1, 2, 3, 4, 5, 6, 7])

plt.savefig('bed_shear_shear_mean_swansea')
plt.show()



cmap = mpl.ListedColormap(["midnightblue", "darkgreen", "palegreen", "gold", "darkorange", 'sienna', 'darkred'])
bounds=[1, 2, 3, 4, 5, 6, 7,8]
# cmap = mpl.ListedColormap(["midnightblue" , "gold", 'darkred'])
# bounds=[0,3, 4,7]
norm = mpl.BoundaryNorm(bounds, cmap.N)

fig22 = plt.imshow(fieldbinned2_plot, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
                   cmap = cmap, norm = norm)

plt.colorbar(fig22,cmap=cmap,
                norm=norm,boundaries=bounds,ticks=[1, 2, 3, 4, 5, 6, 7])

plt.savefig('bed_shear_mean')
plt.show()



# cmap = mpl.ListedColormap(["blue", "darkgray", "red"])
# bounds=[-1,1]
# norm = mpl.BoundaryNorm(bounds, cmap.N)

# cmap = plt.cm.get_cmap('seismic', 3)  # mpl.ListedColormap(["darkgreen", "limegreen", "violet", "orange", "darkorange", 'navy', 'purple', 'red'])
# bounds=[-1,0,1]
# norm = mpl.BoundaryNorm(bounds, cmap.N)
#
# figd = plt.imshow(fieldd_binned_plot, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
#                   cmap=cmap, norm=norm)
# plt.colorbar(figd, cmap=cmap,
#                 norm=norm,boundaries=bounds,ticks=[-1,0,1])
# plt.savefig('fieldd_binned_plot')
# plt.show()


# fig1 = plt.imshow(fieldbinned2_plot, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
#                   cmap='coolwarm')
# plt.colorbar()
# plt.savefig('field1_plot')
# plt.show()

# #diff binned
# cmap = plt.cm.get_cmap('seismic', 11)  # mpl.ListedColormap(["darkgreen", "limegreen", "violet", "orange", "darkorange", 'navy', 'purple', 'red'])
# bounds=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# norm = mpl.BoundaryNorm(bounds, cmap.N)
#
# figdd = plt.imshow(fieldd_binned_plot, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
#                    cmap = cmap, norm = norm)
#
# plt.colorbar(figdd,cmap=cmap,
#                 norm=norm,boundaries=bounds,ticks=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
#
# plt.savefig('fieldbinned_diff_plot')
# plt.show()


# fig1 = plt.imshow(field1_plot, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
#                   cmap='coolwarm')
# plt.colorbar()
# plt.savefig('field1_plot')
# plt.show()
#
#
#
# fig2 = plt.imshow(field2_plot, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
#                   cmap='coolwarm')
# plt.colorbar()
# plt.savefig('field2_plot')
# plt.show()
#
#
# cmap = plt.cm.get_cmap('seismic', 16)  # mpl.ListedColormap(["darkgreen", "limegreen", "violet", "orange", "darkorange", 'navy', 'purple', 'red'])
# bounds=[0,1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# norm = mpl.BoundaryNorm(bounds, cmap.N)
#
# figd = plt.imshow(fieldd_plot, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
#                    cmap = cmap, norm = norm)
#
# plt.colorbar(figd,cmap=cmap,
#                 norm=norm,boundaries=bounds,ticks=[0,1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
#
# plt.savefig('field_diff_plot')
# plt.show()
#
figd = plt.imshow(fieldd_binned_plot, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
                  cmap='seismic')
plt.colorbar()
plt.savefig('bd_shear_mean_diff_swansea')
plt.show()