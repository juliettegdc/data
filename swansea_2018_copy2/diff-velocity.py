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

output_dir = 'comparison'
create_directory(output_dir)
file_location1 = '/data/swansea_2018_copy/outputs/hdf5'  # location of the Elevation2d output files
file_location2 = 'outputs/hdf5'

# t_end =  2724300 #2377 hdf5 files
t_end = 2000  # 2377 hdf5 files
start_file = 0
t_export = 500

t_n = int(t_end / t_export + 1)
thetis_times = t_export * np.arange(t_n) + t_export

mesh1 = Mesh('/data/swansea_2018_copy/inputs/swansea_2018_2.msh')
mesh2 = Mesh('/data/swansea_2018_copy2/inputs/severn_mesh_2.msh')

xyvector1 = mesh1.coordinates.dat.data  # UTM51, nodes
xyvector2 = mesh2.coordinates.dat.data

P1DG_1 = VectorFunctionSpace(mesh1, "DG", 1)
P1DG_2 = VectorFunctionSpace(mesh2, "DG", 1)

x = xyvector1[:, 0]
y = xyvector1[:, 1]

max_x = 545141.9572129188#5739115.327273086 #344135.7783181509
min_x = 359906.8240550338  #322639.8772734722
max_y = 5741909.44544 #474965.0681569836
min_y = 5621942.00438#458630.5856613767

n = 400
gridxy = np.mgrid[min_x:max_x:400j, min_y:max_y:400j].T
gridxy = np.reshape(gridxy, (n * n, 2))

'''
n = 400 # sqrt(no. points)
gridxy = np.mgrid[min(x):max(x):400j, min(y):max(y):400j].T # Ny x Nx x 2
gridxy = np.reshape(gridxy, (n*n, 2)) # N x 2
'''

gridx, gridy = gridxy.T

######################
uv1 = Function(P1DG_1, name='uv_2d') #is th the write way for a vector?
u1_data_set = np.empty((t_n,uv1.dat.data.shape[0]))
v1_data_set = np.empty((t_n,uv1.dat.data.shape[0]))
# uv2 = Function(P1DG_2, name='uv_2d')
# uv2_data_set = np.empty((t_n, gridxy.shape)) # pb

# for i,xy in enumerate(gridxy):
#    for i in range(len(gridxy)):
#        print (xy)

# make an average velocity out of this from hdf5 files - if works, no need  have ut in the code
# also could just export bss and do average and max here too.
for i in range(start_file, int(t_end / t_export) + 1):

    print('Reading h5 files. Time  ', i, i * t_export)
    checkpoint_file1 = DumbCheckpoint(file_location1 + '/Velocity2d_{:05}'.format(i), mode=FILE_READ)
    checkpoint_file1.load(uv1)
    checkpoint_file1.close()
    #print (uv1.dat.data[:,0])
    f1 = uv1.at(gridxy, dont_raise=True)  # interpolate uv1 at gridxy points
    f1 = [np.array([np.nan]) if x == None else x for x in f1] # doesn't work fr 2d array
    f1 = np.array(f1)
    # u1_data_set[i, :] = f1[:]
    print(f1.shape)
    #uv1_data_set[i, :, :] = f1[:,:] # problem

    # checkpoint_file2 = DumbCheckpoint(file_location2 + '/Velocity2d_{:05}'.format(i), mode=FILE_READ)
    # checkpoint_file2.load(elev2)
    # checkpoint_file2.close()
    # f2 = elev2.at(gridxy, dont_raise=True)  # interpolate elev2 at gridxy points
    # f2 = [np.array([np.nan]) if x == None else x for x in f2]
    # f2 = np.array(f2)
    # elev2_data_set[i, :] = f2[:]



# detector_amp1 = []
# detector_phase1 = []
#
# detector_amp2 = []
# detector_phase2 = []
#
# detector_ampd = []
# detector_phased = []

#print(elev1_data_set.shape)
#print('\nHarmonic analysis.\n')
# for i in range(elev1_data_set.shape[1]):  # loop over each point
#     # Mesh1
#     thetis_elev1 = elev1_data_set[:, i] * 10  # take elevation of point at all timesteps
#
#     tide = uptide.Tides(constituents)
#     tide.set_initial_time(datetime.datetime(2018, 1, 1, 0, 0))
#
#     thetis_elev1 = thetis_elev1 - thetis_elev1.mean()
#     thetis_amplitudes1, thetis_phases1 = uptide.analysis.harmonic_analysis(tide, thetis_elev1[int(start_file):],
#                                                                            thetis_times[int(start_file):])
#
#     detector_amp1.append(thetis_amplitudes1)
#     detector_phase1.append(thetis_phases1)
#
#     # Mesh2
#     thetis_elev2 = elev2_data_set[:, i] * 10
#
#     thetis_elev2 = thetis_elev2 - thetis_elev2.mean()
#     thetis_amplitudes2, thetis_phases2 = uptide.analysis.harmonic_analysis(tide, thetis_elev2[int(start_file):],
#                                                                            thetis_times[int(start_file):])
#
#     detector_amp2.append(thetis_amplitudes2)
#     detector_phase2.append(thetis_phases2)
#
#     # Difference
#     amp_diff = thetis_amplitudes2 - thetis_amplitudes1
#     phase_diff = thetis_phases2 - thetis_phases1
#     detector_ampd.append(amp_diff)
#     detector_phased.append(phase_diff)
#
# M2_amp_1 = np.array(detector_amp1)[:, constituents.index('M2')]
# M2_amp_2 = np.array(detector_amp2)[:, constituents.index('M2')]
# M2_amp_d = np.array(detector_ampd)[:, constituents.index('M2')]  # shape = N = Nx . Ny
#
# M2_amp_1 = M2_amp_1.reshape(n, n)
# M2_amp_2 = M2_amp_2.reshape(n, n)
# M2_amp_d = M2_amp_d.reshape(n, n)  # shape = Nx x Ny

# M2_amp_sihwa = M2_amp_d[:int(n/2),:int(n/2)]

# if COMM_WORLD.rank == 0:
#     figd = plt.imshow(M2_amp_d, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
#                       cmap='coolwarm')
#     plt.colorbar()
#     plt.savefig('M2_amp_diff')
#     plt.show()
#
#     fig1 = plt.imshow(M2_amp_1, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
#                       cmap='coolwarm')
#     plt.colorbar()
#     plt.savefig('M2_amp_1')
#     plt.show()
#
#     fig2 = plt.imshow(M2_amp_1, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower',
#                       cmap='coolwarm')
#     plt.colorbar()
#     plt.savefig('M2_amp_2')
#     plt.show()
#