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
file_location1 = 'modified/tidal_model/output/hdf5/' #location of the Elevation2d output files
file_location2 = '1995/output/hdf5/'
constituent_order = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', 'M4', 'MS4', 'MN4'] #LEAVE THIS 
constituents = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', 'M4', 'MS4', 'MN4' ]  #CHOSEN CONSTITUENTS MUST BE IN THE ORDER ABOVE

#t_end =  2724300 #2377 hdf5 files
t_end =  450000 #2377 hdf5 files
start_file = 192
t_export = 900

t_n = int(t_end/t_export + 1)
thetis_times = t_export*np.arange(t_n) + t_export

mesh1 = Mesh('modified/tidal_model/runfiles/mesh.msh')
mesh2 = Mesh('1995/runfiles/mesh.msh')

xyvector1 = mesh1.coordinates.dat.data # UTM51, nodes
xyvector2 = mesh2.coordinates.dat.data

P1DG_1 = FunctionSpace(mesh1, "DG", 1)
P1DG_2 = FunctionSpace(mesh2, "DG", 1)

x = xyvector1[:,0]
y = xyvector1[:,1]

max_x = 344135.7783181509
min_x = 322639.8772734722
max_y = 474965.0681569836
min_y = 458630.5856613767

n = 400 
gridxy = np.mgrid[min_x:max_x:400j, min_y:max_y:400j].T 
gridxy = np.reshape(gridxy, (n*n, 2)) 

'''
n = 400 # sqrt(no. points)
gridxy = np.mgrid[min(x):max(x):400j, min(y):max(y):400j].T # Ny x Nx x 2
gridxy = np.reshape(gridxy, (n*n, 2)) # N x 2
'''

gridx, gridy = gridxy.T

elev1 = Function(P1DG_1, name='elev_2d')
elev1_data_set = np.empty((t_n, gridxy.shape[0]))

elev2 = Function(P1DG_2, name='elev_2d')
elev2_data_set = np.empty((t_n, gridxy.shape[0]))

#for i,xy in enumerate(gridxy):
#    for i in range(len(gridxy)):
#        print (xy)

# Reading in elevation data
for i in range(start_file,int(t_end/t_export)+1):
    print('Reading h5 files. Time  ',i,i*t_export)
    checkpoint_file1 = DumbCheckpoint(file_location1 + '/Elevation2d_{:05}'.format(i), mode=FILE_READ)
    checkpoint_file1.load(elev1)
    checkpoint_file1.close()
    f1 = elev1.at(gridxy, dont_raise=True) # interpolate elev1 at gridxy points 
    f1 = [np.array([np.nan]) if x == None else x for x in f1]   
    f1 = np.array(f1)
    elev1_data_set[i, :] = f1[:]

    checkpoint_file2 = DumbCheckpoint(file_location2 + '/Elevation2d_{:05}'.format(i), mode=FILE_READ)
    checkpoint_file2.load(elev2)
    checkpoint_file2.close()
    f2 = elev2.at(gridxy, dont_raise=True) # interpolate elev2 at gridxy points
    f2 = [np.array([np.nan]) if x == None else x for x in f2]   
    f2 = np.array(f2)
    elev2_data_set[i, :] = f2[:]
    
    #empty_array[:] = f1[:] # work around an error
    #f1 = empty_array[:]

    #f1 = np.array(f1) 
    #f1 = f1.reshape(n,n)

    #plt.imshow(f1, origin='lower', cmap='coolwarm')
    #plt.colorbar()
    #plt.show()

detector_amp1 = []
detector_phase1 = []

detector_amp2 = []
detector_phase2 = []

detector_ampd = []
detector_phased = []

print(elev1_data_set.shape)
print('\nHarmonic analysis.\n')
for i in range(elev1_data_set.shape[1]): # loop over each point
    # Mesh1 
    thetis_elev1 = elev1_data_set[:, i]*10 # take elevation of point at all timesteps
    
    tide = uptide.Tides(constituents)
    tide.set_initial_time(datetime.datetime(2018,1,1,0,0))

    thetis_elev1 = thetis_elev1 - thetis_elev1.mean()
    thetis_amplitudes1, thetis_phases1 = uptide.analysis.harmonic_analysis(tide, thetis_elev1[int(start_file):], thetis_times[int(start_file):])
    
    detector_amp1.append(thetis_amplitudes1)
    detector_phase1.append(thetis_phases1)

    # Mesh2
    thetis_elev2 = elev2_data_set[:, i]*10
    
    thetis_elev2 = thetis_elev2 - thetis_elev2.mean()
    thetis_amplitudes2, thetis_phases2 = uptide.analysis.harmonic_analysis(tide, thetis_elev2[int(start_file):], thetis_times[int(start_file):])
 
    detector_amp2.append(thetis_amplitudes2)
    detector_phase2.append(thetis_phases2)

    # Difference
    amp_diff = thetis_amplitudes2 - thetis_amplitudes1
    phase_diff = thetis_phases2 - thetis_phases1
    detector_ampd.append(amp_diff)
    detector_phased.append(phase_diff)

M2_amp_1 = np.array(detector_amp1)[:,constituents.index('M2')]
M2_amp_2 = np.array(detector_amp2)[:,constituents.index('M2')]
M2_amp_d = np.array(detector_ampd)[:,constituents.index('M2')] # shape = N = Nx . Ny

M2_amp_1 = M2_amp_1.reshape(n,n)
M2_amp_2 = M2_amp_2.reshape(n,n)
M2_amp_d = M2_amp_d.reshape(n,n) # shape = Nx x Ny

#M2_amp_sihwa = M2_amp_d[:int(n/2),:int(n/2)]

if COMM_WORLD.rank == 0:

    figd = plt.imshow(M2_amp_d, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower', cmap='coolwarm')
    plt.colorbar()
    plt.savefig('M2_amp_diff')
    plt.show()

    fig1 = plt.imshow(M2_amp_1, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower', cmap='coolwarm')
    plt.colorbar()
    plt.savefig('M2_amp_1')
    plt.show()

    fig2 = plt.imshow(M2_amp_1, extent=[min(gridx), max(gridx), min(gridy), max(gridy)], origin='lower', cmap='coolwarm')
    plt.colorbar()
    plt.savefig('M2_amp_2')
    plt.show()