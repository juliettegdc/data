import numpy as np
import inputs.input_file_paths
import datetime
import uptide
import math
from firedrake import *
from tools.processing_support_scripts import output_field_h5

from netCDF4 import Dataset as NetCDFFile
from scipy.interpolate import RegularGridInterpolator

from matplotlib import pyplot as plt

outputdir = inputs.input_file_paths.paraview_output_folder

mesh2d = Mesh(inputs.input_file_paths.mesh_file)
xvector = mesh2d.coordinates.dat.data

P1_2D = FunctionSpace(mesh2d, "CG", 1)
elev = Function(P1_2D, name='elev_CG')

M2_amp = Function(P1_2D, name='M2_amp')
S2_amp = Function(P1_2D, name='S2_amp')
M2_phase = Function(P1_2D, name='M2_phase')
S2_phase = Function(P1_2D, name='S2_phase')

dt = inputs.input_file_paths.elevation_output_interval

# Set time periods of harmonic analysis
t_start = 0
t_end = t_start + 2000#29 * 24 * 3600
t_n = int((t_end - t_start)/ dt + 1)

thetis_times = t_start + dt * np.arange(t_n) + dt

constituents = ['M2', 'S2', 'Q1', 'O1', 'P1', 'K1', 'N2',  'K2']
# constituents = ['M2', 'S2']

print(inputs.input_file_paths.paraview_output_folder +
                                               '/elev_' +str(int(t_start+dt)) + '.0.h5')
checkpoint_file = checkpointing.DumbCheckpoint(inputs.input_file_paths.paraview_output_folder +
                                               '/elev_' +str(int(t_start+dt)) + '.0', mode=FILE_READ)
checkpoint_file.load(elev)
checkpoint_file.close()
elev_data_set = np.empty((t_n, elev.dat.data.shape[0]))

for i in range(int(t_n)):
    print('Reading h5 files. Time  ', i, len(range(t_n)))
    checkpoint_file = checkpointing.DumbCheckpoint(inputs.input_file_paths.paraview_output_folder +
                                                   '/elev_' + str(t_start+(i + 1) * dt), mode=FILE_READ)
    checkpoint_file.load(elev)
    checkpoint_file.close()
    elev_data_set[i, :] = elev.dat.data[:]

    # if i+1==4:
    #     break
detector_amplitudes = []
detector_phases = []

for i in range(elev.dat.data.shape[0]):
    thetis_elev = elev_data_set[:, i]

    tide = uptide.Tides(constituents)
    tide.set_initial_time(datetime.datetime(2003, 5, 6, 8, 0))

    # Subtract mean
    thetis_elev = thetis_elev - thetis_elev.mean()
    thetis_amplitudes, thetis_phases = uptide.analysis.harmonic_analysis(tide, thetis_elev[:], thetis_times[:])

    detector_amplitudes.append(thetis_amplitudes)
    detector_phases.append(thetis_phases)

M2_amp.dat.data[:] = np.array(detector_amplitudes)[:, 0]
M2_phase.dat.data[:] = np.array(detector_phases)[:, 0]

S2_amp.dat.data[:] = np.array(detector_amplitudes)[:, 1]
S2_phase.dat.data[:] = np.array(detector_phases)[:, 1]


File('outputs/amp.pvd').write(M2_amp, S2_amp)
File('outputs/phase.pvd').write(M2_phase, S2_phase)
# M2_phase.dat.data[:] = np.arcsin(np.sin(M2_phase.dat.data[:]))
M2_phase.dat.data[:] = np.remainder(M2_phase.dat.data[:],2*math.pi)*360/(2*math.pi)
# S2_phase.dat.data[:] = np.arcsin(np.sin(S2_phase.dat.data[:]))
S2_phase.dat.data[:] = np.remainder(S2_phase.dat.data[:],2*math.pi)*360/(2*math.pi)
File('outputs/phase_mod_pi.pvd').write(M2_phase, S2_phase)


output_field_h5(outputdir,M2_amp,"M2_amp")
output_field_h5(outputdir,M2_phase,"M2_phase")
output_field_h5(outputdir,S2_amp,"S2_amp")
output_field_h5(outputdir,S2_phase,"S2_phase")

#
# import utm
#
# # Read TPXO dataset to compare to
# hRe = NetCDFFile('netcdf/hf.ES2008.nc').variables['hRe'][:]
# hIm = NetCDFFile('netcdf/hf.ES2008.nc').variables['hIm'][:]
# lon = NetCDFFile('netcdf/hf.ES2008.nc').variables['lon_z'][:]
# lat = NetCDFFile('netcdf/hf.ES2008.nc').variables['lat_z'][:]
#
# TPXO_M2_amplitude = np.sqrt(np.square(hRe) + np.square(hIm))
# TPXO_M2_amp_interp = RegularGridInterpolator((lon[:, 0], lat[0, :]), TPXO_M2_amplitude[0, :, :])
#
# # Interpolate onto mesh
# TPXO_M2 = Function(P1_2D)
# TPXO_M2_data = TPXO_M2.dat.data
# for i, xy in enumerate(xvector):
#     lat, lon = utm.to_latlon(xy[0], xy[1], 30, 'U', strict=False)
#     try:
#         TPXO_M2_data[i] = TPXO_M2_amp_interp((lon, lat))
#     except ValueError:
#         TPXO_M2_data[i] = 0.0
#
# File('tpxo_amp.pvd').write(TPXO_M2)
#
# abs_diff = Function(P1_2D)
# abs_diff.dat.data[:] = TPXO_M2_data - M2_amp.dat.data[:]
#
# File('absolute_difference.pvd').write(abs_diff)
#
# rel_diff = Function(P1_2D)
# rel_diff.dat.data[:] = (TPXO_M2_data - M2_amp.dat.data[:]) / np.maximum(TPXO_M2_data, 0.01)
#
# File('relative_difference.pvd').write(rel_diff)