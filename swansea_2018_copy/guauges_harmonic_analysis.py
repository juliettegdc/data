import numpy as np
import inputs.input_file_paths
import datetime
import uptide
import math

# load elevation data of file over month - from csv file

file = open('/data/BODCdata/Newport_20030506_20030604.csv', 'r')

for _ in range(34):
    next(file)

datestrings = []
data = []
t=[]
time= 0

add=0
cnt=0
for line in file:
    words = line.split(',')
    add+= float(words[11])
    cnt+=1

mean= add/cnt
print(mean)

file.close()

file = open('/data/BODCdata/Newport_20030506_20030604.csv', 'r')
for _ in range(34):
    next(file)

datestrings = []
data = []
t=[]
time= 0

for line in file:
    words = line.split(',')
    datestrings.append(words[9])
    data.append(float(words[11])-mean)
    t.append(time)
    time += (15 * 60) #time in seconds, very mesurmnt tkn every 15 min
    if time > 30*24*60*60:
        break

file.close()

# thetis_elev will be the array with the elevations

constituents = ['M2', 'S2', 'Q1', 'O1', 'P1', 'K1', 'N2',  'K2']
tide = uptide.Tides(constituents)
tide.set_initial_time(datetime.datetime(2003, 5, 6, 8, 0))


thetis_amplitudes, thetis_phases = uptide.analysis.harmonic_analysis(tide, data[:], t[:])

M2_amp = thetis_amplitudes[0]
M2_phase = thetis_phases[0]

S2_amp = thetis_amplitudes[1]
S2_phase = thetis_phases[0]

print(M2_amp, M2_phase, S2_amp,S2_phase)




import h5py
from pylab import *
import inputs.input_file_paths

outputdir = inputs.input_file_paths.paraview_output_folder
df=h5py.File( 'outputs_copy/diagnostic_detectors.hdf5', 'r')

t2=df['time'][:]


eta=df['Newport'][:,0]

#t2 = t2[0:432000] # convert to hours, only show 5 days

ion()
close('all')

data = eta[:]


thetis_amplitudes, thetis_phases = uptide.analysis.harmonic_analysis(tide, data[:], t2[:])

M2_amp = thetis_amplitudes[0]
M2_phase = thetis_phases[0]

S2_amp = thetis_amplitudes[1]
S2_phase = thetis_phases[0]

print(M2_amp, M2_phase, S2_amp,S2_phase)
