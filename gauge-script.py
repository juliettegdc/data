import matplotlib.pyplot as plt
import matplotlib.dates as md
import dateutil

"""
file = open('BODCdata/Hinkley_Point_20030506_20030604.csv', 'r')

file.readline()
file.readline()

datestrings = []
data = []
for line in file:
    words = line.split(',')
    datestrings.append(words[9])
    data.append(float(words[11]))

#print (datestrings)
#print(data)

dates = [dateutil.parser.parse(s) for s in datestrings]

#plt.subplots_adjust(bottom=0.2)
#plt.xticks( rotation=25 )

ax=plt.gca()
ax.set_xticks(dates)

xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
ax.xaxis.set_major_formatter(xfmt)
plt.plot(dates, data, "o-")
plt.show()
"""
from thetis import *
import pyproj
import swansea_2018_copy.utm as utm
import numpy as np

"""
tidegauge_file = "gauges.txt"

UTM_ZONE30 = pyproj.Proj(
        proj='utm',
        zone=30,
        datum='WGS84',
        units='m',
        errcheck=True)
LL_WGS84 = pyproj.Proj(proj='latlong', datum='WGS84', errcheck=True)


file = open("gauges.txt", "r")

file.readline()

gauge_names = []
transform = []
gauge_xy = []

for line in file:
        words = line.split(',')
        gauge_names.append(words[0])
        #gauges_latlon.append([float(words[1]), float(words[2])]) #lat, lon
        transform = (utm.from_latlon(float(words[1]), float(words[2])))
        gauge_xy.append([transform[0],transform[1]])

mesh2d = Mesh("swansea_2018_copy2/inputs/severn_mesh_2.msh")

locations, names = select_and_move_detectors(mesh2d, gauge_xy, gauge_names, maximum_distance=10e3)

print (locations, names)
"""
#############################################################################
"""
#run through x and y o find th index at which = to gauge location
mesh2d = Mesh("swansea_2018_copy2/inputs/severn_mesh_2.msh")
P1_2d = FunctionSpace(mesh2d, 'CG', 1)
x,y = SpatialCoordinate(mesh2d)
x_vector, y_vector = interpolate(x, Function(P1_2d)).dat.data, interpolate(y, Function(P1_2d)).dat.data

loc =[501019.24653257156, 5710932.153126809] #[382595.6057414686, 5663883.469540156]#[519779.97957765067, 5706678.7916577]

print(len(x_vector), len(y_vector))

print(np.any(x_vector == 501019.24653257156))
"""
################################
t_start = 0
t_end = t_start + 29 * 24 * 3600
t_n = int((t_end - t_start)/ dt + 1)

thetis_times = t_start + dt * np.arange(t_n) + dt

print(inputs.input_file_paths.paraview_output_folder +
                                               '/elev_-' +str(int(t_start+dt)) + '.0.h5')
checkpoint_file = checkpointing.DumbCheckpoint(inputs.input_file_paths.paraview_output_folder +
                                               '/elev_-' +str(int(t_start+dt)) + '.0', mode=FILE_READ)
checkpoint_file.load(elev)
checkpoint_file.close()
elev_data_set = np.empty((t_n, elev.dat.data.shape[0]))

for i in range(t_n):
    if i==499:
        break
    else :

        print('Reading h5 files. Time  ', i, len(range(t_n)))
        checkpoint_file = checkpointing.DumbCheckpoint(inputs.input_file_paths.paraview_output_folder +
                                                   '/elev_-' + str(t_start+(i+1) * dt), mode=FILE_READ) #i+1
        checkpoint_file.load(elev)
        checkpoint_file.close()
        elev_data_set[i, :] = elev.dat.data[:]

