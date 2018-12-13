import matplotlib.pyplot as plt
import matplotlib.dates as md
import dateutil

file = open('/data/BODCdata/Hinkley_Point_20030506_20030604.csv', 'r')

for _ in range(34):
    next(file)

datestrings = []
data = []
t=[]
time= 0
for line in file:
    words = line.split(',')
    datestrings.append(words[9])
    data.append(float(words[11]))
    t.append(time)
    time += (15 * 60)/(60*60) #time in hours
    if time > 30*24: #*60*60:
        break




# t=t[:int(432000/(60*60))] #converts to hours and only shows 5 days

# plt.plot(t, data[:len(t)])
# plt.show()



import h5py
from pylab import *
import inputs.input_file_paths

outputdir = inputs.input_file_paths.paraview_output_folder
df=h5py.File(outputdir + '/diagnostic_detectors.hdf5', 'r')

t2=df['time'][:]


eta=df['HINKLEY'][:,0]

t2 = t2[0:432000]/(60*60) # convert to hours, only show 5 days

ion()
close('all')

plot(t2,eta[:len(t2)], t, data[:len(t)])
legend(['HINKLEY'])

show(block=True)