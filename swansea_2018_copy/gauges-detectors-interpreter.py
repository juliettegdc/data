import h5py
from pylab import *
import inputs.input_file_paths

outputdir = inputs.input_file_paths.paraview_output_folder
df=h5py.File(outputdir + '/diagnostic_detectors.hdf5', 'r')

df.keys()
print (list(df.keys()))

#print(df['PLYMOUTH'][:]) # [:,0]
# print(df['PLYMOUTH'].shape)
# print(df['PLYMOUTH'][:,0].shape)

#print(df['time'][:])
t=df['time'][:]


eta=df['Avonmouth'][:,0]
eta2=df['Dover'][:,0]
eta3=df['Hinkley_Point'][:,0]
eta4=df['Ilfracombe'][:,0]
eta5=df['Mumbles'][:,0]
eta6=df['Newport'][:,0]


t = t[0:36000]/(60*60) # convert to hours, only show 5 days
#432000
ion()
close('all')

plot(t,eta[:len(t)],t,eta2[:len(t)],t,eta3[:len(t)],t,eta4[:len(t)],t,eta5[:len(t)],t,eta6[:len(t)]) # plots based off the current length of t, allows you to plot mid simulation when lengths are changing
legend(['Avonmouth', 'Dover', 'Hinkley_Point', 'Ilfracombe', 'Mumbles', 'Newport'])

show(block=True)


