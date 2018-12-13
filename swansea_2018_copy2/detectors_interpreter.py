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

#eta=df['PLYMOUTH'][:,0]
eta=df['HINKLEY'][:,0]
#eta2=df['ST.CATHERINE'][:,0]
#eta3=df['BARFLEUR'][:,0]
#eta4=df['ST'][:,0]
#eta5=df['ST.'][:,0]
#eta6=df['ST_0'][:,0]
#eta7=df['ST_1'][:,0]

t = t[0:432000]/(60*60) # convert to hours, only show 5 days

ion()
close('all')

plot(t,eta[:len(t)])#,t,eta2[:len(t)],t,eta3[:len(t)],t,eta4[:len(t)],t,eta5[:len(t)],t,eta6[:len(t)],t,eta7[:len(t)]) # plots based off the current length of t, allows you to plot mid simulation when lengths are changing
legend(['HINKLEY']) #Plymouth','St Catherine','Barfleur', 'ST', 'ST.', 'ST_0', 'ST_1'])

show(block=True)


