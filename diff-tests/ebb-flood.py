#load elev/free surface perturbation data through time at stations
# check when elev goes from + to -
# if + to - => ebb
#if - to + => flood
import h5py
from pylab import *
#import inputs.input_file_paths
from thetis import *

outputdir = '/data/swansea_2018_copy2/outputs'
df=h5py.File(outputdir + '/diagnostic_detectors.hdf5', 'r')


#print(df['PLYMOUTH'][:]) # [:,0]
# print(df['PLYMOUTH'].shape)
# print(df['PLYMOUTH'][:,0].shape)

#print(df['time'][:])
t=df['time'][:]


eta=df['Avonmouth'][:,0]
eta3=df['Hinkley_Point'][:,0]
eta4=df['Ilfracombe'][:,0]
eta5=df['Mumbles'][:,0]
eta6=df['Newport'][:,0]

tees_flood=[]
tees3_flood=[]
tees4_flood=[]
tees5_flood=[]
tees6_flood=[]

tees_ebb=[]
tees3_ebb=[]
tees4_ebb=[]
tees5_ebb=[]
tees6_ebb=[]

#t = t[0:432000]/(60*60) # convert to hours, only show 5 days

# print (len(t))
for i in range(len(t)-1):
    #check when flood at ech gauge station
    if eta[i+1]>0 and eta[i]<0:
        tees_flood.extend(t[i])

    if eta3[i+1]>0 and eta3[i]<0:
        tees3_flood.extend(t[i])

    if eta4[i+1]>0 and eta4[i]<0:
        tees4_flood.extend(t[i])

    if eta5[i+1]>0 and eta5[i]<0:
        tees5_flood.extend(t[i])

    if eta6[i+1]>0 and eta6[i]<0:
        tees6_flood.extend(t[i])

    # chck when ebb at each gauge station
    if eta[i+1]<0 and eta[i]>0:
        tees_ebb.extend(t[i])

    if eta3[i+1]<0 and eta3[i]>0:
        tees3_ebb.extend(t[i])

    if eta4[i+1]<0 and eta4[i]>0:
        tees4_ebb.extend(t[i])

    if eta5[i+1]<0 and eta5[i]>0:
        tees5_ebb.extend(t[i])

    if eta6[i+1]<0 and eta6[i]>0:
        tees6_ebb.extend(t[i])

# when all station sampled, take the average of these t
# divided by 1000 because the h5 files are exported every 1000 time steps, but now that we kno when we would need data,
# we could just ask the code when it runs to take a pic at those times
av_t_flood = []
for i in range(len(tees_flood)):
    a = int((tees_flood[i] + tees3_flood[i] + tees4_flood[i] + tees5_flood[i] + tees6_flood[i]) / (5*1000))
    av_t_flood.append(a)
print(av_t_flood)

av_t_ebb = []
for i in range(len(tees_ebb)):
    a = int((tees_ebb[i] + tees3_ebb[i] + tees4_ebb[i] + tees5_ebb[i] + tees6_ebb[i]) / (5*1000))
    av_t_ebb.append(a)
print(av_t_ebb)

#print (len(av_t_flood), len(av_t_ebb))
#load the velocity data at these time

file_location = '/data/swansea_2018_copy2/outputs/hdf5'
mesh = Mesh('/data/swansea_2018_copy2/inputs/severn_mesh_2.msh')
xyvector = mesh.coordinates.dat.data

V = VectorFunctionSpace(mesh, "DG", 1)
P1CG = FunctionSpace(mesh, "CG", 1)

field_2d = Function(V, name='uv_2d')
magnitude_ebb = Function(P1CG, name= 'mag')
magnitude_flood = Function(P1CG, name= 'mag')

checkpoint_file = checkpointing.DumbCheckpoint(file_location+
                                               '/Velocity2d_000' +str(av_t_ebb[0]), mode=FILE_READ)
checkpoint_file.load(field_2d)
checkpoint_file.close()

#at the moment do the code on the magnitude, but could try to adapt to get it on the velocity vector too?
ex = VelocityMagnitudeSolver(magnitude_ebb, field_2d, 0,  min_val=0, solver_parameters= {'snes_type': 'newtonls'})
ex.solve()

velocity_ebb_data_set = np.empty((len(av_t_ebb), magnitude_ebb.dat.data.shape[0]))
velocity_flood_data_set = np.empty((len(av_t_ebb), magnitude_ebb.dat.data.shape[0]))

print (velocity_ebb_data_set.shape, velocity_flood_data_set.shape)


for i in range(len(av_t_ebb)): #at every time of ebb
    print('Reading h5 files. Time  ', i, len(range(len(av_t_ebb))))

    if len(str(av_t_ebb[i]))==2:
        # ebb
        checkpoint_file = checkpointing.DumbCheckpoint(file_location +
                                                       '/Velocity2d_000' + str(av_t_ebb[i]), mode=FILE_READ)
        checkpoint_file.load(field_2d)
        checkpoint_file.close()
        ex = VelocityMagnitudeSolver(magnitude_ebb, field_2d, 0, min_val=0, solver_parameters={'snes_type': 'newtonls'})
        ex.solve()
        velocity_ebb_data_set[i, :] = magnitude_ebb.dat.data[:]

        # flood
        checkpoint_file = checkpointing.DumbCheckpoint(file_location +
                                                       '/Velocity2d_000' + str(av_t_flood[i]), mode=FILE_READ)
        checkpoint_file.load(field_2d)
        checkpoint_file.close()
        ex = VelocityMagnitudeSolver(magnitude_flood, field_2d, 0, min_val=0,
                                     solver_parameters={'snes_type': 'newtonls'})
        ex.solve()
        velocity_flood_data_set[i, :] = magnitude_flood.dat.data[:]
    elif len(str(av_t_ebb[i]))==3:
        # ebb
        checkpoint_file = checkpointing.DumbCheckpoint(file_location +
                                                       '/Velocity2d_00' + str(av_t_ebb[i]), mode=FILE_READ)
        checkpoint_file.load(field_2d)
        checkpoint_file.close()
        ex = VelocityMagnitudeSolver(magnitude_ebb, field_2d, 0, min_val=0, solver_parameters={'snes_type': 'newtonls'})
        ex.solve()
        velocity_ebb_data_set[i, :] = magnitude_ebb.dat.data[:]

        # flood
        checkpoint_file = checkpointing.DumbCheckpoint(file_location +
                                                       '/Velocity2d_00' + str(av_t_flood[i]), mode=FILE_READ)
        checkpoint_file.load(field_2d)
        checkpoint_file.close()
        ex = VelocityMagnitudeSolver(magnitude_flood, field_2d, 0, min_val=0,
                                     solver_parameters={'snes_type': 'newtonls'})
        ex.solve()
        velocity_flood_data_set[i, :] = magnitude_flood.dat.data[:]

    elif len(str(av_t_ebb[i]))==4:
        # ebb
        checkpoint_file = checkpointing.DumbCheckpoint(file_location +
                                                       '/Velocity2d_0' + str(av_t_ebb[i]), mode=FILE_READ)
        checkpoint_file.load(field_2d)
        checkpoint_file.close()
        ex = VelocityMagnitudeSolver(magnitude_ebb, field_2d, 0, min_val=0, solver_parameters={'snes_type': 'newtonls'})
        ex.solve()
        velocity_ebb_data_set[i, :] = magnitude_ebb.dat.data[:]

        # flood
        checkpoint_file = checkpointing.DumbCheckpoint(file_location +
                                                       '/Velocity2d_0' + str(av_t_flood[i]), mode=FILE_READ)
        checkpoint_file.load(field_2d)
        checkpoint_file.close()
        ex = VelocityMagnitudeSolver(magnitude_flood, field_2d, 0, min_val=0,
                                     solver_parameters={'snes_type': 'newtonls'})
        ex.solve()
        velocity_flood_data_set[i, :] = magnitude_flood.dat.data[:]
    elif len(str(av_t_ebb[i]))==5:
        # ebb
        checkpoint_file = checkpointing.DumbCheckpoint(file_location +
                                                       '/Velocity2d_000' + str(av_t_ebb[i]), mode=FILE_READ)
        checkpoint_file.load(field_2d)
        checkpoint_file.close()
        ex = VelocityMagnitudeSolver(magnitude_ebb, field_2d, 0, min_val=0, solver_parameters={'snes_type': 'newtonls'})
        ex.solve()
        velocity_ebb_data_set[i, :] = magnitude_ebb.dat.data[:]

        # flood
        checkpoint_file = checkpointing.DumbCheckpoint(file_location +
                                                       '/Velocity2d_' + str(av_t_flood[i]), mode=FILE_READ)
        checkpoint_file.load(field_2d)
        checkpoint_file.close()
        ex = VelocityMagnitudeSolver(magnitude_flood, field_2d, 0, min_val=0,
                                     solver_parameters={'snes_type': 'newtonls'})
        ex.solve()
        velocity_flood_data_set[i, :] = magnitude_flood.dat.data[:]

mean_ebb1 = Function(P1CG, name='mean_ebb_velocity')
mean_ebb = []

mean_flood1 = Function(P1CG, name='mean_ebb_velocity')
mean_flood = []

print(magnitude_ebb.dat.data.shape[0], magnitude_flood.dat.data.shape[0])

for i in range(magnitude_ebb.dat.data.shape[0]):
    # ebb
    thetis_velo_ebb = velocity_ebb_data_set[:, i]
    thetis_velo_ebb = thetis_velo_ebb.mean()

    mean_ebb.append(thetis_velo_ebb)

    # flood
    thetis_velo_flood = velocity_flood_data_set[:, i]
    thetis_velo_flood = thetis_velo_flood.mean()

    mean_flood.append(thetis_velo_flood)

mean_ebb1.dat.data[:] = np.array(mean_ebb)[:]
File('mean_ebb.pvd').write(mean_ebb1)

mean_flood1.dat.data[:] = np.array(mean_flood)[:]
File('mean_flood.pvd').write(mean_flood1)
# take an average and max field? for ebb flows and flood flows?


