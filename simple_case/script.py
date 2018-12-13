from thetis import *
#import utm

mesh2d = Mesh('inputs/mesh.msh')


# create a vectorfunction space for uv
W = VectorFunctionSpace(mesh2d, "DG", 1)
uv = Function(W, name='vel_2d')
print(uv.dat.data.shape)

# create arrays
u_data_set = np.empty((int(t_end/t_export)+1,uv.dat.data.shape[0]))
v_data_set = np.empty((int(t_end/t_export)+1,uv.dat.data.shape[0]))
elev_data_set = np.empty((int(t_end/t_export)+1,uv.dat.data.shape[0]))


count = 0
for i in range(0,int(t_end/t_export)+1): # going through the exported data
      #print('Reading h5 files. Time ',i,i*t_export)
    solver_obj.load_state(i)
    u_data_set[count, :] = solver_obj.fields.uv_2d.dat.data[:, 0]  # extract u componnt at time-step i
    v_data_set[count, :] = solver_obj.fields.uv_2d.dat.data[:, 1]  #extract v component at time-step i
    elev_data_set[count, :] = solver_obj.fields.elev_2d.dat.data[:]  # extract elevation at timestep i


    #print(u_data_set[1], v_data_set[1], u[1])
    count += 1



DG_2d = FunctionSpace(mesh2d, 'DG', 1)
x = SpatialCoordinate(mesh2d)
x_vector, y_vector = interpolate(x[0], Function(DG_2d)).dat.data, interpolate(x[1], Function(DG_2d)).dat.data


#x_0, y_0, utm_zone, zone_letter = utm.from_latlon(x_vector[1], y_vector[1])

print(x_vector[1], x_0, x_vector[1], y_0)