# Mesh
#mesh_file = "inputs/severn_mesh_2.msh"  # Elevation boundaries - 1
mesh_file = "inputs/swansea_2018_2.msh"  # with lagoon, Elevation boundaries - 5


# Bathymetry
bathymetry_file = "model_data/severn2.nc"

# Forcing
grid_forcing_file = 'model_data/gridES2008.nc'
hf_forcing_file = 'model_data/hf.ES2008.nc'
range_forcing_coords = ((-5.,-3.),(50,52))

# Detectors
tidegauge_file = 'inputs/useful_gauges.csv'
elevation_detectors = [[432596, 5713648],
                       [422307, 5674612],
                       [490627, 5674006],
                       [501829, 5708798],
                       [519445, 5706866]]

# Outputs folder
paraview_output_folder = 'outputs'

# Constituent field calculation output
spatial_harmonics_distribution = True
elevation_output_interval = 500.

