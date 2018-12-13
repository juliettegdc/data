import inputs.input_file_paths
from thetis import *
import pyproj
import utm
import numpy as np
import inputs.input_file_paths
import datetime
import uptide
import math
from tools.processing_support_scripts import output_field_h5


tidegauge_file = inputs.input_file_paths.tidegauge_file

UTM_ZONE30 = pyproj.Proj(
        proj='utm',
        zone=30,
        datum='WGS84',
        units='m',
        errcheck=True)
LL_WGS84 = pyproj.Proj(proj='latlong', datum='WGS84', errcheck=True)


def get_detectors(mesh2d):
    """
    gauge_names = np.loadtxt(tidegauge_file, skiprows=1, usecols=(0,), dtype=str, delimiter=',')
    gauge_xy = np.loadtxt(tidegauge_file, skiprows=1, usecols=(3,4), delimiter=',')
    ind = np.argsort(gauge_names)
    gauge_names = list(gauge_names[ind])
    gauge_xy = list(gauge_xy[ind])

    #make names unique
    unique_names = []
    last_name = ''; ctr = 0
    for name in gauge_names:
        if name==last_name:
            unique_names.append(name + '_' + str(ctr))
            ctr += 1
        else:
            unique_names.append(name)
            ctr = 0
        last_name = name
    """
    # add custom detectors examples:
    # gauge_latlon.append([58.659477538216585, -3.141400563744444])
    # unique_names.append('Exeter2013')
    # xy = [pyproj.transform(LL_WGS84, crs, lon, lat) for lat, lon in gauge_latlon]


    mesh2d = Mesh(inputs.input_file_paths.mesh_file)
    P1_2d = FunctionSpace(mesh2d, 'CG', 1)

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
    gauges_latlon = []

    #i = 0
    for line in file:
        words = line.split(',')
        gauge_names.append(words[0])
        gauges_latlon.append([float(words[1]), float(words[2])])  # lat, lon
        transform = (utm.from_latlon(float(words[1]), float(words[2])))
        gauge_xy.append([transform[0], transform[1]])

        #print(gauge_names[i], gauges_latlon[i], gauge_xy[i])
        #i += 1


    return select_and_move_detectors(mesh2d, gauge_xy, gauge_names, maximum_distance=10e3)

if __name__ == "__main__":
    mesh2d = Mesh(inputs.input_file_paths.mesh_file)

    locations, names = get_detectors(mesh2d)
    if mesh2d.comm.rank == 0: # only processor 0
        print_output("Found detectors: {}".format(names))

        # write out shape-file
        import shapely.geometry
        import fiona
        import fiona.crs

        schema = {'geometry': 'Point', 'properties': {'name': 'str'}}
        crs = fiona.crs.from_string(UTM_ZONE30.srs)
        with fiona.collection("data/detectors.shp", "w", "ESRI Shapefile", schema, crs=crs) as output:
            for xy, name in zip(locations, names):
                point = shapely.geometry.Point(xy[0], xy[1])
                output.write({'properties': {'name': name}, 'geometry': shapely.geometry.mapping(point)})


