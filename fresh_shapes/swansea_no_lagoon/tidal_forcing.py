import utm
import uptide
import uptide.tidal_netcdf
import datetime
import inputs.input_file_paths
from math import tanh

utm_zone = 30
utm_band = 'V'

# def initial_forcing(t_start):
constituents = ['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2']
tide = uptide.Tides(constituents)
tide.set_initial_time(datetime.datetime(2003,5,6,8,0))   # year, month, day, hour, minute, second
tnci = uptide.tidal_netcdf.OTPSncTidalInterpolator(tide,
     inputs.input_file_paths.grid_forcing_file,
     inputs.input_file_paths.hf_forcing_file, ranges=inputs.input_file_paths.range_forcing_coords)

def set_tidal_field(elev, t):
  tnci.set_time(t)
  mesh2d = elev.function_space().mesh()
  xvector = mesh2d.coordinates.dat.data
  evector = elev.dat.data
  for i,xy in enumerate(xvector):
    lat, lon = utm.to_latlon(xy[0], xy[1], utm_zone, utm_band)
    try:
      evector[i] = tnci.get_val((lon, lat))    # Adding initial a correction depth for LAT
    except uptide.netcdf_reader.CoordinateError:
      evector[i] = 0.    # Adding initial a correction depth for LAT



