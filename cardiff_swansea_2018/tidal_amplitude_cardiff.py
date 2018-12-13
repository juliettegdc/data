import utm
import uptide
import uptide.tidal_netcdf
import datetime
import numpy
import scipy.interpolate
import inputs.input_file_paths

utm_zone = 30
utm_band = 'V'

# def initial_forcing(t_start):
constituents = ['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2']
tide = uptide.Tides(constituents)
tide.set_initial_time(datetime.datetime(2003,5,6,8,0))
tnci = uptide.tidal_netcdf.OTPSncTidalInterpolator(tide,
                                                   inputs.input_file_paths.grid_forcing_file,
                                                   inputs.input_file_paths.hf_forcing_file,
                                                   ranges=inputs.input_file_paths.range_forcing_coords)

def get_lowest_astronomical_tide(elev):
  amp = numpy.sqrt(tnci.real_part**2 + tnci.imag_part**2)
  val = amp.sum(axis=0)
  tnci.interpolator = uptide.netcdf_reader.Interpolator(tnci.nci.origin, tnci.nci.delta, val, tnci.nci.mask)
  mesh2d = elev.function_space().mesh()
  xvector = mesh2d.coordinates.dat.data
  evector = elev.dat.data
  data = numpy.loadtxt('inputs/lat.txt')

  intp = scipy.interpolate.NearestNDInterpolator(data[:,0:2], data[:,2])
  for i,xy in enumerate(xvector):
    evector[i] = intp(xy)
