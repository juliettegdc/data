import utm
from netCDF4 import Dataset as NetCDFFile
import scipy.interpolate
from firedrake import *

minimum_depth = -10.0

utm_zone=30
utm_band='V'

def get_val(X,t):
  lat, lon = utm.to_latlon(X[0], X[1], utm_zone, utm_band)
  return max(minimum_depth, -nci.get_val((lat, lon)))

def get_bathymetry(bathymetry_file, mesh2d):
  nc = NetCDFFile(bathymetry_file)
  lat = nc.variables['lat'][:]
  lon = nc.variables['lon'][:]
  values = nc.variables['z'][:,:]
  values = values.filled(9999.)
  interpolator = scipy.interpolate.RegularGridInterpolator((lat, lon), values)
  P1_2d = FunctionSpace(mesh2d, 'CG', 1)
  bathymetry2d = Function(P1_2d, name="bathymetry")
  xvector = mesh2d.coordinates.dat.data
  bvector = bathymetry2d.dat.data
  assert xvector.shape[0]==bvector.shape[0]
  for i,xy in enumerate(xvector):
      lat, lon = utm.to_latlon(xy[0], xy[1], utm_zone, utm_band)
      bvector[i] = max(-interpolator((lat, lon)), minimum_depth)

#  bathymetry2d.assign(50)

  return bathymetry2d

def smoothen_bathymetry(bathymetry2d):
  v = TestFunction(bathymetry2d.function_space())
  massb = assemble(v * bathymetry2d *dx)
  massl = assemble(v*dx)
  with massl.dat.vec as ml, massb.dat.vec as mb, bathymetry2d.dat.vec as sb:
      ml.reciprocal()
      sb.pointwiseMult(ml, mb)

