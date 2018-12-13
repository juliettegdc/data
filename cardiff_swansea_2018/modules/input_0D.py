import math
from scipy.interpolate import interp1d
import numpy as np
import os

from modules.parameterisations import turbine_parametrisation


def read_outer_elevations(file):
    """
    Reading elevation time series.
    :param file: input_file with array of format (t, h)
    :return: elevation function
    """

    time_elevation_series = np.load(file)
    return interp1d(time_elevation_series[:,0],time_elevation_series[:,1])

def read_twin_outer_elevations(file):
    """
    Reading elevation time series.
    :param file: input_file with array of format (t, h)
    :return: elevation function
    """
    time_elevation_series = np.load(file)

    HW = interp1d(time_elevation_series[:, 0], time_elevation_series[:, 2])
    LW = interp1d(time_elevation_series[:,0],time_elevation_series[:,1])

    return {"LW": LW, "HW": HW}

def read_area_elevation_curve(file, depth_correction  = 4, area_penalty = 1.0):
    """
    Reading depth-area curve.
    :param file: input_file with array of format (h, area)
    :param depth_correction : amplitude correction in m (depending on datum)
    :return: h-area curve function
    """

    depth_area = np.load(file)
    return interp1d(depth_area[:,0] - depth_correction,depth_area[:,1] * area_penalty)

def read_twin_area_elevation_curve(file, depth_correction  = 4):
    """
    Reading depth-area curve.
    :param file: input_file with array of format (h, area)
    :param depth_correction : amplitude correction in m (depending on datum)
    :return: h-area curve function
    """

    depth_area = np.load(file)
    HW = interp1d(depth_area[:, 0] - depth_correction, depth_area[:, 1])
    LW = interp1d(depth_area[:, 0] - depth_correction, depth_area[:, 2])

    return {"LW":LW, "HW": HW}

def sinusoidal_outer_elevation(amplitude = 5.0, period = 12.42 * 3600):
    """
    Sinusoidal outer elevation based on time in h.
    :param amplitude: amplitude of tidal wave
    :param period: period of tidal wave, assuming M2 period of approximately 12.42 h
    :return:
    """
    omega = (2*math.pi) / period
    elevation = lambda t : amplitude * math.sin(omega * t)
    return elevation

def lagoon_system_idealised_area_case(area = 50):
    """
    Constant area case for twin-basin lagoon
    :param area: total inner lagoon area - assuming they are split evenly.
    :return: h-area curve function
    """
    LW = lambda h : area / 2
    HW = lambda h : area / 2

    return {"LW":LW, "HW": HW}

def lagoon_idealised_area_case(area = 50):
    """
    Constant area case for standard lagoon
    :param area: total inner lagoon area
    :return: h-area curve function
    """
    return lambda h: area