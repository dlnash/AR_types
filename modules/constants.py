"""Collection of convenient functions, constants, and unit labels for atmospheric science

"""

import numpy as np

# Adapter.py when remapping coordinate variables
X_LABELS = ['x', 'lon', 'lons', 'longitude', 'longitudes', 'rlon', 'nlon', 'i', 'nav_lon']
Y_LABELS = ['y', 'lat', 'lats', 'latitude', 'latitudes', 'rlat', 'nlat', 'j', 'nav_lat']
Z_LABELS = ['lev', 'plev', 'level']
T_LABELS = ['time', 'times', 'date', 'dates', 'julian']

# VARIABLES (dicts)

constants = {'earth_angular_velocity': 7.2921e-5,
             'gravity': 9.80665,
            }

units = {'ivt': 'kg m$^{-1}$ s$^{-1}$',
        }

# FUNCTIONS
def coriolis_parameter(latitude):   
    omega = constants['earth_angular_velocity']
    lats_rad = np.deg2rad(latitude)
    f = 2.0 * omega * np.sin(lats_rad)
    
    return f