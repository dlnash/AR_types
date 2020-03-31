"""Collection of convenient functions, constants, and unit labels for atmospheric science

"""

import numpy as np


# VARIABLES (dicts)

constants = {'earth_angular_velocity': 7.292e-5,
             'gravity': 9.80665,
            }


units = {'ivt': 'kg m$^{-1}$ s$^{-1}$',
        }



# FUNCTIONS

def coriolis_parameter(latitude):   
    
    lats_rad = np.deg2rad(latitude)
    f = 2.0 * omega * np.sin(lats_rad)
    
    return f