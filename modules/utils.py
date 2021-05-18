"""
Filename:    utils.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: Helpful generic functions
"""

## Imports

import os, sys

def check_mkdir(filename):
    '''Checks if directory exists and if not, makes that directory'''
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise