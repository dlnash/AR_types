"""
Filename:    getERA5_sfc_batch_ivt.py
Author:      Tessa Montini, tmontini@ucsb.edu & Deanna Nash, dlnash@ucsb.edu
Description: Download multi-year ERA5 data on single or pressure levels. Use in conjunction with ERA5_config for input variables. 

"""

import cdsapi
import yaml


### IMPORTANT! CHANGE CONFIG NAME FOR DOWNLOAD DICTIONARY ###
config_name = 'ivt'


# import configuration file for season dictionary choice
yaml_doc = 'ERA5_config.yml'
config = yaml.load(open(yaml_doc), Loader=yaml.SafeLoader)
ddict = config[config_name]


# Loop for downloading annual data files
for yr in range(ddict['start_yr'],ddict['end_yr']+1):
    outfile = ddict['datadir'] + "{0}_{1}.nc".format(ddict['fprefix'], yr)
    c = cdsapi.Client()
    c.retrieve(ddict['data_type'], 
               {'product_type'  : 'reanalysis',
                'variable'      : ddict['var_name'],
                'pressure_level': ddict['levels'],
                'year'          : "{0}".format(yr),
                'month'         : ddict['month'],
                'day'           : ddict['year'],
                'time'          : ddict['time'],
                'area'          : ddict['area'],
                'grid'          : ddict['grid'],
                'format'        : 'netcdf'}, 
               outfile)
    print("Download complete: {filename} \n".format(filename=outfile))