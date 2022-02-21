#!/bin/bash
######################################################################
# Filename:    nco_era5_6hr_to_daily_mean.sh
# Author:      Deanna Nash dlnash@ucsb.edu
# Description: Script to get daily mean from 6-hourly ERA5
#
# - Convert 6-hourly ERA5 data to daily
# - To run: change input parameters, then go to NCO_NRCAT directory
# - "conda activate nco-env", then run "bash nco_era5_6hr_to_daily_mean.sh"
######################################################################

# Input parameters
vars='iwv'
start_yr=1979
end_yr=2019
region='hma'
resol='05'
datadir="/home/sbarc/students/nash/data/ERA5/${vars}/6hr/"
outdir="/home/sbarc/students/nash/data/ERA5/${vars}/daily/"

outer=1      # set outer loop counter

# Loop to concatenate all files within month-year into 1 netcdf file
# Begin outer loop (e.g. each year)
cd $datadir
for year in $(seq $start_yr $end_yr)
do
    echo "Pass $outer in outer loop."
    inner=1    # reset inner loop counter
    infile="era5_${region}_${resol}dg_6hr_${vars}_${year}.nc"
    echo "$infile"
    outfile="${outdir}out.era5_${region}_${resol}dg_daily_${vars}_${year}.nc"
    echo "$outfile"
    # calculate daily mean from 6hourly
    cdo daymean ${infile} ${outfile}

    let "outer+=1" # Increment outer loop counter
    echo "$year concatenation complete"
    echo           # Space between output blocks in pass of outer loop

done
# End of outer loop

exit 0