#!/bin/bash
######################################################################
# Filename:    copy_figs.sh
# Author:      Deanna Nash dlnash@ucsb.edu
# Description: Script to copy final figures for manuscript to one folder
#
######################################################################

# Input parameters
ceofid="IVT"
ssn="djfmam"
neofs=2
nk=3
maindir="/home/sbarc/students/nash/repositories/AR_types/figs/" # main figure folder
finaldir="/home/sbarc/students/nash/repositories/AR_types/figs/final_figs/" # final figure folder

# fig names in main folder
array=(
fig1_elevation_one-minute_with_inset
ivt_clim_all_ssn
fig3_prec_clim 
prec_clim
prec_clim_frac
prec_std
composite_climate_indices
)

array2=(
1
2
3
S1
S2
S3
S7
)

for i in ${!array[*]}
do 
    infile="${maindir}${array[$i]}.png"
    outfile="${finaldir}fig${array2[$i]}.png"
#     echo "${infile} to ${outfile}"
    cp -v ${infile} ${outfile}
done

# fig names in mainfol + ceofid + ssn
maindir="${maindir}${ceofid}/${ssn}/"
infile="${maindir}spatial_tcor.png"
outfile="${finaldir}fig5.png"
# echo "${infile} to ${outfile}"
cp -v ${infile} ${outfile}

# fig names in mainfol + ceofid + ssn + neof + nk
maindir="${maindir}/neof${neofs}/k${nk}/"
array=(
ARTYPE_freq_barplot
composite_9pan
lag0_composite_anomalies_9pan
composite_diff_ENSO
composite_diff_AO
lag2_composite_anomalies_9pan
lag-2_composite_anomalies_9pan
composite_diff_MJO
fig4_variance_elbow_kde 
fig5_transition_heatmap 
boxplot
composite_diff_SH
)

array2=(
S4
7
8
9
10
S5
S6
S9
4
6
S8
11
)

for i in ${!array[*]}
do 
    infile="${maindir}${array[$i]}.png"
    outfile="${finaldir}fig${array2[$i]}.png"
#     echo "${infile} to ${outfile}"
    cp -v ${infile} ${outfile}
done