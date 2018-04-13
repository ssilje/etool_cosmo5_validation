#!/bin/bash
#
# run script for calling the ncl scripts
#   which compare different simulations
#  
# needs: 
# - compare_bias.ncl
# - plot_taylor.ncl (+ taylor_diagram.ncl)
#
# Beate Geyer, HZG, and Susanne Brienen, DWD; 27.1.2015
#---------------------------------------------------------------------
# prepare arguments: variable (loop), experiment IDs

echo 'var=@{var}' 'runids=(/"CON611"/)'  > arguments

# Define and create directories for pdf-plots (outpath) and eps-plots (outpath1)
#outpath="/pool/data/CCLM-EVAL/ANA_CON031/plots/"
#outpath1="/pool/data/CCLM-EVAL/ANA_CON031/plots/epsplots/"

outpath="/scratch/snx3000/ssilje/E-Tool-out/ANA_CON611/plots/"
outpath1="/scratch/snx3000/ssilje/E-Tool-out/ANA_CON611/epsplots/"

[[ -d ${outpath} ]]  || mkdir ${outpath}
[[ -d ${outpath1} ]] || mkdir ${outpath1}

# loop over variables
#for var in "T_2M" "TOT_PREC" "PMSL" "TMAX_2M" "TMIN_2M" "CLCT"; do
for var in "T_2M"; do

  if [[ $var == "CLCT" ]]; then
     OBSDatname=CRU
  else
     OBSDatname=EOBS
  fi
  # compare monthly biases
  ncl `sed s/@{var}/\"${var}\"/g arguments ` obsname=\"${OBSDatname}\" \
       outpath=\"${outpath}\" outpath1=\"${outpath1}\" compare_bias.ncl
  # produce taylor plot of spatial variability
  ncl `sed s/@{var}/\"${var}\"/g arguments ` vartype=\"S\" obsname=\"${OBSDatname}\" \
       outpath=\"${outpath}\" outpath1=\"${outpath1}\" plot_taylor.ncl
  # produce taylor plot of temporal variability
  ncl `sed s/@{var}/\"${var}\"/g arguments ` vartype=\"T\" obsname=\"${OBSDatname}\" \
       outpath=\"${outpath}\" outpath1=\"${outpath1}\" plot_taylor.ncl

done




