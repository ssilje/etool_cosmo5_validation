#!/bin/csh

##PBS -l nodes=1:bwgrid:ppn=1,walltime=04:00:00
##PBS -N ETOOL
##PBS -o ETOOL.out
##PBS -e ETOOL.err

###############################################################################
#####################      ETOOL v2.3      ####################################
#####################      2014-03-26      ####################################
###############################################################################
# This program is designed to perform some basic analysis of COSMO-CLM output #
# data and validate to the E-OBS gridded observational data.                  #
# Monthly statistics for temperature (mean, min, max) precipitation and       #
# sea-level-pressure are calculated, and the standard setup validates the     #
# simulations for the PRUDENCE regions of Europe.                             #
# It is also possible to define custom regions.                               #
#                                                                             #
# The program works for all horizontal resolutions, as the validation is      #
# calculated only for area means of the selected sub-regions. Note that only  #
# land points (FR_LAND >= 0.5) are included in the analysis.                  #
# Program can be called by                                                    #
#     ./ETOOL_v2.3.sh <variable> <run-ID>                                     #
#   or using                                                                  #
#     ETOOL_submit.csh                                                        #
###############################################################################
# OUTPUT:                                                                     #
# The output files are named according to the CDO commands used in the        #
# analysis. Files starting with a user defined run id are CLM files, while    #
# E-OBS data starts with "EOBS". The monthly data are saved in a file         #
# containing the whole period, and then the period mean is calculated. The    #
# sub-regions are then extracted ("REGX").                                    #
# The difference to the E-OBS data is then calculated and stored in the files #
# named "<RUNID>_DIFF...".                                                    #
# Note that one can define different periods for the observations and the     #
# model data, but the output files will then carry only information about     #
# the model time-period.                                                      #
###############################################################################
# DEPENDENCIES:                                                               #
# Climate Data Operators (CDO)                                                #
#     only standard tools are used so the program should not be sensitive to  #
#     the exact version used.                                                 #
#     Download at [http://www.mpimet.mpg.de/fileadmin/software/cdo/]          #
# E-OBS data (regular grid) :                                                 #
# http://eca.knmi.nl/download/ensembles/ensembles.php                         #
# Cloud fraction data from CRU: CRU_TS3.22                                    #
#  http://www.cru.uea.ac.uk/cru/data/hrg/                                     #
###############################################################################
# CHANGES:                                                                    #
# 13/08/2014, Susanne Brienen, DWD: extension to TMIN_2M, TMAX_2M, PMSL       #
# 18/08/2014, Klaus Keuler, BTU: code optimization                            #
# 30/09/2014, Susanne Brienen, DWD: refomulation for call per variable,       #
#              restart option (only for crash during reading loop over years!)#
# 05/12/2014, Andrew Ferrone, CRP-GL: sub-directories per variable for        #
#              parallel usage                                                 #
# 30/01/2015, Andrew Ferrone, LIST: introduce options to decreace verbosity   #
#              and introduce check for number of files after untaring         #
# 11/02/2015, Susanne Brienen, DWD: run tar only if data are not yet in the   #
#              working directory; do one more mergetime step in case of hourly#
#              model output; introduce new variable: total cloud cover in     #
#              comparison with CRU_TS3.22 data                                #
# 10/03/2015  Klaus Keuler, BTU: manipulation of dates for correct calculation#
#              of monthly means by using CDO shifttime,-1minute command;      #
#              correction of path allocations in working directories;         #
#              optimization of CDO commands by concatenating operations       #       
# 24/03/2015  Susanne Brienen, DWD: save daily values as well; include option #
#              to use subchain-post-output (monthly time series) instead of   #
#              direct model output; produce seasonal means for Taylor plots;  #
#              small corrections and optimizations                            #
###############################################################################
# AUTHOR: Peter Berg, IMK-TRO, KIT (berg@kit.edu)                             #
###############################################################################

###############################################################################
##### USER INPUT SECTION ######################################################
###############################################################################
# VarName1         - variable name (notation as in CCLM output)               #
# RUNID            - Any identifier for the analyzed simulation               #
# StartYearM/O     - First year to be analysed of the model/observations      #
# EndYearM/O       - Last year to be analysed of the model/observations       #
# EOBSRES          - Resolution of the E-OBS data (0.25 or 0.50 degrees)      #
# EOBSVER          - Version of the E-OBS data                                #
# RefDataT_2M      - E-OBS file for temperature (original version)            #
# RefDataTOT_PREC  - E-OBS file for precipitation (original version)          #
# RefDataTMIN_2M   - E-OBS file for minimum temperature (original version)    #
# RefDataTMAX_2M   - E-OBS file for maximum temperature (original version)    #
# RefDataPMSL      - E-OBS file for sea-level pressure (original version)     #
# CloudPath        - Path for cloud fraction data                             #
# RefDataCLCT      - OBS file for cloud fraction (original version)           #
# CLMPath          - Path to the base of the CLM output directories           #
# CLMLSMFile       - Path to a file containing the FR_LAND parameter          #
# CLMPath_T_2M_AV  - Path to the temperature variable output                  #
# CLMPath_TOT_PREC - Path to the precipitation variable output                #
# CLMPath_T_MIN_2M - Path to the minimum temperature variable output          #
# CLMPath_T_MAX_2M - Path to the maximum temperature variable output          #
# CLMPath_PMSL     - Path to the sea-level pressure variable output           #
# CLMPath_CLCT     - Path to the cloud fraction variable output               #
# RX               - Whether region X should be analysed (T) or not (F).      #
#                    The regions are: R1-BI, R2-IP, R3-FR, R4-ME, R5-SC,      # 
#                    R6-AL, R7-MD, R8-EA, R9-(extra), R0-(extra2)             #
# SAVpath          - Directory for storage of output                          #
# TMPDIR           - Temporary directory (Note: this will be deleted!)        #
# SAVpathT         - Directory for storage of output for Taylor plots         #
# SUBCHAIN_TS      - T for using post-processed time series of subchain       #
#                    scripts instead of direct model output                   #
# HINC_*           - time increment in hours of different variables in CCLM   #
#                    output                                                   #
###############################################################################

# Clean ${TMP_DIR} even if interupted #########################################

#trap ' echo "Caught Signal ... cleaning up." 1>&2 ' 1 2 3 6
#       rm -r ${TMPDIR} ;
#       rm -r ${workspace}/out03 ;
#       rm -r ${workspace}/out04 ;
#       echo "Done cleanup ... quitting." 1>&2 ;
#       exit 1 ' 1 2 3 6

###############################################################################

setenv VarName1 $1
echo 'selected variable: ' $VarName1

set REST=F  # if "T": restart after last completed year
set RUNID=`echo $2 | sed 's/[^a-zA-Z0-9_-]//g'`  #CON516
echo 'for simulation: ' $RUNID   
set StartYearM=1981 
set EndYearM=2000 
eval set StartYearO=${StartYearM} 
eval set EndYearO=${EndYearM}
set EOBSRES=0.25
set EOBSVER=10.0 
setenv workspace ${SCRATCH}/E-Tool/${RUNID}
setenv savespace ${SCRATCH}/E-Tool-out
setenv EOBSPath /project/pr04/observations/eobs_${EOBSRES}deg_reg_v${EOBSVER}
setenv RefDataT_2M tg_${EOBSRES}deg_reg_v${EOBSVER}.nc
setenv RefDataTOT_PREC rr_${EOBSRES}deg_reg_v${EOBSVER}.nc
setenv RefDataTMIN_2M tn_${EOBSRES}deg_reg_v${EOBSVER}.nc
setenv RefDataTMAX_2M tx_${EOBSRES}deg_reg_v${EOBSVER}.nc
setenv RefDataPMSL pp_${EOBSRES}deg_reg_v${EOBSVER}.nc
setenv RefDataElev elev_${EOBSRES}deg_reg_v${EOBSVER}.nc
setenv CloudPath /project/pr04/observations/clct
setenv RefDataCLCT cru_ts3.22.1901.2013.cld.dat.nc
set CLMPath  = /project/pr04/cosmo5_validation/${RUNID}/output
set CLMLSMFile=/project/pr04/cosmo5_validation/lffd1979010100c.nc
setenv CLMInDir_T_2M out04
setenv CLMInDir_TOT_PREC out03
setenv CLMInDir_TMIN_2M out07
setenv CLMInDir_TMAX_2M out07
setenv CLMInDir_PMSL out04
setenv CLMInDir_CLCT out04
set R1=T
set R2=T
set R3=T
set R4=T
set R5=T
set R6=T
set R7=T
set R8=T
set R9=T  # If "T" you will also need to specify the coordinates:
setenv REGION_9 -31.0,61.0,34.0,76.0 
set R0=F  # If "T" you will also need to specify the coordinates:
setenv REGION_0 6.0,15.0,47.2,54.9
set SAVpath=${savespace}/ANA_${RUNID}
if (! -e ${SAVpath}) mkdir -p ${SAVpath}
if ( ${REST} == F ) then
##  if (-e ${workspace}) rm -rf ${workspace}
  mkdir -p ${workspace}
else
  if (! -e ${workspace}) mkdir -p ${workspace}
endif
set TMPDIR=`mktemp -dp ${workspace}`
set SAVpathT=${SAVpath}/taylor
if (! -e ${SAVpathT}) mkdir -p ${SAVpathT}
# Set verbosity of output
# Recommended settings for testing (i.e. more verbose)
#setenv tar_opt -xvf
#setenv cdo_opt 
# Recommended setting for running (i.e. less verbose)
setenv tar_opt -xf
setenv cdo_opt -s

# if you want to use pre-processed time series of subchain-post:
set SUBCHAIN_TS=F

# some variable-specific settings which could be different for other simulations
# time increment of output in hours
set HINC_T_2M=3
set HINC_TOT_PREC=1
set HINC_TMIN_2M=24
set HINC_TMAX_2M=24
set HINC_PMSL=3
set HINC_CLCT=3

# in case of problems, see also comments further down starting with "note:"
#  when e.g. your model output data is packed in a different structure 

###############################################################################
## END OF USER INPUT ##########################################################
############################################################################### 

############################################################################### 
## more settings which normally are not necessary to adapt ####################
############################################################################### 

switch ( $VarName1 )
   case T_2M:
             setenv VarName1proc daymean
             setenv VarName1proc2 monmean
             setenv VarName1proc2alt FALSE
             breaksw
   case T???_2M:
             setenv VarName1proc daymean
             setenv VarName1proc2 monmean
             setenv VarName1proc2alt FALSE
             breaksw
   case TOT_PREC:
             setenv VarName1proc daysum
             setenv VarName1proc2 monmean
             setenv VarName1proc2alt monsum
             breaksw
   case PMSL:
             setenv VarName1proc daymean
             setenv VarName1proc2 monmean
             setenv VarName1proc2alt FALSE
             breaksw
   case CLCT:
   # data are already monthly means, so ${VarName1proc*} not really necessary here!
   #  but keep for consistent file names
             setenv VarName1proc daymean
             setenv VarName1proc2 monmean
             setenv VarName1proc2alt FALSE
             breaksw
   default:
      echo 'variable ' $VarName1 ' not yet considered in ETOOL!! Quit program!'
      exit 1
      breaksw
endsw
echo 'selected operators: '
echo 'VarName1proc ' $VarName1proc
echo 'VarName1proc2 ' $VarName1proc2
echo 'VarName1proc2alt ' $VarName1proc2alt

# specific settings per variable
eval set CLMInDir_var=\${CLMInDir_${VarName1}}
eval set HINC=\${HINC_${VarName1}}
eval set RefData=\${RefData${VarName1}}

set CLMPath_var=${workspace}/${CLMInDir_var}

if ( $VarName1 == 'CLCT' ) then
   set OBSPath=${CloudPath}  # input data path
   set OBSDatname='CRU'      # prefix for output file name
else
   set OBSPath=${EOBSPath}
   set OBSDatname='EOBS'
endif
echo 'obs path: ' $OBSPath

# calculate expected number of file names
# note: this should be changed if incomplete years are processed!!
@ exp_nb_files = (24 / $HINC) * 365    #+ 1
echo 'minimum number of expected files: ' $exp_nb_files

# Variable name in OBS data
switch ( $VarName1 )
   case T_2M*:
             set EOBSVAR='tg'
             breaksw
   case TOT_PREC:
             set EOBSVAR='rr'
             breaksw
   case TMIN_2M:
             set EOBSVAR='tn'
             breaksw
   case TMAX_2M:
             set EOBSVAR='tx'
             breaksw
   case PMSL:
             set EOBSVAR='psl'
             breaksw
   case CLCT:
             set EOBSVAR='cld'
             breaksw
   default:
      echo 'ERROR: variable not available in OBS data! Quit program'
      exit 1
      breaksw
endsw


# This environment variable makes sure that the rotated grid is used by the CDOs
setenv IGNORE_ATT_COORDINATES 0

# BI - British Isles
setenv REGION_1 -10.0,2.0,50.0,59.0
# IP - Iberian Peninsula
setenv REGION_2 -10.0,3.0,36.0,44.0
# FR - France
setenv REGION_3 -5.0,5.0,44.0,50.0
# ME - Middle Europe
setenv REGION_4 2.0,16.0,48.0,55.0
# SC - Scandinavia
setenv REGION_5 5.0,30.0,55.0,70.0
# AL - Alps
setenv REGION_6 5.0,15.0,44.0,48.0
# MD - Mediterranean
setenv REGION_7 3.0,25.0,36.0,44.0
# EA - Eastern Europe
setenv REGION_8 16.0,30.0,44.0,55.0

# determine nuber of years in specified time period
@ ny = ($EndYearM - $StartYearM) + 1
echo 'number of years: ' $ny 

# cdo processing for Taylor plot calculation
setenv VarProc1 yseasmean      # climatological seasonal processing cdo command
setenv VarProc2 timmean        # climatological annual processing cdo command
setenv VarProc3 seasmean       # seasonal processing cdo command
setenv VarProc4 yearmean       # annual processing cdo command

setenv SEASONS 'DJF MAM JJA SON year'
##############################################################################
## PERFORM THE COSMO-CLM CALCULATIONS ########################################
##############################################################################
echo `date`

# Extract the land-sea mask and remap to regular grid #
cdo ${cdo_opt} remapbil,${EOBSPath}/${RefDataElev} -selvar,FR_LAND ${CLMLSMFile} ${TMPDIR}/clmlsm.nc 
cdo ${cdo_opt} setrtomiss,0,.5 ${TMPDIR}/clmlsm.nc ${SAVpath}/clmlsm.nc
# for restart:
if ( ${REST} == T ) then
   # get last year of previous run
   cd ${workspace}
   set Year0="`cat year.dat`"
   @ Year = ${Year0} + 1
else
   set Year=$StartYearM
endif
echo 'start at year: ' $Year
while (${Year} <= ${EndYearM})

 echo 'Processing year ' ${Year}

 @ YearN = ${Year} + 1

 if ( ${SUBCHAIN_TS} == T ) then  # use time series produced by subchain-post

   cdo ${cdo_opt} mergetime ${CLMPath}/${Year}_??/${VarName1}_ts.nc ${TMPDIR}/${VarName1}_${Year}.nc

   # Remove first date of current year
   cdo ${cdo_opt} delete,timestep=1 ${TMPDIR}/${VarName1}_${Year}.nc tmp.nc
   mv tmp.nc ${TMPDIR}/${VarName1}_${Year}.nc

 else # use direct model output

   # Untar files from $CLMPath to $workspace, but only if not yet available
   if ( -e ${CLMPath_var} ) then
    cd ${CLMPath_var}

    # Remove first date of current year (Year)
    if ( -f lffd${Year}010100.nc ) then
      rm -f lffd${Year}010100.nc 
      echo "remove file " lffd${Year}010100.nc 
    endif

    set nb_files=`find -L . -name "lffd${Year}??????.nc" | wc -l`

    # Check if first date of next year (YearN) is available
    if ( -f lffd${YearN}010100.nc ) then
      @ nb_files = ${nb_files} + 1
      echo "last date of year " ${Year} " available: " lffd${YearN}010100.nc 
    else # maybe it is enough just to copy again this one file
      # note: you might need to adapt this tar-command if your data is tared in a different structure!
      #       In the following, the individual data files for each time step per year are 
      #       expected to be in a subfolder of ${workspace} with the name ${CLMPath_var}
      tar ${tar_opt} ${CLMPath}/year${Year}/cclm_${CLMInDir_var}_????12????_??????????.tar \
          ${CLMInDir_var}/lffd${YearN}010100.nc --strip-components=1
      if ( -f lffd${YearN}010100.nc ) then
        @ nb_files = ${nb_files} + 1
        echo "last date of year " ${Year} " now available: " lffd${YearN}010100.nc 
      endif
    endif
   else
    set nb_files = 0
    echo "No files of year " ${Year} " available"
   endif 

   if ( ${nb_files} < ${exp_nb_files} ) then
    echo 'untar CLM data to $workspace...'
    cd ${workspace}
    # note: you might need to adapt this loop if your data is tared in a different structure!
    #       In the following, the individual data files for each time step per year are 
    #       expected to be in a subfolder of ${workspace} with the name ${CLMPath_var}
    foreach fil (${CLMPath}/${CLMInDir_var}/${RUNID}_${CLMInDir_var}_${Year}.tar)
      tar ${tar_opt} ${fil}
    end
    mv ${Year} ${CLMInDir_var}
    cd ${CLMInDir_var}
    if (-l lffd${YearN}010100.nc) rm -f lffd${YearN}010100.nc
    # check if all files are there
    # Remove first date of current year (Year)
    if ( -f lffd${Year}010100.nc ) then
      rm -f lffd${Year}010100.nc 
      echo "remove file " lffd${Year}010100.nc 
    endif
    set nb_files=`find -L . -type f -name "lffd${Year}??????.nc" | wc -l`
    # Check if first date of next year (YearN) is available
    if ( -f lffd${YearN}010100.nc ) then
      @ nb_files = ${nb_files} + 1
      echo "last date of year " ${Year} " available: " lffd${YearN}010100.nc 
    else # maybe it's in the package of next year's January:
      # note: you might need to adapt this loop if your data is tared in a different structure!
      #       In the following, the individual data files for each time step per year are 
      #       expected to be in a subfolder of ${workspace} with the name ${CLMPath_var}
      tar ${tar_opt} ${CLMPath}/${CLMInDir_var}/${RUNID}_${CLMInDir_var}_${YearN}.tar \
          ${YearN}/lffd${YearN}010100.nc --strip-components=1
      if ( -f lffd${YearN}010100.nc ) then
        @ nb_files = ${nb_files} + 1
        echo "last date of year " ${Year} " now available: " lffd${YearN}010100.nc
      else
        echo "still not there..."
      endif
    endif
    if  ( ${nb_files} < ${exp_nb_files} ) then
      echo 'ERROR: Missing files in: ' ${CLMPath_var}} '  Quit program'
      exit 1
    endif
    echo '...done'
   else
    echo 'CLM data already available in ' $workspace
   endif

   # Extract the variables of interest #
   cd ${CLMPath_var}
   echo 'extract the variables of interest...'

   foreach fil (lffd${Year}??????.nc lffd${YearN}010100.nc)
     cdo ${cdo_opt} selvar,${VarName1} $fil ${TMPDIR}/${fil}_${VarName1}
     #echo $fil
   end
   echo '...done'

   # Merge time (this is performed in two steps as the nr. of files can be too large for the OS) #
   echo 'merge time...'
   #
   foreach imon (01 02 03 04 05 06 07 08 09 10 11 12)
    # even one more step for hourly data
    if ( ${HINC} == 1 ) then
      # calculate number of days in month
      set ndays=`cal ${imon} ${Year} | egrep -v '[a-z]' | wc -w`
      foreach iday ( `seq -w ${ndays}` )
        cdo ${cdo_opt} mergetime ${TMPDIR}/lffd${Year}${imon}${iday}??.nc_${VarName1} \
                      ${TMPDIR}/${VarName1}_${Year}_${imon}_${iday}.nc
      end
      cdo ${cdo_opt} mergetime ${TMPDIR}/${VarName1}_${Year}_${imon}_??.nc ${TMPDIR}/${VarName1}_${Year}_${imon}.nc
    else
      cdo ${cdo_opt} mergetime ${TMPDIR}/lffd${Year}${imon}????.nc_${VarName1} ${TMPDIR}/${VarName1}_${Year}_${imon}.nc
    endif
   end
   cdo ${cdo_opt} mergetime ${TMPDIR}/${VarName1}_${Year}_??.nc ${TMPDIR}/lffd${YearN}010100.nc_${VarName1} \
                 ${TMPDIR}/${VarName1}_${Year}.nc
   echo '...done'

 endif # subchain-post output or direct model output


 # some special conversions to the model data #
 switch ( $VarName1 )
    case T*2M*:
             echo '(nothing special to do at the moment for model temperature data)'
             breaksw
    case TOT_PREC:
             echo '(nothing special to do at the moment for model precipitation data)'
             breaksw
    case PMSL:
             # for PMSL: convert model data from Pa to hPa (EOBS data are in hPa)
             echo 'convert model mean sea level pressure data from Pa to hPa'
             cdo ${cdo_opt} divc,100.00 ${TMPDIR}/${VarName1}_${Year}.nc ${TMPDIR}/${VarName1}_${Year}_hPa.nc
             # change also meta data
             cdo ${cdo_opt} setunit,'hPa' ${TMPDIR}/${VarName1}_${Year}_hPa.nc ${TMPDIR}/${VarName1}_${Year}.nc
             rm ${TMPDIR}/${VarName1}_${Year}_hPa.nc
             breaksw
    case CLCT:
             echo '(nothing special to do at the moment for model cloud cover data)'
             breaksw
    default:
      echo 'ERROR: Quit program'
      exit 1
      breaksw
 endsw

 echo 'perform remapping and aggregation of data...'
 # Remap to observational grid in regular coordinates
 cdo ${cdo_opt} remapbil,${EOBSPath}/${RefDataElev} \
     ${TMPDIR}/${VarName1}_${Year}.nc ${TMPDIR}/${VarName1}_${Year}_remapped.nc

 # Set sea points to NaN #
 cdo ${cdo_opt} add ${TMPDIR}/${VarName1}_${Year}_remapped.nc ${SAVpath}/clmlsm.nc ${TMPDIR}/${VarName1}_${Year}.nc2
 cdo ${cdo_opt} sub ${TMPDIR}/${VarName1}_${Year}.nc2 ${SAVpath}/clmlsm.nc ${TMPDIR}/${VarName1}_${Year}.nc

 # Shift time by 1 minute backward to assign the correct calendar date to every timestep at end of day
 cdo ${cdo_opt} -shifttime,-1minute ${TMPDIR}/${VarName1}_${Year}.nc   ${TMPDIR}/${VarName1}_${Year}_shift.nc

 # to make sure that only current year is retained, check again number of time steps
 cdo ${cdo_opt} selyear,${Year} ${TMPDIR}/${VarName1}_${Year}_shift.nc ${TMPDIR}/${VarName1}_${Year}_shift2.nc 
 set nb_files_2=`cdo ${cdo_opt} ntime ${TMPDIR}/${VarName1}_${Year}_shift2.nc`
 if  ( ${nb_files_2} < ${exp_nb_files} ) then
   if ( ${SUBCHAIN_TS} == T ) then  # use time series produced by subchain-post
      echo 'ERROR: Missing files in: ' ${workspace} '; expected number: ' ${exp_nb_files} ' -- Quit program'
   else
      echo 'ERROR: Missing files in: ' ${CLMPath_var} '; expected number: ' ${exp_nb_files} '  Quit program'
   endif
      exit 1
 endif

 # Sum or average data #
 cdo ${cdo_opt} ${VarName1proc} ${TMPDIR}/${VarName1}_${Year}_shift2.nc ${TMPDIR}/${VarName1}_${VarName1proc}_${Year}.nc
 cdo ${cdo_opt} ${VarName1proc2} ${TMPDIR}/${VarName1}_${VarName1proc}_${Year}.nc \
     ${TMPDIR}/${VarName1}_${VarName1proc}_${VarName1proc2}_${Year}.nc
 if (${VarName1proc2alt} != FALSE) cdo ${cdo_opt} ${VarName1proc2alt} ${TMPDIR}/${VarName1}_${VarName1proc}_${Year}.nc \
             ${TMPDIR}/${VarName1}_${VarName1proc}_${VarName1proc2alt}_${Year}.nc

 if (${Year} == ${StartYearM}) then
   mv ${TMPDIR}/${VarName1}_${VarName1proc}_${VarName1proc2}_${Year}.nc ${SAVpath}/${VarName1}_1_merge.tmpOLD
   if (${VarName1proc2alt} != FALSE) mv ${TMPDIR}/${VarName1}_${VarName1proc}_${VarName1proc2alt}_${Year}.nc \
             ${SAVpath}/${VarName1}_2_merge.tmpOLD
   # save also daily values for whole domain
   mv ${TMPDIR}/${VarName1}_${VarName1proc}_${Year}.nc ${SAVpath}/${VarName1}_daily_merge.tmpOLD

 endif
 if (${Year} != ${StartYearM}) then
   cdo ${cdo_opt} mergetime ${SAVpath}/${VarName1}_1_merge.tmpOLD \
       ${TMPDIR}/${VarName1}_${VarName1proc}_${VarName1proc2}_${Year}.nc ${SAVpath}/${VarName1}_1_merge.tmp
   mv ${SAVpath}/${VarName1}_1_merge.tmp ${SAVpath}/${VarName1}_1_merge.tmpOLD
   if (${VarName1proc2alt} != FALSE) then
         cdo ${cdo_opt} mergetime ${SAVpath}/${VarName1}_2_merge.tmpOLD \
                       ${TMPDIR}/${VarName1}_${VarName1proc}_${VarName1proc2alt}_${Year}.nc \
                       ${SAVpath}/${VarName1}_2_merge.tmp
         mv ${SAVpath}/${VarName1}_2_merge.tmp ${SAVpath}/${VarName1}_2_merge.tmpOLD
   endif
   cdo ${cdo_opt} mergetime ${SAVpath}/${VarName1}_daily_merge.tmpOLD \
       ${TMPDIR}/${VarName1}_${VarName1proc}_${Year}.nc ${SAVpath}/${VarName1}_daily_merge.tmp
   mv ${SAVpath}/${VarName1}_daily_merge.tmp ${SAVpath}/${VarName1}_daily_merge.tmpOLD
 endif
 echo '...done'

 cd ${workspace}
 # write last completed year to file for restart #
 echo $Year >! year.dat
 # Clear up the tmpdir in between #
 rm -r ${TMPDIR}
 if ( ${SUBCHAIN_TS} == F ) then
   rm -r ${CLMInDir_var}
 endif
 mkdir ${TMPDIR}
 @ Year = ${Year} + 1

end # year

echo 'calculate statistics of multi-annual data...'
# Rename the merged files #
mv ${SAVpath}/${VarName1}_1_merge.tmpOLD \
   ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearM}-${EndYearM}.nc
if (${VarName1proc2alt} != FALSE)  mv ${SAVpath}/${VarName1}_2_merge.tmpOLD \
          ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearM}-${EndYearM}.nc
mv ${SAVpath}/${VarName1}_daily_merge.tmpOLD ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${StartYearM}-${EndYearM}.nc

# Calculate statistics #
cdo ${cdo_opt} ymonmean ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearM}-${EndYearM}.nc \
    ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearM}-${EndYearM}_ymonmean.nc
if (${VarName1proc2alt} != FALSE) cdo ${cdo_opt} ymonmean \
    ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearM}-${EndYearM}.nc \
    ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearM}-${EndYearM}_ymonmean.nc
echo '...done'

##############################################################################
## PREPARE THE OBSERVATIONAL DATA ############################################
##############################################################################
echo 'prepare OBS data...'

if ( ${VarName1} == 'CLCT' ) then   # CRU data

  # select relevant period and convert to fractions
  cdo ${cdo_opt} -b 4 -setunit,'1' -mulc,0.01 -selyear,${StartYearO}/${EndYearO} ${OBSPath}/${RefData} \
       ${TMPDIR}/${OBSDatname}_${VarName1}_${StartYearO}-${EndYearO}_frac.nc

  # remap to EOBS grid
  cdo ${cdo_opt} remapbil,${EOBSPath}/${RefDataElev} \
        ${TMPDIR}/${OBSDatname}_${VarName1}_${StartYearO}-${EndYearO}_frac.nc \
        ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearO}-${EndYearO}.nc

  # data are already monthly, so ${VarName1proc2} not necessary here, go directly to ymonmean
  cdo ${cdo_opt} ymonmean ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearO}-${EndYearO}.nc \
        ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearO}-${EndYearO}_ymonmean.nc

  echo '*NOTE* no daily observational data available for ' $VarName1 '; File ' \
    ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${StartYearO}-${EndYearO}.nc ' will not be written'

else   # EOBS data

  cdo ${cdo_opt} -b 4 selyear,${StartYearO}/${EndYearO} \
        ${EOBSPath}/${RefData} ${TMPDIR}/${OBSDatname}_${VarName1}_${StartYearO}-${EndYearO}_K.nc

  # special manipulations depending on variable
  switch ( $VarName1 )
    case T*2M*:
      # T_2M, TMIN_2M, TMAX_2M #
      # convert to Kelvin
      cdo ${cdo_opt} -setunit,'K' -addc,273.15 ${TMPDIR}/${OBSDatname}_${VarName1}_${StartYearO}-${EndYearO}_K.nc \
           ${TMPDIR}/${OBSDatname}_${VarName1}_${StartYearO}-${EndYearO}.nc
      breaksw
    case TOT_PREC:
      # TOT_PREC #
      mv ${TMPDIR}/${OBSDatname}_${VarName1}_${StartYearO}-${EndYearO}_K.nc \
         ${TMPDIR}/${OBSDatname}_${VarName1}_${StartYearO}-${EndYearO}.nc 
      echo '(nothing special to do at the moment for precipitation)'
      breaksw
    case PMSL:
      # PMSL #
      mv ${TMPDIR}/${OBSDatname}_${VarName1}_${StartYearO}-${EndYearO}_K.nc \
         ${TMPDIR}/${OBSDatname}_${VarName1}_${StartYearO}-${EndYearO}.nc 
      echo '(nothing special to do at the moment for mean sea level pressure)'
      breaksw
    default:
      echo 'variable ' $VarName1 ' not yet considered in ETOOL!! Quit program!'
      exit 1
      breaksw
  endsw

  cdo ${cdo_opt} ${VarName1proc2} ${TMPDIR}/${OBSDatname}_${VarName1}_${StartYearO}-${EndYearO}.nc \
      ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearO}-${EndYearO}.nc
  cdo ${cdo_opt} ymonmean \
      ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearO}-${EndYearO}.nc \
      ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearO}-${EndYearO}_ymonmean.nc

  if (${VarName1proc2alt} != FALSE) then
    cdo ${cdo_opt} ${VarName1proc2alt} ${TMPDIR}/${OBSDatname}_${VarName1}_${StartYearO}-${EndYearO}.nc \
        ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearO}-${EndYearO}.nc
    cdo ${cdo_opt} ymonmean \
        ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearO}-${EndYearO}.nc \
        ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearO}-${EndYearO}_ymonmean.nc
  endif

  # save also daily data
  cp ${TMPDIR}/${OBSDatname}_${VarName1}_${StartYearO}-${EndYearO}.nc \
     ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${StartYearO}-${EndYearO}.nc

endif

echo '...done'

##############################################################################
## LOOP OVER THE SUB-REGIONS AND CALCULATE THE DIFFERENCES ###################
##############################################################################
echo 'Calculating statistics for the sub-regions'
foreach reg (1 2 3 4 5 6 7 8 9 0)
 set rloop=1
 eval if \( \$R${reg} != T \) set rloop=0
 if ( $rloop == 1) then

  # OBS #
  eval cdo ${cdo_opt} sellonlatbox,\$REGION_${reg} \
           ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearO}-${EndYearO}_ymonmean.nc \
           ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearO}-${EndYearO}_ymonmean_REG${reg}.nc
  if (${VarName1proc2alt} != FALSE) eval cdo ${cdo_opt} sellonlatbox,\$REGION_${reg} \
           ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearO}-${EndYearO}_ymonmean.nc \
           ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearO}-${EndYearO}_ymonmean_REG${reg}.nc
  # for monthly data produce field means over regions
  eval cdo ${cdo_opt} sellonlatbox,\$REGION_${reg} \
           ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearO}-${EndYearO}.nc \
           ${TMPDIR}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearO}-${EndYearO}_REG${reg}.nc
  cdo ${cdo_opt} fldmean ${TMPDIR}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearO}-${EndYearO}_REG${reg}.nc \
           ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearO}-${EndYearO}_REG${reg}mn.nc
  if (${VarName1proc2alt} != FALSE) then
     eval cdo ${cdo_opt} sellonlatbox,\$REGION_${reg} \
           ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearO}-${EndYearO}.nc \
           ${TMPDIR}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearO}-${EndYearO}_REG${reg}.nc
     cdo ${cdo_opt} fldmean ${TMPDIR}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearO}-${EndYearO}_REG${reg}.nc \
           ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearO}-${EndYearO}_REG${reg}mn.nc
  endif
  # cut out regions also for daily data
  if (${VarName1} != "CLCT") then
     eval cdo ${cdo_opt} sellonlatbox,\$REGION_${reg} \
           ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${StartYearO}-${EndYearO}.nc \
           ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${StartYearO}-${EndYearO}_REG${reg}.nc
  endif
  # CLM #
  eval cdo ${cdo_opt} sellonlatbox,\$REGION_${reg} \
           ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearM}-${EndYearM}_ymonmean.nc \
           ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearM}-${EndYearM}_ymonmean_REG${reg}.nc
  if (${VarName1proc2alt} != FALSE) eval cdo ${cdo_opt} sellonlatbox,\$REGION_${reg} \
           ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearM}-${EndYearM}_ymonmean.nc \
           ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearM}-${EndYearM}_ymonmean_REG${reg}.nc
  # for monthly data produce field means over regions
  eval cdo ${cdo_opt} sellonlatbox,\$REGION_${reg} \
           ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearM}-${EndYearM}.nc \
           ${TMPDIR}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearM}-${EndYearM}_REG${reg}.nc
  cdo ${cdo_opt} fldmean \
           ${TMPDIR}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearM}-${EndYearM}_REG${reg}.nc \
           ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearM}-${EndYearM}_REG${reg}mn.nc
  if (${VarName1proc2alt} != FALSE) then
     eval cdo ${cdo_opt} sellonlatbox,\$REGION_${reg} \
           ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearM}-${EndYearM}.nc \
           ${TMPDIR}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearM}-${EndYearM}_REG${reg}.nc
     cdo ${cdo_opt} fldmean \
           ${TMPDIR}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearM}-${EndYearM}_REG${reg}.nc \
           ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearM}-${EndYearM}_REG${reg}mn.nc
  endif
  # cut out regions also for daily data
  eval cdo ${cdo_opt} sellonlatbox,\$REGION_${reg} \
           ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${StartYearM}-${EndYearM}.nc \
           ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${StartYearM}-${EndYearM}_REG${reg}.nc

  #Calculate differences #
  cdo ${cdo_opt} sub \
    ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearM}-${EndYearM}_ymonmean_REG${reg}.nc \
    ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearO}-${EndYearO}_ymonmean_REG${reg}.nc \
    ${SAVpath}/${RUNID}_DIFF_CLM-${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearM}-${EndYearM}_ymonmean_REG${reg}.nc
  if (${VarName1proc2alt} != FALSE) cdo ${cdo_opt} sub \
    ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearM}-${EndYearM}_ymonmean_REG${reg}.nc \
    ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearO}-${EndYearO}_ymonmean_REG${reg}.nc \
    ${SAVpath}/${RUNID}_DIFF_CLM-${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearM}-${EndYearM}_ymonmean_REG${reg}.nc

 endif
end

echo '... done'

##############################################################################
## SPECIAL PREPARATIONS FOR PLOTTING TAYLOR DIAGRAMS #########################
##############################################################################
echo 'Calculating correlations and variances for Taylor diagrams'
foreach reg (1 2 3 4 5 6 7 8 9 0)
 set rloop=1
 eval if \( \$R${reg} != T \) set rloop=0
 if ( $rloop == 1) then

   # prepare seasonal and annual mean data
   cdo ${cdo_opt} ${VarProc1} \
       ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearM}-${EndYearM}_ymonmean_REG${reg}.nc \
       ${TMPDIR}/sim_temp_seas.nc
   cdo ${cdo_opt} ${VarProc2} \
       ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearM}-${EndYearM}_ymonmean_REG${reg}.nc \
       ${TMPDIR}/sim_temp_year.nc
   cdo ${cdo_opt} ${VarProc1} \
       ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearO}-${EndYearO}_ymonmean_REG${reg}.nc \
       ${TMPDIR}/ref_temp_seas.nc
   cdo ${cdo_opt} ${VarProc2} \
       ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearO}-${EndYearO}_ymonmean_REG${reg}.nc \
       ${TMPDIR}/ref_temp_year.nc
   if (${VarName1proc2alt} != FALSE) then
     cdo ${cdo_opt} ${VarProc1} \
         ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearM}-${EndYearM}_ymonmean_REG${reg}.nc \
         ${TMPDIR}/sim_temp_seas_2.nc
     cdo ${cdo_opt} ${VarProc2} \
         ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearM}-${EndYearM}_ymonmean_REG${reg}.nc \
         ${TMPDIR}/sim_temp_year_2.nc
     cdo ${cdo_opt} ${VarProc1} \
         ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearO}-${EndYearO}_ymonmean_REG${reg}.nc \
         ${TMPDIR}/ref_temp_seas_2.nc
     cdo ${cdo_opt} ${VarProc2} \
         ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearO}-${EndYearO}_ymonmean_REG${reg}.nc \
         ${TMPDIR}/ref_temp_year_2.nc
   endif

   # calculate spatial correlation and variances of simulation and reference data
   # and the normalized variance for each season and the annual mean 
   cdo ${cdo_opt} fldcor ${TMPDIR}/sim_temp_seas.nc ${TMPDIR}/ref_temp_seas.nc \
       ${SAVpathT}/${RUNID}_${VarName1}_${VarName1proc2}_SCOR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_seas_REG${reg}.nc
   cdo ${cdo_opt} fldcor ${TMPDIR}/sim_temp_year.nc ${TMPDIR}/ref_temp_year.nc \
       ${SAVpathT}/${RUNID}_${VarName1}_${VarName1proc2}_SCOR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_year_REG${reg}.nc
   cdo ${cdo_opt} -div -fldvar ${TMPDIR}/sim_temp_seas.nc -fldvar ${TMPDIR}/ref_temp_seas.nc  \
       ${SAVpathT}/${RUNID}_${VarName1}_${VarName1proc2}_SNVAR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_seas_REG${reg}.nc
   cdo ${cdo_opt} -div -fldvar ${TMPDIR}/sim_temp_year.nc -fldvar ${TMPDIR}/ref_temp_year.nc  \
       ${SAVpathT}/${RUNID}_${VarName1}_${VarName1proc2}_SNVAR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_year_REG${reg}.nc
   if (${VarName1proc2alt} != FALSE) then
     cdo ${cdo_opt} fldcor ${TMPDIR}/sim_temp_seas_2.nc ${TMPDIR}/ref_temp_seas_2.nc \
         ${SAVpathT}/${RUNID}_${VarName1}_${VarName1proc2alt}_SCOR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_seas_REG${reg}.nc
     cdo ${cdo_opt} fldcor ${TMPDIR}/sim_temp_year_2.nc ${TMPDIR}/ref_temp_year_2.nc \
         ${SAVpathT}/${RUNID}_${VarName1}_${VarName1proc2alt}_SCOR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_year_REG${reg}.nc
     cdo ${cdo_opt} -div -fldvar ${TMPDIR}/sim_temp_seas_2.nc -fldvar ${TMPDIR}/ref_temp_seas_2.nc  \
         ${SAVpathT}/${RUNID}_${VarName1}_${VarName1proc2alt}_SNVAR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_seas_REG${reg}.nc
     cdo ${cdo_opt} -div -fldvar ${TMPDIR}/sim_temp_year_2.nc -fldvar ${TMPDIR}/ref_temp_year_2.nc  \
         ${SAVpathT}/${RUNID}_${VarName1}_${VarName1proc2alt}_SNVAR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_year_REG${reg}.nc
   endif

   rm -f ${TMPDIR}/???_temp_*.nc

   # cut and calculate area mean of aregion
   # calculate (and split) seasonal and yearly means of simulation and reference data
   eval cdo ${cdo_opt} -splitseas -${VarProc3} -fldmean -sellonlatbox,\$REGION_${reg} \
            ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearM}-${EndYearM}.nc \
            ${TMPDIR}/sim_temp_
   eval cdo ${cdo_opt} -splitseas -${VarProc3} -fldmean -sellonlatbox,\$REGION_${reg} \
            ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearO}-${EndYearO}.nc \
            ${TMPDIR}/ref_temp_
   eval cdo ${cdo_opt} -${VarProc4} -fldmean -sellonlatbox,\$REGION_${reg} \
            ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearM}-${EndYearM}.nc \
            ${TMPDIR}/sim_temp_year.nc
   eval cdo ${cdo_opt} -${VarProc4} -fldmean -sellonlatbox,\$REGION_${reg} \
            ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2}_${StartYearO}-${EndYearO}.nc \
            ${TMPDIR}/ref_temp_year.nc
   if (${VarName1proc2alt} != FALSE) then
     eval cdo ${cdo_opt} -splitseas -${VarProc3} -fldmean -sellonlatbox,\$REGION_${reg} \
              ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearM}-${EndYearM}.nc \
              ${TMPDIR}/sim_temp_2_
     eval cdo ${cdo_opt} -splitseas -${VarProc3} -fldmean -sellonlatbox,\$REGION_${reg} \
              ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearO}-${EndYearO}.nc \
              ${TMPDIR}/ref_temp_2_
     eval cdo ${cdo_opt} -${VarProc4} -fldmean -sellonlatbox,\$REGION_${reg} \
              ${SAVpath}/${RUNID}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearM}-${EndYearM}.nc \
              ${TMPDIR}/sim_temp_2_year.nc
     eval cdo ${cdo_opt} -${VarProc4} -fldmean -sellonlatbox,\$REGION_${reg} \
              ${SAVpath}/${OBSDatname}_${VarName1}_${VarName1proc}_${VarName1proc2alt}_${StartYearO}-${EndYearO}.nc \
              ${TMPDIR}/ref_temp_2_year.nc
   endif

   # calculate temporal correlation and variances of simulation and reference data
   # and the normalized variance for each season and the annual mean
   foreach seas (${SEASONS})
      #   exclude last december value from winter season time series
      if ( ${seas} == "DJF" ) then
        cdo ${cdo_opt} splitsel,$ny ${TMPDIR}/sim_temp_DJF.nc ${TMPDIR}/simDJF_
        mv ${TMPDIR}/simDJF_000000.nc ${TMPDIR}/sim_temp_DJF.nc
        cdo ${cdo_opt} splitsel,$ny ${TMPDIR}/ref_temp_DJF.nc ${TMPDIR}/refDJF_
        mv ${TMPDIR}/refDJF_000000.nc ${TMPDIR}/ref_temp_DJF.nc
        if (${VarName1proc2alt} != FALSE) then
          cdo ${cdo_opt} splitsel,$ny ${TMPDIR}/sim_temp_2_DJF.nc ${TMPDIR}/simDJF_2_
          mv ${TMPDIR}/simDJF_2_000000.nc ${TMPDIR}/sim_temp_2_DJF.nc
          cdo ${cdo_opt} splitsel,$ny ${TMPDIR}/ref_temp_2_DJF.nc ${TMPDIR}/refDJF_2_
          mv ${TMPDIR}/refDJF_2_000000.nc ${TMPDIR}/ref_temp_2_DJF.nc
        endif
      endif

      cdo ${cdo_opt} timcor ${TMPDIR}/sim_temp_${seas}.nc ${TMPDIR}/ref_temp_${seas}.nc \
          ${TMPDIR}/${RUNID}_${VarName1}_${VarName1proc2}_TCOR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_${seas}_REG${reg}.nc
      cdo ${cdo_opt} -div -timvar ${TMPDIR}/sim_temp_${seas}.nc -timvar ${TMPDIR}/ref_temp_${seas}.nc  \
          ${TMPDIR}/${RUNID}_${VarName1}_${VarName1proc2}_TNVAR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_${seas}_REG${reg}.nc
      if (${VarName1proc2alt} != FALSE) then
        cdo ${cdo_opt} timcor ${TMPDIR}/sim_temp_2_${seas}.nc ${TMPDIR}/ref_temp_2_${seas}.nc \
            ${TMPDIR}/${RUNID}_${VarName1}_${VarName1proc2alt}_TCOR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_${seas}_REG${reg}.nc
        cdo ${cdo_opt} -div -timvar ${TMPDIR}/sim_temp_2_${seas}.nc -timvar ${TMPDIR}/ref_temp_2_${seas}.nc  \
            ${TMPDIR}/${RUNID}_${VarName1}_${VarName1proc2alt}_TNVAR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_${seas}_REG${reg}.nc
      endif

   end

   # merge seasonal values into one file
   cdo ${cdo_opt} -O mergetime \
       ${TMPDIR}/${RUNID}_${VarName1}_${VarName1proc2}_TCOR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_???_REG${reg}.nc \
       ${SAVpathT}/${RUNID}_${VarName1}_${VarName1proc2}_TCOR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_seas_REG${reg}.nc
   cdo ${cdo_opt} -O mergetime \
       ${TMPDIR}/${RUNID}_${VarName1}_${VarName1proc2}_TNVAR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_???_REG${reg}.nc \
       ${SAVpathT}/${RUNID}_${VarName1}_${VarName1proc2}_TNVAR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_seas_REG${reg}.nc
   if (${VarName1proc2alt} != FALSE) then
     cdo ${cdo_opt} -O mergetime \
         ${TMPDIR}/${RUNID}_${VarName1}_${VarName1proc2alt}_TCOR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_???_REG${reg}.nc \
         ${SAVpathT}/${RUNID}_${VarName1}_${VarName1proc2alt}_TCOR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_seas_REG${reg}.nc
     cdo ${cdo_opt} -O mergetime \
         ${TMPDIR}/${RUNID}_${VarName1}_${VarName1proc2alt}_TNVAR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_???_REG${reg}.nc \
         ${SAVpathT}/${RUNID}_${VarName1}_${VarName1proc2alt}_TNVAR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_seas_REG${reg}.nc
   endif

   mv ${TMPDIR}/${RUNID}_${VarName1}_${VarName1proc2}_TCOR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_year_REG${reg}.nc ${SAVpathT}
   mv ${TMPDIR}/${RUNID}_${VarName1}_${VarName1proc2}_TNVAR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_year_REG${reg}.nc ${SAVpathT}
   if (${VarName1proc2alt} != FALSE) then
     mv ${TMPDIR}/${RUNID}_${VarName1}_${VarName1proc2alt}_TCOR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_year_REG${reg}.nc ${SAVpathT}
     mv ${TMPDIR}/${RUNID}_${VarName1}_${VarName1proc2alt}_TNVAR_CLM-${OBSDatname}_${StartYearM}-${EndYearM}_year_REG${reg}.nc ${SAVpathT}
   endif

   rm -f ${TMPDIR}/???_temp_*.nc

 endif
end

echo '... done'


rm -r ${TMPDIR}

echo 'Calculations finished at:',`date`

