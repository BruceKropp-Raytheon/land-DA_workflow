#!/bin/sh

set -xue

# Set other dates
PTIME=$($NDATE -${DATE_CYCLE_FREQ_HR} $PDY$cyc)

YYYY=${PDY:0:4}
MM=${PDY:4:2}
DD=${PDY:6:2}
HH=${cyc}
YYYP=${PTIME:0:4}
MP=${PTIME:4:2}
DP=${PTIME:6:2}
HP=${PTIME:8:2}

FILEDATE=${YYYY}${MM}${DD}.${HH}0000

JEDI_STATICDIR=${JEDI_PATH}/jedi-bundle/fv3-jedi/test/Data
JEDI_EXECDIR=${JEDI_PATH}/build/bin

case $MACHINE in
  "hera")
    run_cmd="srun"
    ;;
  "orion")
    run_cmd="srun"
    ;;
  "hercules")
    run_cmd="srun"
    ;;
  *)
    run_cmd=`which mpiexec`
    ;;
esac

# copy sfc_data files into work directory
for itile in {1..6}
do
  sfc_fn="${FILEDATE}.sfc_data.tile${itile}.nc"
  if [ -f ${DATA_RESTART}/${sfc_fn} ]; then
    cp -p ${DATA_RESTART}/${sfc_fn} .
  elif [ -f ${WARMSTART_DIR}/${sfc_fn} ]; then
    cp -p ${WARMSTART_DIR}/${sfc_fn} .
  else
    err_exit "Initial sfc_data files do not exist"
  fi
  # copy sfc_data file for comparison
  cp -p ${sfc_fn} "${sfc_fn}_ini"
done
# Copy obserbation file to work directory
mkdir -p ${DATA}/obs
obs_type_lower="${OBS_TYPE,,}"
obs_suffix="${obs_type_lower}_snow.nc"
ln -nsf "${COMIN}/obs/${obs_type_lower}_snow_${PDY}${cyc}.nc" "${DATA}/obs/obs.${cycle}.${obs_suffix}"

# update coupler.res file
settings="\
  'coupler_calendar': ${COUPLER_CALENDAR}
  'yyyp': !!str ${YYYP}
  'mp': !!str ${MP}
  'dp': !!str ${DP}
  'hp': !!str ${HP}
  'yyyy': !!str ${YYYY}
  'mm': !!str ${MM}
  'dd': !!str ${DD}
  'hh': !!str ${HH}
" # End of settings variable

fp_template="${PARMlandda}/templates/template.coupler.res"
fn_namelist="${DATA}/${FILEDATE}.coupler.res"
${USHlandda}/fill_jinja_template.py -u "${settings}" -t "${fp_template}" -o "${fn_namelist}"

# Copy static data files
mkdir -p ${DATA}/Data/fv3files
mkdir -p ${DATA}/Data/fieldmetadata
cp -p ${PARMlandda}/jedi/fv3files/fmsmpp.nml ${DATA}/Data/fv3files/.
cp -p ${JEDI_STATICDIR}/fv3files/field_table_ufs ${DATA}/Data/fv3files/field_table
ln -nsf ${JEDI_STATICDIR}/fv3files/akbk${NPZ}.nc4 ${DATA}/Data/fv3files/akbk.nc4

if [ "${FRAC_GRID}" = "YES" ]; then
  snowdepth_vn="snodl"
  cp -p ${JEDI_STATICDIR}/fieldmetadata/gfs_v17-land.yaml ${DATA}/Data/fieldmetadata/gfs-land.yaml
else
  snowdepth_vn="snwdph"
  cp -p ${PARMlandda}/jedi/fieldmetadata/gfs-land.yaml ${DATA}/Data/fieldmetadata/gfs-land.yaml
fi

# Intermediate/Output directories
mkdir -p ${DATA}/diags
mkdir -p ${DATA}/anl
mkdir -p ${DATA}/bkg
mkdir -p ${DATA}/test

# Set up background
if [ "${JEDI_ALGORITHM}" = "3dvar" ]; then
  for itile in {1..6}
  do
    sfc_fn="${FILEDATE}.sfc_data.tile${itile}.nc"
    sfc_bkg_fn="${FILEDATE}.sfc_data.tile${itile}.nc"
    cp -p ${sfc_fn} "${DATA}/bkg/${sfc_bkg_fn}"
    ln -nsf "${FIXlandda}/FV3_fix_tiled/C${RES}/C${RES}_oro_data.tile${itile}.nc" "${DATA}/bkg/."
  done
  cp -p ${FILEDATE}.coupler.res ${DATA}/bkg

  # Set JEDI executable
  jedi_exe_fn="fv3jedi_var.x"

else # letkf
  for ens in {1..2}
  do
    mkdir -p $DATA/mem${ens}
    cp -p ${FILEDATE}.sfc_data.tile*.nc ${DATA}/mem${ens}
    cp -p ${FILEDATE}.coupler.res ${DATA}/mem${ens}
  done

  ${USHlandda}/letkf_create_ens.py $FILEDATE $snowdepth_vn 30
  if [[ $? != 0 ]]; then
    err_exit "letkf create failed"
  fi

  # Set JEDI executable
  jedi_exe_fn="fv3jedi_letkf.x"
fi

# Copy JEDI input yaml file
jedi_nml_fn="jedi_${JEDI_ALGORITHM}_snow.yaml"
cp -p "${COMIN}/${jedi_nml_fn}" .

################################################
# RUN JEDI
################################################

export pgm="${jedi_exe_fn}"
. prep_step
${run_cmd} -n ${NPROCS_ANALYSIS} ${JEDI_EXECDIR}/$pgm ${jedi_nml_fn} >>$pgmout 2>errfile
export err=$?; err_chk
cp errfile errfile_jedi_letkf
if [[ $err != 0 ]]; then
  err_exit "JEDI DA failed"
fi

# save intermediate sfc_data files
for itile in {1..6}
do
  sfc_fn="${FILEDATE}.sfc_data.tile${itile}.nc"
  cp -p ${sfc_fn} "${sfc_fn}_old"
done

################################################
# Apply Increment to UFS sfc_data files
################################################

# temporary for apply_incr
for itile in {1..6}
do
  ln -nsf ${FILEDATE}.snowinc.sfc_data.tile${itile}.nc ${FILEDATE}.xainc.sfc_data.tile${itile}.nc
done

if [ "${FRAC_GRID}" = "YES" ]; then
  frac_grid=".true."
else
  frac_grid=".false."
fi
orog_path="${FIXlandda}/FV3_fix_tiled/C${RES}"
orog_fn_base="C${RES}_oro_data"

cat << EOF > apply_incr_nml
&noahmp_snow
 date_str=${YYYY}${MM}${DD}
 hour_str=${HH}
 res=${RES}
 frac_grid="${frac_grid}"
 rst_path="${DATA}"
 inc_path="${DATA}"
 orog_path="${orog_path}"
 otype=${orog_fn_base}
 ntiles=6
 ens_size=1
/
EOF

export pgm="apply_incr.exe"
. prep_step
# (n=6): this is fixed, at one task per tile (with minor code change). 
${run_cmd} -n 6 ${EXEClandda}/$pgm >>$pgmout 2>errfile
export err=$?; err_chk
cp errfile errfile_apply_incr
if [[ $err != 0 ]]; then
  err_exit "apply snow increment failed"
fi

for itile in {1..6}
do
  cp -p ${DATA}/${FILEDATE}.snowinc.sfc_data.tile${itile}.nc ${COMOUT}
done 

for itile in {1..6}
do
  cp -p ${DATA}/${FILEDATE}.sfc_data.tile${itile}.nc ${COMOUT}
done

if [ -d diags ]; then
  cp -p diags/* ${COMOUThofx}
  ln -nsf ${COMOUThofx}/* ${DATA_HOFX}
fi


############################################################
# Comparison plot of sfc_data by JEDI increment
############################################################
DO_PLOT_SFC_COMP="YES"
if [ "${DO_PLOT_SFC_COMP}" = "YES" ]; then

  fn_sfc_base="${FILEDATE}.sfc_data.tile"
  fn_inc_base="${FILEDATE}.snowinc.sfc_data.tile"
  out_title_base="Land-DA::SFC-DATA::${PDY}::"
  out_fn_base="landda_comp_sfc_${PDY}_"
  # zlevel_number is valid only for 3-D fields such as stc/smc/slc
  zlevel_number="1"

  cat > plot_comp_sfc.yaml <<EOF
work_dir: '${DATA}'
fn_sfc_base: '${fn_sfc_base}'
fn_inc_base: '${fn_inc_base}'
jedi_exe: '${JEDI_ALGORITHM}'
orog_path: '${orog_path}'
orog_fn_base: '${orog_fn_base}'
out_title_base: '${out_title_base}'
out_fn_base: '${out_fn_base}'
fix_dir: '${FIXlandda}'
zlevel_number: '${zlevel_number}'
EOF

  ${USHlandda}/plot_comp_sfc_data.py
  if [ $? -ne 0 ]; then
    err_exit "sfc_data comparison plot failed"
  fi

  # Copy result file to COMOUT
  cp -p ${out_fn_base}* ${COMOUTplot}

fi


###########################################################
# WE2E test
###########################################################
if [ "${WE2E_TEST}" == "YES" ]; then
  path_fbase="${FIXlandda}/test_base/we2e_com/${RUN}.${PDY}"
  fn_sfc="${FILEDATE}.sfc_data.tile"
  fn_inc="${FILEDATE}.snowinc.sfc_data.tile"
  fn_hofx="diag.ghcn_snow_${PDY}${cyc}.nc"
  we2e_log_fp="${LOGDIR}/${WE2E_LOG_FN}"
  if [ ! -f "${we2e_log_fp}" ]; then
    touch ${we2e_log_fp}
  fi
  # surface data tiles
  for itile in {1..6}
  do
    ${USHlandda}/compare.py "${path_fbase}/${fn_sfc}${itile}.nc" "${COMOUT}/${fn_sfc}${itile}.nc" ${WE2E_ATOL} ${we2e_log_fp} "ANALYSIS" ${FILEDATE} "sfc_data.tile${itile}"
  done
  # increment tiles
  for itile in {1..6}
  do
    ${USHlandda}/compare.py "${path_fbase}/${fn_inc}${itile}.nc" "${COMOUT}/${fn_inc}${itile}.nc" ${WE2E_ATOL} ${we2e_log_fp} "ANALYSIS" ${FILEDATE} "snowinc.tile${itile}"
  done
  # H(x)
  ${USHlandda}/compare.py "${path_fbase}/hofx/${fn_hofx}" "${COMOUT}/hofx/${fn_hofx}" ${WE2E_ATOL} ${we2e_log_fp} "ANALYSIS" ${FILEDATE} "HofX"
fi
