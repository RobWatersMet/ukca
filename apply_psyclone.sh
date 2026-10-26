#!/usr/bin/bash

set -e

CONFIG=$CONDA_PREFIX/share/psyclone/psyclone.cfg
OPTS="-l output --config ${CONFIG} -s $(pwd)/psyclone/omp_cpu_trans.py"

# Paths to all Fortran files to be parsed
FILES="
  src/science/core/chemistry/asad/asad_cdrive.F90
  src/science/core/chemistry/asad/asad_diffun.F90
  src/science/core/chemistry/asad/asad_ftoy.F90
  src/science/core/chemistry/asad/asad_fyfixr.F90
  src/science/core/chemistry/asad/asad_fyinit.F90
  src/science/core/chemistry/asad/asad_fyself.F90
  src/science/core/chemistry/asad/asad_inicnt_col_mod.F90
  src/science/core/chemistry/asad/asad_inicnt.F90
  src/science/core/chemistry/asad/asad_prls.F90
  src/science/core/chemistry/asad/asad_sparse_vars.F90
  src/science/core/chemistry/asad/asad_spimpmjp.F90
  src/science/core/chemistry/asad/asad_spmjpdriv.F90
  src/science/core/chemistry/asad/asad_steady.F90
"

# Files that include preprocessing statements
PREPROC="
  src/science/core/chemistry/asad/asad_sparse_vars.F90
"

# Problematic files
EXCLUDED="
"

for FILE in ${FILES}; do

  # Skip problematic files
  tmp=0
  for EXC in ${EXCLUDED}; do
    if [ "${FILE}" == "${EXC}" ]; then
      echo "Skipping ${FILE}"
      tmp=1
    fi
  done
  if [ "${tmp}" == "1" ]; then
    continue
  fi

  # Apply default preprocessing
  for PP in ${PREPROC}; do
    if [ "${FILE}" == "${PP}" ]; then
      fcm revert ${FILE}
      gfortran ${FILE} -cpp -E -C | grep -v "#" >tmp.F90
      mv tmp.F90 ${FILE}
    fi
  done

  # Apply PSyclone with no optimisation script
  echo "Applying PSyclone to ${FILE}"
  # psyclone ${OPTS} -o /dev/null ${FILE}
  psyclone ${OPTS} -o $(pwd)/psyclone/processed/${FILE} ${FILE}
  fcm revert ${FILE}
done
