#!/usr/bin/bash

set -e

CONFIG=$CONDA_PREFIX/share/psyclone/psyclone.cfg
OPTS="-l output --config ${CONFIG} -s $(pwd)/psyclone/omp_cpu_trans.py"
PREPROC_OPTS='-DARRAYCHECK -DLFRIC'

SOURCE_DIRECTORY='/home/rw796/ukca/r5776_ukca_gpuport_nvtx_calipers/'

# Paths to all Fortran files to be parsed
FILES="
  src/science/core/chemistry/asad/asad_sparse_vars.F90
"

# Need to process the files
for FILE in ${FILES}; do
  out_file=$FILE + '_prep'

  source_path="${SOURCE_DIRECTORY}${FILE}"

  gfortran -E -cpp $PREPROC_OPTS $source_path -o $out_file
  rm -f $FILE
  mv $out_file $FILE

  # Apply PSyclone with no optimisation script
  echo "Applying PSyclone to ${FILE}"
  # psyclone ${OPTS} -o /dev/null ${FILE}
  psyclone ${OPTS} -o "${FILE}_psyclone" ${FILE}
  rm -f $FILE
  mv "${FILE}_psyclone" $FILE

done


