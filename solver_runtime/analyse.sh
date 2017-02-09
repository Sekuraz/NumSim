#!/bin/bash

###################
#
#  analyse.sh
#
#  A script to measure the runtime of the current source tree comparing the different solvers for the Driven-Cavity.
#  It is using ${REP} repetition to get more reliable results. This is the sequential version.
#  The output is written to runtime_seq_${date}/${METHOD}.csv
#
###################

REP=5
DIR=runtime_seq_$(date +%Y-%m-%d_%H-%M-%S)
mkdir ${DIR}

rm -f drivencavity*.log drivencavity*.err

cmake -H.. -B. -DWRITE_VTK=off -DDEBUG_VISU=off -DCMAKE_BUILD_TYPE=Release
make -j

for METHOD in CG MG RBSOR SOR; do
  echo "Method ${METHOD}:"

  for N in 16 32 64 128 256 #8 # not used because too small
  do
    echo "  Size ${N}:"
    FILE=drivencavity${N}
    for i in $(seq 1 ${REP}); do
      echo "  run $i: /usr/bin/time -f\"${N} %e\" -ao ${FILE}.log "\
          "./numsim -s ${METHOD} -I ${FILE} >/dev/null 2>${FILE}_${METHOD}_${i}.err"
      /usr/bin/time -f"${N} %e" -ao ${FILE}_${METHOD}.log \
          ./numsim -s ${METHOD} -I ${FILE} >/dev/null 2>${FILE}_${METHOD}_${i}.err
    done
    awk -F' ' '{k+=$1; sum+=$2;}; END { print k/NR, sum/NR, NR }' ${FILE}_${METHOD}.log >>${DIR}/${METHOD}.csv
  done
done
