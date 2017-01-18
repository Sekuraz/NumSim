#!/bin/bash
REP=4

for METHOD in CG MG RBSOR SOR; do
  echo "  Method ${METHOD}:"
  rm -f drivencavity_${METHOD}.csv

  for N in 8 16 32 64 #128
  do
    echo "Size ${N}:"
    FILE=drivencavity${N}
    rm -f ${FILE}.log
    for i in $(seq 1 ${REP}); do
      echo "  run $i: /usr/bin/time -f\"${N} %e\" -ao ${FILE}.log "\
          "../numsim -s ${METHOD} -I ${FILE} >/dev/null 2>${FILE}_${METHOD}.err"
      /usr/bin/time -f"${N} %e" -ao ${FILE}_${METHOD}.log \
          ../numsim -s ${METHOD} -I ${FILE} >/dev/null 2>${FILE}_${METHOD}.err
    done
    awk -F' ' '{k+=$1; sum+=$2;}; END { print k/NR, sum/NR, NR }' ${FILE}_${METHOD}.log >>drivencavity_${METHOD}.csv
  done
done
