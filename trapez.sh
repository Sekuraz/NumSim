#!/bin/bash

N=200

DIR=Trapez$(date +%Y-%m-%d_%H-%M-%S)
mkdir ${DIR}

SIGMA=166.666666666667
MU=1500
START=$(bc <<< "${MU} - 3 * ${SIGMA}")
END=$(bc <<< "${MU} + 3 * ${SIGMA}")
STEP=$(echo "(${END}-${START}) / $N" | bc)

echo "Simulating Re from ${START} to ${END} in steps of ${STEP}"

for n in $(seq 1 $N); do
  RE=$(bc <<< "${START} + $n * ${STEP}")
  ./creator -re ${RE} -tau 0.75 -dt 0.002 -iter 1000
  ./numsim -I drivencavity &>${DIR}/run${n}.log
done

rm -f drivencavity.geom drivencavity.param

