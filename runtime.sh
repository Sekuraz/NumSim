#!/bin/bash

###################
#
#  runtime.sh
#
#  A script to run the current numsim executable with 1 to 14 processes for runtime measurements.
#  It uses the default driven cavity and writes output to time.log and s128x128n${proccount}
#
###################

for N in $(seq 14)
do
    /usr/bin/time -a -f "N = ${N}: \t%E real,\t%U user,\t%S sys" -o time.log mpirun -n ${N} ./numsim &>s128x128n${N}.log
done
