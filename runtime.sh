#!/bin/bash

for N in $(seq 14)
do
/usr/bin/time -a -f "N = ${N}: \t%E real,\t%U user,\t%S sys" -o time.log mpirun -n ${N} ./numsim &>s128x128n${N}.log
done
