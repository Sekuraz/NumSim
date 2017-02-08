#!/bin/bash

###################
#
#  profiler.sh
#
#  A script to get the call graph for the different solvers.
#  It accepts the same arguments as the numsim executable, except for the solver which is set here.
#  The output files are profiler_${METHOD}.svg
#
###################

for METHOD in CG MG RBSOR SOR; do
  if [ -z "$@" ]; then
    PARAM="-s ${METHOD} -I profiler"
  else
    PARAM="-s ${METHOD} $@"
  fi
  valgrind --tool=callgrind --callgrind-out-file=profiler.out ./numsim ${PARAM}
  ./gprof2dot.py --format=callgrind profiler.out -o profiler
  dot -Tsvg -o profiler_${METHOD}.svg profiler
done
rm -f profiler{,.out}
