#!/bin/bash

for METHOD in CG MG RBSOR SOR; do
  if [ -z "$@" ]; then
    PARAM="-s ${METHOD} -I profiler"
  else
    PARAM="-s ${METHOD} $@"
  fi
  valgrind --tool=callgrind --callgrind-out-file=profiler.out ./numsim ${PARAM}
  ./gprof2dot.py --format=callgrind profiler.out -o profiler
  dot -Tsvg -o profiler_${METHOD}$(date  +%Y-%m-%d_%H-%M-%S).svg profiler
done
rm -f profiler{,.out}
