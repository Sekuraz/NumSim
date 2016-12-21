#!/bin/bash

valgrind --tool=callgrind --callgrind-out-file=profiler.out ./numsim $@
./gprof2dot.py --format=callgrind profiler.out -o profiler
dot -Tsvg -O profiler
rm -f profiler
