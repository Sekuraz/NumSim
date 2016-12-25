#!/bin/bash

DIR=Montecarlo$(date +%Y-%m-%d_%H-%M-%S)

mkdir ${DIR}

for I in $(seq 1 500); do
  ./creator -rand 1500x166.666666666667 -tau 0.75 -dt 0.002 -iter 1000
  ./numsim -I drivencavity &>${DIR}/run${I}.log
done

rm -f drivencavity.geom drivencavity.param

