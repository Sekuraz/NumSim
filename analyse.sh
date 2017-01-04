#!/bin/bash

awk -f mc_mean-sd.awk Montecarlo*/run*.log > mc_mean-sd.csv
awk '(NR-1) % 50 == 0' mc_mean-sd.csv > mc_mean-sd_small.csv
echo 'MC: mean-sd done'

tail -qn 1 Montecarlo*/run*.log > mc_hist.csv
echo 'MC: hist done'

awk -f mc_conv.awk mc_hist.csv > mc_conv.csv
echo 'MC: conv done'

awk -f tr_mean-sd.awk Trapez*/run{0..200..4}.log > tr50_mean-sd.csv
awk '(NR-1) % 50 == 0' tr50_mean-sd.csv > tr50_mean-sd_small.csv
awk -f tr_mean-sd.awk Trapez*/run{0..200..2}.log > tr100_mean-sd.csv
awk '(NR-1) % 50 == 0' tr100_mean-sd.csv > tr100_mean-sd_small.csv
awk -f tr_mean-sd.awk Trapez*/run*.log > tr200_mean-sd.csv
awk '(NR-1) % 50 == 0' tr200_mean-sd.csv > tr200_mean-sd_small.csv
echo 'TR: mean-sd done'

tail -qn 1 Trapez*/run*.log | awk -f tr_conv.awk > tr_conv.csv
echo 'TR: conv done'

pdflatex mc_analyse.tex
pdflatex tr50_analyse.tex
pdflatex tr100_analyse.tex
pdflatex tr200_analyse.tex
rm -f mc_analyse.aux tr*_analyse.aux mc_analyse.log tr*_analyse.log
