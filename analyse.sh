#!/bin/bash

awk -f mc_mean-sd.awk Montecarlo*/run*.log > mc_mean-sd_2000.csv
awk -f mc_mean-sd.awk Montecarlo*/run{1..500..2}.log > mc_mean-sd_1000.csv
awk -f mc_mean-sd.awk Montecarlo*/run{1..500..4}.log > mc_mean-sd_500.csv
echo 'MC: mean-sd done'

tail -qn 1 Montecarlo*/run*.log > mc_hist.csv
echo 'MC: hist done'

awk -f mc_conv.awk mc_hist.csv > mc_conv.csv
echo 'MC: conv done'

awk -f tr_mean-sd.awk Trapez*/run*.log > tr_mean-sd_200.csv
awk -f tr_mean-sd.awk Trapez*/run{0..200..2}.log > tr_mean-sd_100.csv
awk -f tr_mean-sd.awk Trapez*/run{0..200..4}.log > tr_mean-sd_50.csv
echo 'TR: mean-sd done'

tail -qn 1 Trapez*/run*.log | awk -f tr_conv.awk > tr_conv.csv
echo 'TR: conv done'

pdflatex mc_analyse.tex
pdflatex tr_analyse.tex
rm -f mc_analyse.aux tr_analyse.aux mc_analyse.log tr_analyse.log
