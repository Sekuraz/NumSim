#!/bin/bash

awk -f mc_mean-sd.awk Montecarlo-*/run*.log > mc_mean-sd.csv
awk '(NR-1) % 50 == 0' mc_mean-sd.csv > mc_mean-sd_small.csv

tail -qn 1 Montecarlo-*/run*.log > mc_hist.csv

head -n 500 mc_hist.csv | awk -f mc_conv.awk > mc_conv.csv
head -n 1000 mc_hist.csv | awk -f mc_conv.awk >> mc_conv.csv
awk -f mc_conv.awk mc_hist.csv >> mc_conv.csv

#awk -f tr50_mean-sd.awk Montecarlo-*/run*.log > tr50_mean-sd.csv
#awk '(NR-1) % 50 == 0' tr50_mean-sd.csv > tr50_mean-sd_small.csv

#awk -f tr100_mean-sd.awk Montecarlo-*/run*.log > tr100_mean-sd.csv
#awk '(NR-1) % 50 == 0' tr100_mean-sd.csv > tr100_mean-sd_small.csv

#awk -f tr200_mean-sd.awk Montecarlo-*/run*.log > tr200_mean-sd.csv
#awk '(NR-1) % 50 == 0' tr200_mean-sd.csv > tr200_mean-sd_small.csv

#awk -f tr_conv.awk Montecarlo-*/run*.log > tr_conv.csv
#tail -qn 1 Montecarlo-*/run*.log > tr_hist.csv

pdflatex mc_analyse.tex
#pdflatex tr50_analyse.tex
#pdflatex tr100_analyse.tex
#pdflatex tr200_analyse.tex
rm -f mc_analyse.aux #tr*_analyse.aux mc_analyse.log tr*_analyse.log
