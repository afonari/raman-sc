#!/bin/bash

gnuplot << EOF
reset
set terminal postscript enhanced "Helvetica" 20
set output "Raman.ps"
set xrange[0:4000]
set style line 1 lt 1 lc -1 lw 4 # black-solid
#set style line 2 lt 2 lc 1 lw 4 # red-dashed
#set style line 3 lt 3 lc 2 lw 4 # red-dashed
# set title "Off-resonance Raman:\nSolid - AO; Dashed - PAW\n@ PBE level."
#set nokey
set xlabel "Energy, cm^{-1}"
set ylabel "Intensity, a.u."

plot "crystal_raman.dat-broaden.dat" u 1:2 w l ls 1
EOF