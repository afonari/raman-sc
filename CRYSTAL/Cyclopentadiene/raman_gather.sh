#!/bin/bash

backup_dir="Raman_Backup"
mkdir -p ${backup_dir}

rm -f crystal_raman.dat

for dir in ./Modes*
do
    cat "${dir}"/crystal_raman.dat >> crystal_raman.tmp
    cp "${dir}"/OUTCAR.* "${backup_dir}"
done

head -n 1 crystal_raman.tmp > crystal_raman.dat
sed '/^#/d' crystal_raman.tmp > crystal_raman.unsorted
sort -k 2 -n crystal_raman.unsorted >> crystal_raman.dat

rm -f crystal_raman.tmp
rm -f crystal_raman.unsorted
