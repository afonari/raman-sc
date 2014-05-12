#!/bin/bash

nderiv_step_size='2_0.01'

export CRYSTAL_RAMAN_RUN='runmpi09 INCAR'

# manually split modes
Modes[0]='03_10'
Modes[1]='11_20'
Modes[2]='21_33'

for m in ${Modes[*]}
do
    CRYSTAL_RAMAN_PARAMS="${m}_${nderiv_step_size}"
    export CRYSTAL_RAMAN_PARAMS=${CRYSTAL_RAMAN_PARAMS}
    #
    FOLDER="Modes_${CRYSTAL_RAMAN_PARAMS}"
    #
    echo "Running CRYSTAL in ${FOLDER}"
    #
    if [ -d "${FOLDER}" ]; then
        cd "${FOLDER}"
        qsub raman.sub
        cd ..
    else
        mkdir "${FOLDER}"
        cd "${FOLDER}"
        ln -s ../INCAR.d12 ./
        ln -s ../FORT34.phon ./
        ln -s ../OUTCAR.phon ./
        ln -s ../raman.sub ./
        qsub raman.sub
        cd ..
    fi
done
