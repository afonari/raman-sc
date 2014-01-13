## Cyclopentadiene example: Split modes

Cyclopentadiene molecule has 11 atoms: C5H6. Thus, there are 33 modes (lowest three acoustic should be close to zero). Note that VASP prints the modes in the *decreasing* eigenvalue (frequency) order. In the current example, the lowest *four* modes are imaginary (probably due to the optimization issues), these will be skipped in the Raman calculation.

The following procedure can be used to obtain Raman activities for the G-point phonons. Change the folder to one containing `POSCAR` and `OUTCAR` from a phonon calculation. `INCAR` should have `NWRITE=3` variable set. The directory should contain:
```
INCAR.raman  - should contain LEPSILON=.TRUE. or LCALCEPS=.TRUE. because we want 'MACROSCOPIC STATIC DIELECTRIC TENSOR' in the OUTCAR
KPOINTS      - just kpoints
OUTCAR.phon  - should contain 'Eigenvectors after division by SQRT(mass)' either from VASP itself or from Henkelman's script
POSCAR.phon  - at this point only positive scales are supported.
POTCAR       - hard PAWs usually give better agreement with the experiment (`ATOMTYPE_h` in the potpaw folders)
raman.sub    - see below
raman_sub.sh - see below
```

Contents of the `raman_sub.sh`:
```bash
#!/bin/bash

DIR=$(pwd)

nderiv_step_size='2_0.01'

export VASP_RAMAN_RUN='mpirun /opt/VASP/bin/vasp5.2.12_p > job.out'

# manually split modes
Modes[0]='01_10'
Modes[1]='11_20'
Modes[2]='21_29'

for m in ${Modes[*]}
do
    VASP_RAMAN_PARAMS="${m}_${nderiv_step_size}"
    export VASP_RAMAN_PARAMS=${VASP_RAMAN_PARAMS}
    #
    FOLDER="Modes_${VASP_RAMAN_PARAMS}"
    #
    echo "Running VASP in ${FOLDER}"
    #
    if [ -d "${FOLDER}" ]; then
        cd ${FOLDER}
        qsub raman.sub
        cd ..
    else
        mkdir ${FOLDER}
        cd ${FOLDER}
        ln -s ../OUTCAR.phon ./OUTCAR.phon
        ln -s ../POSCAR.phon ./POSCAR.phon
        ln -s ../POTCAR ./POTCAR
        ln -s ../INCAR.raman ./INCAR
        ln -s ../raman.sub ./raman.sub
        ln -s ../KPOINTS ./KPOINTS
        qsub raman.sub
        cd ..
    fi
done
```

Contents of the `raman.sub`:
```bash
#!/bin/bash
#PBS -j oe
#PBS -N CPD
#PBS -l nodes=1:ppn=8
#PBS -l pmem=1000mb
#PBS -l walltime=10:00:00
#PBS -V

# environment for the libraries
# [...]

# run the job
#
cd ${PBS_O_WORKDIR}

python27 vasp_raman.py > vasp_raman.out

```

Submit all calculations using: `bash ./raman_sub.sh`. Three folders will be created where calculations will run: `Modes_01_10_2_0.01`, `Modes_11_20_2_0.01` and `Modes_21_29_2_0.01`.

After all calculations are done, it may be a good idea to save OUTCAR files from all the `OUTCAR_*` folders for the future reference. Also, `vasp_raman.dat` will contain Raman activities for the future processing (for example with `broaden.py`).