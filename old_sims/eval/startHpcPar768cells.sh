#!/usr/bin/env zsh

python genHpcPar.py -c par768cells -np 768 -n 5 -horizondir /home/ms926272/repos/horizonParsim -sampledir=/home/ms926272/repos/abstractLTE -t 00:15 -sync cLbtsSync
for x in {0..4}
do
    bsub < hpc_par768cells$x.sh
    rm hpc_par768cells$x.sh
done
