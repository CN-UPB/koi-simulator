#!/usr/bin/env zsh

python genHpcPar.py -c par192cells -np 192 -n 5 -horizondir /home/ms926272/repos/horizonParsim -sampledir=/home/ms926272/repos/abstractLTE -t 00:05 -sync cLbtsSync
for x in {0..4}
do
    bsub < hpc_par192cells$x.sh
    rm hpc_par192cells$x.sh
done
