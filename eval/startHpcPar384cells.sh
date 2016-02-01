#!/usr/bin/env zsh

python genHpcPar.py -c par384cells -np 384 -n 5 -horizondir /home/ms926272/repos/horizonParsim -sampledir=/home/ms926272/repos/abstractLTE -t 00:15 -sync cLbtsSync
for x in {0..4}
do
    bsub < hpc_par384cells$x.sh
    rm hpc_par384cells$x.sh
done
