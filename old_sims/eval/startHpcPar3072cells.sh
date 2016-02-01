#!/usr/bin/env zsh

python genHpcPar.py -c par3072cells -np 3072 -n 5 -horizondir /home/ms926272/repos/horizonParsim -sampledir=/home/ms926272/repos/abstractLTE -t 01:00 -sync cLbtsSync
for x in {0..4}
do
    bsub < hpc_par3072cells$x.sh
    rm hpc_par3072cells$x.sh
done
