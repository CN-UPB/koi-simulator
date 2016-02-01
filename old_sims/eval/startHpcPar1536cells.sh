#!/usr/bin/env zsh

python genHpcPar.py -c par1536cells -np 1536 -n 5 -horizondir /home/ms926272/repos/horizonParsim -sampledir=/home/ms926272/repos/abstractLTE -t 00:30 -sync cLbtsSync
for x in {0..4}
do
    bsub < hpc_par1536cells$x.sh
    rm hpc_par1536cells$x.sh
done
