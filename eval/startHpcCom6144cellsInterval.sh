#!/usr/bin/env zsh

python genHpcCom.py -c com6144cells -np 512 -n 5 -horizondir /home/ms926272/repos/horizonParsim -sampledir=/home/ms926272/repos/abstractLTE -t 00:30 -sync cIntervalBasedSync -suffix Int
for x in {0..4}
do
    bsub < hpc_com6144cellsInt$x.sh
    rm hpc_com6144cellsInt$x.sh
done
