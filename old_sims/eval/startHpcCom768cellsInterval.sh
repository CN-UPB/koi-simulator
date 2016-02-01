#!/usr/bin/env zsh

python genHpcCom.py -c com768cells -np 64 -n 5 -horizondir /home/ms926272/repos/horizonParsim -sampledir=/home/ms926272/repos/abstractLTE -t 00:15 -sync cIntervalBasedSync -suffix Int
for x in {0..4}
do
    bsub < hpc_com768cellsInt$x.sh
    rm hpc_com768cellsInt$x.sh
done
