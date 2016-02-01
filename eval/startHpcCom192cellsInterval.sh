#!/usr/bin/env zsh

python genHpcCom.py -c com192cells -np 16 -n 5 -horizondir /home/ms926272/repos/horizonParsim -sampledir=/home/ms926272/repos/abstractLTE -t 00:05 -sync cIntervalBasedSync -suffix Int
for x in {0..4}
do
    bsub < hpc_com192cellsInt$x.sh
    rm hpc_com192cellsInt$x.sh
done
