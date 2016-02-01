#!/usr/bin/env zsh

python genHpcCom.py -c com384cells -np 32 -n 5 -horizondir /home/ms926272/repos/horizonParsim -sampledir=/home/ms926272/repos/abstractLTE -t 00:10 -sync cLbtsSync
for x in {0..4}
do
    bsub < hpc_com384cells$x.sh
    rm hpc_com384cells$x.sh
done
