#!/usr/bin/env zsh

python genHpcCom.py -c com3072cells -np 256 -n 5 -horizondir /home/ms926272/repos/horizonParsim -sampledir=/home/ms926272/repos/abstractLTE -t 01:00 -sync cLbtsSync
for x in {0..4}
do
    bsub < hpc_com3072cells$x.sh
    rm hpc_com3072cells$x.sh
done
