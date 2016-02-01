#!/usr/bin/env zsh

python genHpcCom.py -c com1536cells -np 128 -n 5 -horizondir /home/ms926272/repos/horizonParsim -sampledir=/home/ms926272/repos/abstractLTE -t 00:15 -sync cLbtsSync
for x in {0..4}
do
    bsub < hpc_com1536cells$x.sh
    rm hpc_com1536cells$x.sh
done
