#!/usr/bin/env zsh

python genHpcPar.py -c par6144cells -np 6144 -n 5 -horizondir /home/ss136247/thesis/horizon-4 -sampledir=/home/ss136247/thesis/horizon-4/samples/abstractLTE -t 00:20 -sync cIntervalBasedSync -suffix Int
for x in {0..4}
do
    bsub < hpc_par6144cellsInt$x.sh
    rm hpc_par6144cellsInt$x.sh
done
