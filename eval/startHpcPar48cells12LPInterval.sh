#!/usr/bin/env zsh

python genHpcPar48cells12LPInterval.py -n 5
for x in {0..4}
do
    bsub < hpcPar48cells12LPInterval$x.sh
    rm hpcPar48cells12LPInterval$x.sh
done
