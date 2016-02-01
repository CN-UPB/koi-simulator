#!/usr/bin/env zsh

python genHpcHor48cells.py -n 5
for x in {0..4}
do
    bsub < hpcHor48cells$x.sh
    rm hpcHor48cells$x.sh
done
