#!/usr/bin/env zsh

file=lte.csv

echo "Runname; Runtime; " > $file
for x in {0..29}
do
    tmp=`cat hor48cells$x | grep "StopWatch" | cut -d ' ' -f 3`
    echo "hor48cells$x; $tmp; " >> $file
done

for i in seq48cells par48cells par48cellsInterval com48cells com48cellsInterval par48cells12LP par48cells12LPInterval
do
    for x in {0..29}
    do
        tmp=`cat results/$i$x.0 | grep "StopWatch" | cut -d ' ' -f 3`
        echo "$i$x; $tmp; " >> $file
    done
done
