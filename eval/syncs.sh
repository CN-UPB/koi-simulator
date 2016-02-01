#!/usr/bin/env zsh

for i in par48cells par48cellsInterval com48cells com48cellsInterval par48cells12LP par48cells12LPInterval
do
    result="0"
    for x in {0..29}
    do

        float events=0
        res=`cat results/$i$x.* | grep "Simulation time limit reached" | cut -d ' ' -f 11 | cut -c 2- | sed 's/.$//'`
        lines=("${(f)res}")
        for l in $lines
        do
            let events+=l
        done

        float syncs=`cat results/$i$x.0 | grep "syncs" | cut -d ' ' -f 8`
        tmp=events/syncs
        let result+=events/syncs
    done
    let result/=30
    echo "$i: $result"
done
