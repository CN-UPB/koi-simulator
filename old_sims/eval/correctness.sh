#!/usr/bin/env zsh

for filename in 48cellsTMP1.txt 48cellsTMP2.txt
do
    if [ -s $filename ]; then
        echo "TMP file $filename exists!"
        exit 1
    fi
done

echo "hor48cells"
for x in {0..29}
do
    `cat results/seq48cells$x.0 | grep " packets!" | sort > 48cellsTMP1.txt`
    `cat hor48cells$x | grep " packets!" | sort > 48cellsTMP2.txt`
    res=`diff 48cellsTMP1.txt 48cellsTMP2.txt`
    if [ -n "$res" ];
    then
        echo "hor48cells$x not equal!"
    fi
    `rm 48cellsTMP1.txt`
    `rm 48cellsTMP2.txt`
done

for i in par48cells par48cellsInterval com48cells com48cellsInterval par48cells12LP par48cells12LPInterval
do
    echo $i
    for x in {0..29}
    do
        `cat results/seq48cells$x.0 | grep " packets!" | sort > 48cellsTMP1.txt`
        `cat results/$i$x.* | grep " packets!" | sort > 48cellsTMP2.txt`
        res=`diff 48cellsTMP1.txt 48cellsTMP2.txt`
        if [ -n "$res" ];
        then
            echo "$i$x not equal!"
        fi
        `rm 48cellsTMP1.txt`
        `rm 48cellsTMP2.txt`
    done
done
