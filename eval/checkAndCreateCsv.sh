#!/bin/bash

function checkForErrors
{
	if ! test -f $2
	then
		echo $1 NOT DONE YET!!!
		return 1
	fi
	if grep -i warning $2
	then
		echo $1 HAS WARNINGS!!!
		return 1
	fi
	if ! grep -i stopwatch $2 > /dev/null
	then
		echo $1 HAS NOT COMPLETED!!!
		return 1
	else
		time=`grep -i stopwatch $2 | cut -f 2 -d ":" | tr -d " "`
		echo "$3;$4;$5;$time" >> lteLargeScale.csv
	fi
	return 0
}

function check
{
	for((rep=0;rep<5;rep++))
	do
		checkForErrors ${1} ${1}${rep}.0 $2 $3 $4 $5 $6 $7 || return
	done
	#echo $1 ok
}

function combineAndCheck
{
	if test $pdes = "par" -o $pdes = "com"
	then
		for sync in "" "Int"
		do
			syncStr=$sync
			if test -z "$syncStr"; then syncStr="--"; fi
			check results/$pdes${cells}cells$sync $pdes $cells $syncStr
		done
	else
		check results/$pdes${cells}cells $pdes $cells --
	fi
}

echo "PDES;size;SYNC;time" > lteLargeScale.csv
for pdes in seq par com
do
	for cells in 48 96 192 384 768 1536
	do
		combineAndCheck
	done
done

