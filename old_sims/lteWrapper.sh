#!/usr/bin/env zsh

name=$1
rep=$2
rank=$3
shift
shift
shift
./abstractLTE $@ -r $rep > eval/results/$name$rep.$rank
