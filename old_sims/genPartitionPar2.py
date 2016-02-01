#!/usr/bin/python

#u need to delete the previous partition value before execution this script
#for example with the vim command :%s/ = [0-9]*$/ = /g
#generates out of genPartition1.py output the right par config

import os
import sys

dat = open("test.ini", "r+")
out = open("par-par.ini", "w+")
counter = 0
lines = dat.read().splitlines()
for line in lines:
    if(line != ''):
        out.write(line + str(counter) + "\n")
        counter+=1
    else:
        out.write("\n")
dat.close()
