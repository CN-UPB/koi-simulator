#!/usr/bin/python
import sys

#currently the parameters are set for the 768 cell example

#lines without offset
cell = 0
lineOffset = 0
maxLine = 1
lineCounter = 0

for line in range(0, 64):
    for row in range(0, 48):
        rowOffset = int(row) / 6
        print "**.cell[" + str(cell) + "].partition-id = " + str(rowOffset + lineOffset)
        cell += 1
    print ""

    if lineCounter >= maxLine-1:
        lineCounter = 0
        lineOffset += 8
        if maxLine == 2:
            maxLine = 1
        elif maxLine == 1:
            maxLine = 1
    else:
        lineCounter += 1

#lines with offset
cell = 3072
lineOffset = 0
maxLine = 1
lineCounter = 0

for line in range(0, 64):
    for row in range(0, 48):
        rowOffset = int(row) / 6
        print "**.cell[" + str(cell) + "].partition-id = " + str(rowOffset + lineOffset)
        cell += 1
    print ""

    if lineCounter >= maxLine-1:
        lineCounter = 0
        lineOffset += 8
        if maxLine == 2:
            maxLine = 1
        elif maxLine == 1:
            maxLine = 1
    else:
        lineCounter += 1
