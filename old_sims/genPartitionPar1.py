#!/bin/python

import os

#generates out of the com partition the data for genPartitionPar2.py script

#u need to adjust the range to the number of lp in the com scenerio
for i in range(0, 512):
    os.system("cat tmp.ini | grep \" = " + str(i) + "$\" >> test.ini")
    os.system("echo \"\" >> test.ini")
