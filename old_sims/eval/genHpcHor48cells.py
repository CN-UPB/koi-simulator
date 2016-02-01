#!/bin/python
import sys
import argparse

parser = argparse.ArgumentParser(description='Creates hpc script for the eval.')
parser.add_argument('-n', required=True, type=int,
                   help='the number of runs')

args = parser.parse_args()


for i in range(0, args.n):
    f = open("hpcHor48cells" + str(i) + ".sh", "w")

    f.write("#!/usr/bin/env zsh\n")
    f.write("\n")

    f.write("### Job name\n")
    f.write("#BSUB -J hor48cells" + str(i) + "\n")
    f.write("\n")

    f.write("### File / path where output will be written, the %J is the job id\n")
    f.write("#BSUB -o hor48cells" + str(i) + "\n")
    f.write("\n")

    f.write("### (OFF) Different file for STDERR, if not to be merged with STDOUT\n")
    f.write("# #BSUB -e hor48cells" + str(i) + ".e\n")
    f.write("\n")

    f.write("### Request the time you need for execution in minutes\n")
    f.write("### The format for the parameter is: [hour:]minute,\n")
    f.write("### that means for 80 minutes you could also use this: 1:20\n")
    f.write("#BSUB -W 00:03\n")
    f.write("\n")

    f.write("### Request vitual memory you need for your job in MB\n")
    f.write("#BSUB -M 1024\n")
    f.write("\n")

    f.write("### (OFF) Specify your mail address\n")
    f.write("# #BSUB -u sascha.schmerling@rwth-aachen.de\n")
    f.write("\n")

    f.write("### Send a mail when job is scheduled and done\n")
    f.write("##BSUB -N\n")
    f.write("\n")

    f.write("### Request a certaion node type\n")
    f.write("#BSUB -m mpi-s\n")
    f.write("### Use node exclusive\n")
    f.write("#BSUB -x \n")
    f.write("\n")

    f.write("### Multi threaded\n")
    f.write("#BSUB -a openmp\n")
    f.write("### with this number of threads\n")
    f.write("#BSUB -n 12\n")
    f.write("\n")

    f.write("### Set up omnet\n")
    f.write("export PATH=/home/ss136247/thesis/horizon-4/bin:$PATH\n")
    f.write("export LD_LIBRARY_PATH=/home/ss136247/thesis/horizon-4/lib:$LD_LIBRARY_PATH\n")
    f.write("export HOSTNAME\n")
    f.write("export HOST\n")
    f.write("\n")

    f.write("### Change to the work directory\n")
    f.write("cd /home/ss136247/thesis/horizon-4/samples/abstractLTE\n")
    f.write("\n")

    f.write(" ### Execute your application\n")
    f.write("time ./abstractLTE -c hor48cells --use-thread-pool=true --parallel-simulation=false -r " + str(i) + "\n")

    f.close()

