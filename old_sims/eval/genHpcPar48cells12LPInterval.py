#!/bin/python
import sys
import argparse

parser = argparse.ArgumentParser(description='Creates hpc script for the eval.')
parser.add_argument('-n', required=True, type=int,
                   help='the number of runs')

args = parser.parse_args()


for i in range(0, args.n):
    f = open("hpcPar48cells12LPInterval" + str(i) + ".sh", "w")

    f.write("#!/usr/bin/env zsh\n")
    f.write("\n")

    f.write("### Job name\n")
    f.write("#BSUB -J par48cells12LPInterval" + str(i) + "\n")
    f.write("\n")

    f.write("### File / path where output will be written, the %J is the job id\n")
    f.write("#BSUB -o par48cells12LPInterval" + str(i) + ".%J\n")
    f.write("\n")

    f.write("### (OFF) Different file for STDERR, if not to be merged with STDOUT\n")
    f.write("# #BSUB -e par48cells12LPInterval" + str(i) + ".e%J2\n")
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

    f.write("### Hybrid Job with N MPI Processes in groups to M processes per node\n")
    f.write("#BSUB -n 12\n")
    f.write("#BSUB -R \"span[ptile=12]\"\n")
    f.write("### Request a certaion node type\n")
    f.write("#BSUB -m mpi-s\n")
    f.write("\n")

    f.write("### Choose a MPI: either Open MPI or Intel MPI\n")
    f.write("### Use esub for Open MPI\n")
    f.write("#BSUB -a intelmpi\n")
    f.write("module switch openmpi intelmpi\n")
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

    f.write("### Execute your application\n")
    f.write("cmd=\"\"\n");
    f.write("for x in {0..11}\n");
    f.write("do\n");
    f.write("  tmp=\"-n 1 ./lteWrapper.sh par48cells12LPInterval " + str(i) + " $x -c par48cells12LP --use-thread-pool=false --parallel-simulation=true --parsim-synchronization-class=cIntervalBasedSync\"\n");
    f.write("  if [ \"$x\" -ne 0 ];\n");
    f.write("   then\n");
    f.write("   cmd=\"$cmd : $tmp\"\n");
    f.write("   else\n");
    f.write("   cmd=\"$tmp\"\n");
    f.write("  fi\n");
    f.write("done\n");
    f.write("MYFLAGS=`echo $FLAGS_MPI_BATCH | sed \'s/-n[p]* [0-9]*//g\'`\n");
    f.write("echo \"-------------\"\n");
    f.write("echo $MPIEXEC\n");
    f.write("echo $FLAGS_MPI_BATCH\n");
    f.write("echo $MYFLAGS\n");
    f.write("echo $cmd\n");
    f.write("echo \"-------------\"\n");
    f.write("$MPIEXEC $MYFLAGS $cmd\n");

    f.close()
