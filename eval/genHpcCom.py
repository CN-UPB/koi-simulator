#!/bin/python
import sys
import argparse

parser = argparse.ArgumentParser(description='Creates hpc script for the eval.')
parser.add_argument('-np', required=True, type=int,
                   help='Number of MPI Processes')
parser.add_argument('-t', required=True, type=str,
                   help='Runtime e.g. 00:05 for 5 min')
parser.add_argument('-c', required=True, type=str,
                   help='Configname')
parser.add_argument('-n', required=True, type=int,
                   help='the number of runs')
parser.add_argument('-horizondir', required=True, type=str,
                   help='horizon directory')
parser.add_argument('-sampledir', required=True, type=str,
                   help='Sample directory')
parser.add_argument('-sync', required=True, type=str,
                   help='synchronization class that will be used')
parser.add_argument('-suffix', required=False, type=str,
                   help='added to the output file name', default='')

args = parser.parse_args()


for i in range(0, args.n):
    f = open("hpc_" + args.c + args.suffix + str(i) + ".sh", "w")

    f.write("#!/usr/bin/env zsh\n")
    f.write("\n")

    f.write("### Job name\n")
    f.write("#BSUB -J " + args.c + args.suffix + str(i) + "\n")
    f.write("\n")

    f.write("### File / path where output will be written, the %J is the job id\n")
    f.write("#BSUB -o " + args.c + args.suffix + str(i) + ".%J\n")
    f.write("\n")

    f.write("### (OFF) Different file for STDERR, if not to be merged with STDOUT\n")
    f.write("# #BSUB -e " + args.c + args.suffix + str(i) + ".e%J\n")
    f.write("\n")

    f.write("### Request the time you need for execution in minutes\n")
    f.write("### The format for the parameter is: [hour:]minute,\n")
    f.write("### that means for 80 minutes you could also use this: 1:20\n")
    f.write("#BSUB -W " + args.t + "\n")
    f.write("\n")

    f.write("### Request vitual memory you need for your job in MB\n")
    f.write("#BSUB -M 1024\n")
    f.write("\n")

    f.write("### (OFF) Specify your mail address\n")
    f.write("#BSUB -u stoffers@comsys.rwth-aachen.de\n")
    f.write("\n")

    f.write("### Send a mail when job is scheduled and done\n")
    f.write("#BSUB -N\n")
    f.write("\n")

    f.write("### Hybrid Job with N MPI Processes in groups to M processes per node\n")
    f.write("#BSUB -n " + str(args.np) + "\n")
    f.write("#BSUB -R \"span[ptile=1]\"\n")
    f.write("### Request a certaion node type\n")
    f.write("#BSUB -m mpi-s\n")
    f.write("### Use nodes exclusive\n")
    f.write("#BSUB -x \n")
    f.write("### Each MPI process with T Threads\n")
    f.write("export OMP_NUM_THREADS=12\n")
    f.write("\n")

    f.write("### Choose a MPI: either Open MPI or Intel MPI\n")
    f.write("### Use esub for Open MPI\n")
    f.write("#BSUB -a intelmpi\n")
    f.write("module switch openmpi intelmpi\n")
    f.write("\n")

    f.write("### Set up omnet\n")
    f.write("export PATH=" + args.horizondir + "/bin:$PATH\n")
    f.write("export LD_LIBRARY_PATH=/home/ss136247/thesis/horizon-4/lib:$LD_LIBRARY_PATH\n")
    f.write("export HOSTNAME\n")
    f.write("export HOST\n")
    f.write("\n")

    f.write("### Change to the work directory\n")
    f.write("cd " + args.sampledir + "\n")
    f.write("\n")

    f.write("### Execute your application\n")
    f.write("cmd=\"\"\n");
    f.write("for x in {0.." + str(args.np-1) + "}\n");
    f.write("do\n");
    f.write("  tmp=\"-n 1 ./lteWrapper.sh " + args.c + args.suffix + " " + str(i) + " $x -c " + args.c + " --use-thread-pool=true --parallel-simulation=true --parsim-synchronization-class=" + args.sync + "\"\n");
    f.write("  if [ \"$x\" -ne 0 ];\n");
    f.write("	then\n");
    f.write("	cmd=\"$cmd : $tmp\"\n");
    f.write("	else\n");
    f.write("	cmd=\"$tmp\"\n");
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

