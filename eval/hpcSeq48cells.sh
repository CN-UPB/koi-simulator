#!/usr/bin/env zsh

### Job name
#BSUB -J seq48cells

### File / path where output will be written, the %J is the job id
#BSUB -o seq48cells.%J

### (OFF) Different file for STDERR, if not to be merged with STDOUT
# #BSUB -e seq48cells.e%J

### Request the time you need for execution in minutes
### The format for the parameter is: [hour:]minute,
### that means for 80 minutes you could also use this: 1:20
#BSUB -W 00:10

### Request vitual memory you need for your job in MB
#BSUB -M 1024

### (OFF) Specify your mail address
# #BSUB -u sascha.schmerling@rwth-aachen.de

### Send a mail when job is scheduled and done
##BSUB -N

### Hybrid Job with N MPI Processes in groups to M processes per node
#BSUB -n 30
### Request a certain node type
#BSUB -m mpi-s

### Choose a MPI: either Open MPI or Intel MPI
### Use esub for Open MPI
#BSUB -a intelmpi
module switch openmpi intelmpi

### Set up omnet
export PATH=/home/ss136247/thesis/horizon-4/bin:$PATH
export LD_LIBRARY_PATH=/home/ss136247/thesis/horizon-4/lib:$LD_LIBRARY_PATH
export HOSTNAME
export HOST

### Change to the work directory
cd /home/ss136247/thesis/horizon-4/samples/abstractLTE

### Execute your application
cmd=""
for x in {0..29}
do
  tmp="-n 1 ./lteWrapper.sh seq48cells $x 0 -c seq48cells --use-thread-pool=false --parallel-simulation=false"
  if [ "$x" -ne 0 ];
    then
      cmd="$cmd : $tmp"
    else
      cmd="$tmp"
  fi
done
MYFLAGS=`echo $FLAGS_MPI_BATCH | sed 's/-n[p]* [0-9]*//g'`
echo "-------------"
echo $MPIEXEC
echo $FLAGS_MPI_BATCH
echo $MYFLAGS
echo $cmd
echo "-------------"
$MPIEXEC $MYFLAGS $cmd
