#!/usr/bin/env Rscript
# This script will execute evaluation scripts for the local experiment.
# It expects to find the KoI scripts/ directory in it's parent directory.
# The script library can be found in https://git.cs.upb.de/koi/simulations.git.

source("../scripts/kbest_delays.R")

#numMs <- c(20,30,40,50,60,70,80)
numMs <- c(20,30)
numCells <- 18
numFloors <- 2

delays_ms_multifloor(numMs,1,"./results_exploration/",numCells,9,figname="delays_num_ms_exploration.pdf",tablename="delays_smaller_1ms_exploration.tab")
