#!/usr/bin/env Rscript
# This script will execute evaluation scripts for the local experiment.
# It expects to find the KoI scripts/ directory in it's parent directory.
# The script library can be found in https://git.cs.upb.de/koi/simulations.git.

source("../scripts/kbest_delays.R")

numMs <- c(20,30,40,50,60,70,80)
numCells <- 9

delays_ms_exploration(numMs,10,"./results_i2i/",numCells,figname="delays_num_ms_i2i.pdf",tablename="delays_smaller_1ms_i2i.tab")
#delays_ms_exploration(numMs,1,"./results_o2o/",numCells,figname="delays_num_ms_o2o.pdf",tablename="delays_smaller_1ms_o2o.tab")
