#!/usr/bin/env Rscript
# This script will execute evaluation scripts for the local experiment.
# It expects to find the KoI scripts/ directory in it's parent directory.
# The script library can be found in https://git.cs.upb.de/koi/simulations.git.

source("../scripts/kbest_delays.R")

numCells <- 9
numMs <- c(20,30,40,50,60,70,80,90,100)

delays_ms_exploration(numMs,1,"./results_exploration/",numCells,
                      "delays_20_to_100_ms_exploration.pdf",
                      "delays_smaller_1ms_20_to_100_ms_exploration.tab")
delays_ms_exploration(numMs,30,"./results_multirun/",numCells,
                      "delays_20_to_100_ms_30_repeats.pdf",
                      "delays_smaller_1ms_20_to_100_ms_30_repeats.tab")