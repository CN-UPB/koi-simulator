#!/usr/bin/env zsh

#bsub < hpcSeq96cells.sh
./startHpcCom96cells.sh
./startHpcCom96cellsInterval.sh
./startHpcPar96cells.sh
./startHpcPar96cellsInterval.sh
