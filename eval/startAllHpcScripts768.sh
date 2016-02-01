#!/usr/bin/env zsh

#bsub < hpcSeq768cells.sh
./startHpcCom768cells.sh
./startHpcCom768cellsInterval.sh
./startHpcPar768cells.sh
./startHpcPar768cellsInterval.sh
