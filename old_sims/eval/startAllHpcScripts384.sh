#!/usr/bin/env zsh

#bsub < hpcSeq384cells.sh
./startHpcCom384cells.sh
./startHpcCom384cellsInterval.sh
./startHpcPar384cells.sh
./startHpcPar384cellsInterval.sh
