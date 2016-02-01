#!/usr/bin/env zsh

bsub < hpcSeq48cells.sh
./startHpcCom48cells.sh
./startHpcCom48cellsInterval.sh
./startHpcHor48cells.sh
./startHpcPar48cells.sh
./startHpcPar48cells12LP.sh
./startHpcPar48cells12LPInterval.sh
./startHpcPar48cellsInterval.sh
