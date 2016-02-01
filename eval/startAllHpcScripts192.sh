#!/usr/bin/env zsh

bsub < hpcSeq192cells.sh
./startHpcCom192cells.sh
./startHpcCom192cellsInterval.sh
./startHpcPar192cells.sh
./startHpcPar192cellsInterval.sh
