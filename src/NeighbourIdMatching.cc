/*
 * NeigbourIdMatching.cc
 *
 *  Created on: Jul 17, 2013
 *      Author: Sascha Schmerling
 */

#include "NeighbourIdMatching.h"

using namespace omnetpp;

NeighbourIdMatching::NeighbourIdMatching(int ownBsId, int maxNumberOfNeighbours, cModule *cell)  {
    setupNeighbours(ownBsId, maxNumberOfNeighbours, cell);
    this->cell = cell;
}

void NeighbourIdMatching::setupNeighbours(int ownBsId, int maxNumberOfNeighbours, cModule *cell)  {
    //find the neighbours and store the pair (bsId, position in data structures) in a map
    this->ownBsId = ownBsId;
    insert(ownBsId, 0, -1); //store the own bs in the map
    for(int i = 0; i < maxNumberOfNeighbours; ++i)  {
        cGate *outGate = cell->gate(cell->findGate("toCell", i));
        if(outGate->isConnected())  {
            //get the id of the other cell
            cGate *targetGate = outGate->getNextGate();
            int bsId = targetGate->getOwnerModule()->par("bsId");
            insert(bsId, numberOfNeighbours(), i);
        }
    }
}

void NeighbourIdMatching::insert(int bsId, int dataStrId, int gateId)  {
    matching.insert(std::make_pair(bsId, std::make_pair(dataStrId, gateId)));
}

int NeighbourIdMatching::getDataStrId(int bsId)  {
    return (matching.at(bsId)).first;
}

int NeighbourIdMatching::getGateId(int bsId)  {
    return (matching.at(bsId)).second;
}

NeighbourMap *NeighbourIdMatching::getNeighbourMap()  {
    return &matching;
}

int NeighbourIdMatching::getNumberOfMS(int bsId){
	if(bsId==ownBsId){
		return cell->par("numberOfMobileStations");
	}
	else{
		int bsGate = getGateId(bsId);
		cGate *out = cell->gate(cell->findGate("toCell",bsGate));
		cGate *targetGate = out->getNextGate();
		return targetGate->getOwnerModule()->par("numberOfMobileStations");
	}
}

int NeighbourIdMatching::numberOfNeighbours()  {
    return matching.size();
}

void NeighbourIdMatching::showMatching()  {
    for(auto& i:matching){
			EV << "BsId " << i.first << ": DataStrId " << (i.second).first <<
				" GateId " << (i.second).second << endl;
		}
}
