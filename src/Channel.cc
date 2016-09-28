/*
 * Channel.cc
 *
 *  Created on: Jul 15, 2014
 *      Author: Thomas Prinz
 * 
 */
#include "Channel.h"
#include "MessageTypes.h"

using std::vector;

const double Channel::speedOfLight = 299792458.0;

bool Channel::init(cSimpleModule* module,
		const vector<vector<Position>>& msPositions, 
		std::map<int,Position>& neighbourPositions){
	bsId = module->par("bsId");
	downRBs = module->par("downResourceBlocks");
	initModule = module;
	maxNumberOfNeighbours = module->par("maxNumberOfNeighbours");
	neighbourPositions = neighbourPositions;
	numberOfMobileStations = module->par("numberOfMobileStations");
	tti = module->par("tti");
	upRBs = module->par("upResourceBlocks");

	// Find the neighbours and store the pair (bsId, position in data structures) in a map
	cModule *cell = module->getParentModule()->getParentModule();
	neighbourIdMatching = new NeighbourIdMatching(bsId, maxNumberOfNeighbours, cell);
	// Get Playground size from cell module:
	sizeX = cell->par("playgroundSizeX");
	sizeY = cell->par("playgroundSizeY");
	// Resize the vectors for the per-RB trans infos to the number of 
	// resource blocks.
	transInfo.first.resize(upRBs);
	transInfo.second.resize(downRBs);
	// There are a number of dynamically allocated member variables which 
	// are only allocated in this init method, which need only be freed 
	// in the destructor iff init has actually been called.
	initialized = true;
}

void Channel::addTransInfo(TransInfo* trans){
	int dir = trans->getMessageDirection();
	if(dir==MessageDirection::up || dir==MessageDirection::d2dUp){
		transInfo.first[trans->getRb()].push_front(trans);
	}
	else if(dir==MessageDirection::down || dir==MessageDirection::d2dDown){
		transInfo.second[trans->getRb()].push_front(trans);
	}
}
