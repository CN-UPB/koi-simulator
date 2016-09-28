/*
 * Channel.cc
 *
 *  Created on: Jul 15, 2014
 *      Author: Thomas Prinz
 * 
 */
#include "Channel.h"
#include "MessageTypes.h"
#include "NeighbourIdMatching.h"
#include "Position.h"
#include "TransInfo_m.h"
#include <forward_list>
#include <vector>

using std::forward_list;
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
	SINRcounter = 0;
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
	return initialized;
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

double Channel::calcInterference(forward_list<TransInfo*>& interferers,
		int rb,
		int receiverId,
		int SINRCounter,
		MessageDirection dir){
	double interference = 0.0;
	forward_list<TransInfo*>::iterator prev(interferers.before_begin());
	for(auto it = interferers.begin(); it!=interferers.end(); prev=it++){
		if((*it)->getCreationTime()>=simTime()-tti){
			if(receiverId==-1){
				// Interference at the local base station
				interference += (*it)->getPower() 
					* coeffUpTable[(*it)->getBsId()][0][(*it)->getMsId()][SINRCounter][rb];
			}
			else{
				// Interference at local mobile stations
				if((*it)->getMessageDirection()==MessageDirection::down){
					// Interference from a neighbouring BS
					interference += (*it)->getPower() * coeffDownTable[receiverId][(*it)->getBsId()][SINRCounter][rb];
				}
				else if((*it)->getMessageDirection()==MessageDirection::d2dDown){
					// Interference from a MS transmitting D2D on the same down resource 
					// block
					interference += (*it)->getPower() 
						* coeffDownD2DTable[(*it)->getBsId()][receiverId][(*it)->getMsId()][SINRCounter][rb];
				}
				else if((*it)->getMessageDirection()==MessageDirection::d2dUp
						|| (*it)->getMessageDirection()==MessageDirection::up){
					interference += (*it)->getPower() 
						* coeffUpD2DTable[(*it)->getBsId()][receiverId][(*it)->getMsId()][SINRCounter][rb];
				}
			}
		}
		else{
			delete *it;
			interferers.erase_after(prev);
			it=prev;
		}
	}
	interference += getTermalNoise(300,180000);
	return interference;
}

double Channel::calcUpSINR(int RB, 
		int msId,
		double transPower){
	int SINRCounter = 3; //originally set to std::round( simTime().dbl() * 1000.0* 4.0) - 1 
	double received = 0;
	double interference = calcInterference(transInfo.first[RB],
			RB,-1,SINRCounter,MessageDirection::up);
	received = transPower * coeffUpTable[bsId][0][msId][SINRcounter][RB];
	// Convert to db scale
	return 10 * log10( received / interference );
}

double Channel::calcDownSINR(int RB, 
		int msId,
		double transPower){
	int SINRCounter = 3; //originally set to std::round( simTime().dbl() * 1000.0* 4.0) - 1 
	double received = 0;
	double interference = calcInterference(transInfo.second[RB],
			RB,msId,SINRCounter,MessageDirection::down);
	received = transPower * coeffDownTable[msId][bsId][SINRcounter][RB];
	// Convert to db scale
	return 10 * log10( received / interference );
}

double Channel::calcD2DSINR(int RB, 
		int sendMsID,
		int receiveMsId,
		MessageDirection direction,
		double transPower){
	int SINRCounter = 3; //originally set to std::round( simTime().dbl() * 1000.0* 4.0) - 1 
	double received = 0;
	double interference = 0.0;
	if(direction==MessageDirection::d2dDown){
		received = transPower * coeffDownD2DTable[bsId][receiveMsId][sendMsID][SINRCounter][RB];
		interference = calcInterference(transInfo.second[RB],
				RB,receiveMsId,SINRCounter,MessageDirection::down);
	}
	else{
		received = transPower * coeffUpD2DTable[bsId][receiveMsId][sendMsID][SINRCounter][RB];
		interference = calcInterference(transInfo.first[RB],
				RB,receiveMsId,SINRCounter,MessageDirection::up);
	}
	// Convert to db scale
	return 10 * log10( received / interference );
}

double Channel::calcAvgUpSINR(int RB, 
            int msId,
            double transPower){
  // Compute SINR for MS->BS communication
  double res = calcUpSINR(RB,msId,transPower);
  // Add the average over all possible D2D connections to all local MS
  for(int i=0; i<numberOfMobileStations;++i){
    if(i!=msId){
      res += calcD2DSINR(RB,msId,i,MessageDirection::d2dUp,transPower);
    }
  }
  // Return the average over all SINR values
  // numberOfMobileStations is the number of values from the D2D computations,
  // +1 from the MS->BS computation.
  return res/(numberOfMobileStations+1);
}

double Channel::calcAvgDownSINR(int RB, 
            double transPower){
  // Compute SINR for DOWN resource blocks for base stations.
  double res = 0.0;
  // Average over SINR values for all local mobile stations
  for(int i=0; i<numberOfMobileStations;++i){
    res += calcDownSINR(RB,i,transPower);
  }
  return res/numberOfMobileStations;
}

double Channel::calcAvgD2DDownSINR(int RB, 
            int msId,
            double transPower){
  // Calculate average SINR for DOWN resource blocks when used for D2D by MS
  double res = 0.0;
  // Average over SINR values for all local mobile stations
  for(int i=0; i<numberOfMobileStations;++i){
    if(msId!=i){
      res += calcD2DSINR(RB,msId,i,MessageDirection::d2dDown,transPower);
    }
  }
  return res/numberOfMobileStations;
}

// Johnson Nyquist Noise
double Channel::getTermalNoise(double temp, double bandwidth){
	return (temp * bandwidth * 1.3806488e-23);
}

