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

#include <stdexcept>
#include <forward_list>
#include <vector>

using std::forward_list;

const double Channel::speedOfLight = 299792458.0;

bool Channel::init(cSimpleModule* module,
		const std::vector<std::vector<Position>>& msPositions, 
		const std::map<int,Position>& neighbourPositions){
	chnBandwidth = module->par("chnBandwidth");
	rbBandwidth = module->par("bandwidthPerRB");
	bsId = module->par("bsId");
	considerInterference = module->par("considerInterference");
	downRBs = module->par("downResourceBlocks");
	d2dActive = module->par("d2dActive");
	debug = module->par("debug");
	initModule = module;
	maxNumberOfNeighbours = module->par("maxNumberOfNeighbours");
	this->neighbourPositions = neighbourPositions;
	numberOfMobileStations = module->par("numberOfMobileStations");
	tti = module->par("tti");
	upRBs = module->par("upResourceBlocks");
	msPos = msPositions;

	// Find the neighbours and store the pair (bsId, position in data structures) in a map
	cModule *cell = module->getParentModule()->getParentModule();
	neighbourIdMatching = new NeighbourIdMatching(bsId, maxNumberOfNeighbours, cell);
	// Get Playground size from cell module:
	sizeX = cell->par("playgroundSizeX");
	sizeY = cell->par("playgroundSizeY");
	// Resize the vectors for the per-RB trans infos to the number of 
	// resource blocks.
	if(considerInterference){
		transInfo.first.resize(upRBs);
		transInfo.second.resize(downRBs);
	}
	// There are a number of dynamically allocated member variables which 
	// are only allocated in this init method, which need only be freed 
	// in the destructor iff init has actually been called.
	initialized = true;
	return initialized;
}

void Channel::addTransInfo(TransInfo* trans){
	if(considerInterference){
		int dir = trans->getMessageDirection();
		if(dir==MessageDirection::up || dir==MessageDirection::d2dUp){
			transInfo.first[trans->getRb()].push_front(trans);
		}
		else if(dir==MessageDirection::down || dir==MessageDirection::d2dDown){
			transInfo.second[trans->getRb()].push_front(trans);
		}
	}
	else{
		// We don't consider interference, so delete the message
		delete trans;
	}
}

double Channel::calcInterference(const forward_list<TransInfo*>& interferers,
		int rb,
		int receiverId){
	double interference = 0.0;
	if(considerInterference){
		for(auto interferer:interferers){
			if(receiverId==-1){
				// Interference at the local base station
				interference += interferer->getPower() 
					* coeffUpTable[interferer->getBsId()][0][interferer->getMsId()][rb];
			}
			else{
				// Interference at local mobile stations
				if(interferer->getMessageDirection()==MessageDirection::down){
					// Interference from a neighbouring BS
					interference += interferer->getPower() * coeffDownTable[receiverId][interferer->getBsId()][rb];
				}
				else if(interferer->getMessageDirection()==MessageDirection::d2dDown){
					// Interference from a MS transmitting D2D on the same down resource 
					// block
					interference += interferer->getPower() 
						* coeffDownD2DTable[interferer->getBsId()][receiverId][interferer->getMsId()][rb];
				}
				else if(interferer->getMessageDirection()==MessageDirection::d2dUp
						|| (interferer->getMessageDirection()==MessageDirection::up 
							&& d2dActive)){
					interference += interferer->getPower() 
						* coeffUpD2DTable[interferer->getBsId()][receiverId][interferer->getMsId()][rb];
				}
			}
		}
	}
	interference += getTermalNoise(300,chnBandwidth);
	return interference;
}

double Channel::calcUpSINR(int RB, 
		int msId,
		double transPower){
	double received = 0;
	double interference = calcInterference(transInfo.first[RB],
			RB,-1);
	received = transPower * coeffUpTable[bsId][0][msId][RB];
	// Convert to db scale
	return 10 * log10( received / interference );
}

double Channel::calcDownSINR(int RB, 
		int msId,
		double transPower){
	double received = 0;
	double interference = calcInterference(transInfo.second[RB],
			RB,msId);
	received = transPower * coeffDownTable[msId][bsId][RB];
	// Convert to db scale
	return 10 * log10( received / interference );
}

double Channel::senseUpSINR(int RB, 
		int msId,
		double transPower){
	double received = 0;
	received = transPower * coeffUpTable[bsId][0][msId][RB];
	// Convert to db scale
	return 10 * log10( received / getTermalNoise(300,chnBandwidth) );
}

double Channel::senseDownSINR(int RB, 
		int msId,
		double transPower){
	double received = 0;
	received = transPower * coeffDownTable[msId][bsId][RB];
	// Convert to db scale
	return 10 * log10( received / getTermalNoise(300,chnBandwidth) );
}

double Channel::calcD2DSINR(int RB, 
		int sendMsID,
		int receiveMsId,
		MessageDirection direction,
		double transPower){
	double received = 0;
	double interference = 0.0;
	if(direction==MessageDirection::d2dDown){
		received = transPower * coeffDownD2DTable[bsId][receiveMsId][sendMsID][RB];
		interference = calcInterference(transInfo.second[RB],
				RB,receiveMsId);
	}
	else{
		received = transPower * coeffUpD2DTable[bsId][receiveMsId][sendMsID][RB];
		interference = calcInterference(transInfo.first[RB],
				RB,receiveMsId);
	}
	// Convert to db scale
	return 10 * log10( received / interference );
}

double Channel::calcAvgUpSINR(int RB, 
            int msId,
            double transPower){
  // Compute SINR for MS->BS communication
  double res = calcUpSINR(RB,msId,transPower);
	if(d2dActive){
		// Add the average over all possible D2D connections to all local MS
		for(int i=0; i<numberOfMobileStations;++i){
			if(i!=msId){
				res += calcD2DSINR(RB,msId,i,MessageDirection::d2dUp,transPower);
			}
		}
		// Compute the average over all SINR values
		// numberOfMobileStations is the number of values from the D2D computations,
		// +1 from the MS->BS computation.
		res = res/(numberOfMobileStations+1);
	}
  return res;
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
	if(d2dActive){
		// Average over SINR values for all local mobile stations
		for(int i=0; i<numberOfMobileStations;++i){
			if(msId!=i){
				res += calcD2DSINR(RB,msId,i,MessageDirection::d2dDown,transPower);
			}
		}
	}
  return res/numberOfMobileStations;
}

void Channel::clearTransInfo(){
	for(auto& currList:transInfo.first){
		for(auto inf:currList){
			delete inf;
		}
		currList.clear();
	}
	for(auto& currList:transInfo.second){
		for(auto inf:currList){
			delete inf;
		}
		currList.clear();
	}
}

double Channel::calcLongtermUpSINR(int rb, int msId, double transPower){
	throw std::runtime_error("Long term SINR calculation not implemented.");
}

double Channel::calcLongtermDownSINR(int rb, int msId, double transPower){
	throw std::runtime_error("Long term SINR calculation not implemented.");
}

// Johnson Nyquist Noise
double Channel::getTermalNoise(double temp, double bandwidth){
	return (temp * bandwidth * 1.3806488e-23);
}

std::ostream& Channel::printCoeffUpTables(std::ostream& out){
	out << "BS" << "\t" << "MS" << "\t" << "RB" << "\t" << "Coeff" << std::endl;
	for(size_t i=0; i<coeffUpTable.size(); i++){
		for(size_t j=0; j<coeffUpTable[i][0].size(); j++){
			for(size_t r=0; r<downRBs; r++){
				out << i << "\t" << j << "\t" << r 
					<< "\t" << coeffUpTable[i][0][j][r] << std::endl;
			}
		}
	}
	return out;
}

std::ostream& Channel::printCoeffDownTables(std::ostream& out){
	out << "MS" << "\t" << "BS" << "\t" << "RB" 
		<< "\t" << "Coeff" << std::endl;
	for(size_t i=0; i<coeffDownTable.size(); i++){
		for(size_t j=0; j<coeffDownTable[i].size(); j++){
			for(size_t r=0; r<downRBs; r++){
				out << i << "\t" << j << "\t" << r
					<< "\t" << coeffDownTable[i][j][r] << std::endl;
			}
		}
	}
	return out;
}

void Channel::recomputePerTTIValues(){
	// Empty by default. Override this method if a channel needs any values 
	// recomputed at the beginning of a TTI.
}

Channel::~Channel(){
	if(initialized){
		// All of the following member variables are only allocated 
		// in the init method, not the constructor. As a result,
		// the destructor might be called with all of the following 
		// variables uninitialized, leading to errors. Thus, they 
		// are only freed when init has been called, as indicated by 
		// the initialized variable.
		delete neighbourIdMatching; 
	}
}
