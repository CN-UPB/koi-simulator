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
using std::vector;

const double Channel::speedOfLight = 299792458.0;

bool Channel::init(cSimpleModule* module,
		const vector<vector<Position>>& msPositions, 
		std::map<int,Position>& neighbourPositions){
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
	SINRcounter = 0;
	tti = module->par("tti");
	upRBs = module->par("upResourceBlocks");

	// Find the neighbours and store the pair (bsId, position in data structures) in a map
	cModule *cell = module->getParentModule()->getParentModule();
	neighbourIdMatching = new NeighbourIdMatching(bsId, maxNumberOfNeighbours, cell);
	// Get Playground size from cell module:
	sizeX = cell->par("playgroundSizeX");
	sizeY = cell->par("playgroundSizeY");
	// Get position resend interval (Stationary MS assumed during this interval)
	timeSamples = module->par("positionResendInterval");
	// 4 Samples per TTI; for smooth Fourier transform
	// scale according to positionResendInterval; timeSamples should be comparable
	// to sim-time-limit
	timeSamples *= 4; 
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

double Channel::calcInterference(forward_list<TransInfo*>& interferers,
		int rb,
		int receiverId,
		int SINRCounter,
		MessageDirection dir){
	double interference = 0.0;
	forward_list<TransInfo*>::iterator prev(interferers.before_begin());
	if(considerInterference){
		for(auto it = interferers.begin(); it!=interferers.end(); prev=it++){
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
						|| ((*it)->getMessageDirection()==MessageDirection::up 
							&& d2dActive)){
					interference += (*it)->getPower() 
						* coeffUpD2DTable[(*it)->getBsId()][receiverId][(*it)->getMsId()][SINRCounter][rb];
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

double Channel::senseUpSINR(int RB, 
		int msId,
		double transPower){
	int SINRCounter = 3; //originally set to std::round( simTime().dbl() * 1000.0* 4.0) - 1 
	double received = 0;
	received = transPower * coeffUpTable[bsId][0][msId][SINRcounter][RB];
	// Convert to db scale
	return 10 * log10( received / getTermalNoise(300,chnBandwidth) );
}

double Channel::senseDownSINR(int RB, 
		int msId,
		double transPower){
	int SINRCounter = 3; //originally set to std::round( simTime().dbl() * 1000.0* 4.0) - 1 
	double received = 0;
	received = transPower * coeffDownTable[msId][bsId][SINRcounter][RB];
	// Convert to db scale
	return 10 * log10( received / getTermalNoise(300,chnBandwidth) );
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
	for(auto iterRb = transInfo.first.begin(); iterRb!=transInfo.first.end();
			++iterRb){
		forward_list<TransInfo*>& currList(*iterRb);
		for(auto inf:currList){
			delete inf;
		}
		currList.clear();
	}
	for(auto iterRb = transInfo.second.begin(); iterRb!=transInfo.second.end();
			++iterRb){
		forward_list<TransInfo*>& currList(*iterRb);
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
	out << "BS" << "\t" << "MS" << "\t" << "Time" << "\t" << "RB"
		<< "\t" << "Coeff" << std::endl;
	for(size_t i=0; i<coeffUpTable.size(); i++){
		for(size_t j=0; j<coeffUpTable[i][0].size(); j++){
			for(size_t t=0; t<timeSamples; t++){
				for(size_t r=0; r<downRBs; r++){
					out << i << "\t" << j << "\t" << t << "\t" << r 
						<< "\t" << coeffUpTable[i][0][j][t][r] << std::endl;
				}

			}
		}
	}
	return out;
}

std::ostream& Channel::printCoeffDownTables(std::ostream& out){
	out << "MS" << "\t" << "BS" << "\t" << "Time" << "\t" << "RB" 
		<< "\t" << "Coeff" << std::endl;
	for(size_t i=0; i<coeffDownTable.size(); i++){
		for(size_t j=0; j<coeffDownTable[i].size(); j++){
			for(size_t t=0; t<timeSamples; t++){
				for(size_t r=0; r<downRBs; r++){
					out << i << "\t" << j << "\t" << t << "\t" << r
						<< "\t" << coeffDownTable[i][j][t][r] << std::endl;
				}

			}
		}
	}
	return out;
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
