/**
 * @file ExpChannel.cc
 * @brief Implementation of the ExpChannel class.
 */

#include "ExpChannel.h"
#include "includes.h"
#include "Position.h"
#include <fstream>
#include <cmath>
#include <vector>

using std::vector;

bool ExpChannel::init(cSimpleModule* module,
		const std::vector<std::vector<Position>>& msPositions, 
		std::map<int,Position>& neighbourPositions){
	// First, execute the parent init method of the Channel class
	Channel::init(module,msPositions,neighbourPositions);
	expMean = module->par("expMean");
	plExp = module->par("plExp");
	int run = std::stoi(ev.getConfig()->substituteVariables("${runnumber}"));
	std::string fname("./results/run_"+std::to_string(run)+"_coeff_table_down_"
			+std::to_string(bsId)+".dat");
	downValues.open(fname,std::ofstream::trunc);
	downValues << "TTI\t" << "BS\t" << "MS\t" << "RB\t" << "PL\t" << "Exp\t" << "Coeff" << "\n"; 
	fname = "./results/run_"+std::to_string(run)+"_coeff_table_up_"+std::to_string(bsId)+".dat";
	upValues.open(fname,std::ofstream::trunc);
	upValues << "TTI\t" << "Cell\t" << "MS\t" << "BS\t" << "RB\t" << "PL\t" << "Exp\t" << "Coeff" << "\n"; 
	recomputeCoefficients(msPositions);
	return true;
}

double ExpChannel::pathloss(Position sender, Position receiver){
	double distX = pow(sender.x-receiver.x,2);
	double distY = pow(sender.y-receiver.y,2);
	double dist = sqrt(distX+distY);
	return 10*plExp*log10(dist);
}

void ExpChannel::recomputeCoefficients(
		const vector<vector<Position>>& msPositions){
	auto val = simTime().dbl()/tti;
	// +1 because when this method is called, it computes the values for the 
	// NEXT TTI, not the current one.
	int tti = std::floor(val);
	size_t numBs = msPositions.size();
	// If the mobile station positions have not yet been stored locally, do so 
	// now.
	if(msPos.empty()){
		this->msPos = msPositions;
	}
	// Compute DOWN RB coefficients
	coeffDownTable.resize(numberOfMobileStations,
			vector<vector<vector<double>>>(numBs,
				vector<vector<double>>(timeSamples,
					vector<double>(upRBs))));
	double pl = 0.0;
	double exp = 0.0;
	for(size_t msIds=0; msIds<numberOfMobileStations; ++msIds){
		for(size_t bsIds=0; bsIds<numBs; ++bsIds){
			pl = pathloss(neighbourPositions[bsIds],msPositions[bsId][msIds]);
			for(size_t t=0; t<timeSamples; ++t){
				for(size_t rb=0; rb<downRBs; ++rb){
					exp = exponential(expMean);
					coeffDownTable[msIds][bsIds][t][rb] = pl * exp;
				}
			}
			for(size_t rb=0; rb<downRBs; ++rb){
				downValues << tti << "\t" << bsIds << "\t"
					<< msIds << "\t"
					<< rb << "\t"
					<< pl << "\t"
					<< coeffDownTable[msIds][bsIds][timeSamples-1][rb]/pl << "\t"
					<< coeffDownTable[msIds][bsIds][timeSamples-1][rb]
					<< std::endl;
			}
		}
	}
	// Compute UP RB coefficients
	coeffUpTable.resize(numBs,
			vector<vector<vector<vector<double>>>>(1));
	size_t numMs;
	for(size_t bsIds=0; bsIds<numBs; ++bsIds){
		numMs = msPositions[bsIds].size();
		coeffUpTable[bsIds][0].resize(numMs,
				vector<vector<double>>(timeSamples,
					vector<double>(upRBs)));
		for(size_t msIds=0; msIds<numMs; ++msIds){
			pl = pathloss(msPositions[bsIds][msIds],neighbourPositions[bsId]);
			for(size_t t=0; t<timeSamples; ++t){
				for(size_t rb=0; rb<upRBs; ++rb){
					exp = exponential(expMean);
					coeffUpTable[bsIds][0][msIds][t][rb] = pl * exp;
				}
			}
			for(size_t rb=0; rb<upRBs; ++rb){
				upValues << tti << "\t" << bsIds << "\t"
					<< msIds << "\t"
					<< bsId << "\t"
					<< rb << "\t"
					<< pl << "\t"
					<< coeffUpTable[bsIds][0][msIds][timeSamples-1][rb]/pl << "\t"
					<< coeffUpTable[bsIds][0][msIds][timeSamples-1][rb]
					<< std::endl;
			}
		}
	}
	// Compute D2D UP RB coefficients
	coeffUpD2DTable.resize(numBs,
			vector<vector<vector<vector<double>>>>(numberOfMobileStations));
	for(size_t bsIds=0; bsIds<numBs; ++bsIds){
		numMs = msPositions[bsIds].size();
		for(size_t recMsId=0; recMsId<numberOfMobileStations; ++recMsId){
			coeffUpD2DTable[bsIds][recMsId].resize(numMs,
					vector<vector<double>>(timeSamples,
						vector<double>(upRBs)));
			for(size_t msIds=0; msIds<numMs; ++msIds){
				for(size_t t=0; t<timeSamples; ++t){
					for(size_t rb=0; rb<upRBs; ++rb){
						pl = pathloss(msPositions[bsIds][msIds],msPositions[bsId][recMsId]);
						coeffUpD2DTable[bsIds][recMsId][msIds][t][rb] = pl * exponential(expMean);
					}
				}
			}
		}
	}
	// Compute D2D DOWN RB coefficients
	coeffDownD2DTable.resize(numBs,
			vector<vector<vector<vector<double>>>>(numberOfMobileStations));
	for(size_t bsIds=0; bsIds<numBs; ++bsIds){
		numMs = msPositions[bsIds].size();
		for(size_t recMsId=0; recMsId<numberOfMobileStations; ++recMsId){
			coeffDownD2DTable[bsIds][recMsId].resize(numMs,
					vector<vector<double>>(timeSamples,
						vector<double>(upRBs)));
			for(size_t msIds=0; msIds<numMs; ++msIds){
				for(size_t t=0; t<timeSamples; ++t){
					for(size_t rb=0; rb<downRBs; ++rb){
						pl = pathloss(msPositions[bsIds][msIds],msPositions[bsId][recMsId]);
						coeffDownD2DTable[bsIds][recMsId][msIds][t][rb] = pl * exponential(expMean);
					}
				}
			}
		}
	}
}

void ExpChannel::updateChannel(const vector<vector<Position>>& msPos){
	recomputeCoefficients(msPos);
}

void ExpChannel::handleMessage(cMessage* msg){
}

void ExpChannel::clearTransInfo(){
	Channel::clearTransInfo();
	recomputeCoefficients(msPos);
}

ExpChannel::~ExpChannel(){
	downValues.close();
	upValues.close();
}
