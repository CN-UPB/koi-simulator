/**
 * @file ExpChannel.cc
 * @brief Implementation of the ExpChannel class.
 */

#include "ExpChannel.h"
#include "includes.h"
#include "Position.h"
#include "util.h"
#include "VecNd.h"

#include <cmath>
#include <fstream>
#include <vector>

using std::vector;

bool ExpChannel::init(cSimpleModule* module,
		const std::vector<std::vector<Position>>& msPositions, 
		std::map<int,Position>& neighbourPositions){
	// First, execute the parent init method of the Channel class
	Channel::init(module,msPositions,neighbourPositions);
	expMean = module->par("expMean");
	plExp = module->par("plExp");
	initOffset = module->par("initOffset");
	std::string fname("coeff_table_down-"+std::to_string(bsId));
	downValues = std::move(getResultFile(fname));
	downValues << "TTI\t" << "BS\t" << "MS\t" << "RB\t" << "PL\t" << "Exp\t" << "Coeff" << "\n"; 
	fname = "coeff_table_up-"+std::to_string(bsId);
	upValues = std::move(getResultFile(fname));
	upValues << "TTI\t" << "Cell\t" << "MS\t" << "BS\t" << "RB\t" << "PL\t" << "Exp\t" << "Coeff" << "\n"; 
	recomputeCoefficients(msPositions);
	return true;
}

double ExpChannel::pathgain(Position sender, Position receiver){
	double distX = pow(sender.x-receiver.x,2);
	double distY = pow(sender.y-receiver.y,2);
	double dist = sqrt(distX+distY);
	double pl = (10*plExp*log10(dist));
	// Convert pathloss to linear scale and return gain instead of loss
	return std::pow(10,-pl/10);
}

void ExpChannel::recomputeCoefficients(
		const vector<vector<Position>>& msPositions){
	auto val = ((simTime().dbl()-initOffset)/tti).dbl();
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
			VectorNd<double,3>(numBs,
				VectorNd<double,2>(timeSamples,vector<double>(upRBs))));
	double pg = 0.0;
	double exp = 0.0;
	for(size_t msIds=0; msIds<numberOfMobileStations; ++msIds){
		for(size_t bsIds=0; bsIds<numBs; ++bsIds){
			pg = pathgain(neighbourPositions[bsIds],msPositions[bsId][msIds]);
			for(size_t t=0; t<timeSamples; ++t){
				for(size_t rb=0; rb<downRBs; ++rb){
					exp = exponential(expMean);
					coeffDownTable[msIds][bsIds][t][rb] = pg * exp;
				}
			}
			if(simTime()>initOffset){
				for(size_t rb=0; rb<downRBs; ++rb){
					downValues << tti << "\t" << bsIds << "\t"
						<< msIds << "\t"
						<< rb << "\t"
						<< pg << "\t"
						<< coeffDownTable[msIds][bsIds][timeSamples-1][rb]/pg << "\t"
						<< coeffDownTable[msIds][bsIds][timeSamples-1][rb]
						<< std::endl;
				}
			}
		}
	}
	// Compute UP RB coefficients
	coeffUpTable.resize(numBs,
			VectorNd<double,4>(1));
	size_t numMs;
	for(size_t bsIds=0; bsIds<numBs; ++bsIds){
		numMs = msPositions[bsIds].size();
		coeffUpTable[bsIds][0].resize(numMs,
				vector<vector<double>>(timeSamples,
					vector<double>(upRBs)));
		for(size_t msIds=0; msIds<numMs; ++msIds){
			pg = pathgain(msPositions[bsIds][msIds],neighbourPositions[bsId]);
			for(size_t t=0; t<timeSamples; ++t){
				for(size_t rb=0; rb<upRBs; ++rb){
					exp = exponential(expMean);
					coeffUpTable[bsIds][0][msIds][t][rb] = pg * exp;
				}
			}
			if(simTime()>initOffset){
				for(size_t rb=0; rb<upRBs; ++rb){
					upValues << tti << "\t" << bsIds << "\t"
						<< msIds << "\t"
						<< bsId << "\t"
						<< rb << "\t"
						<< pg << "\t"
						<< coeffUpTable[bsIds][0][msIds][timeSamples-1][rb]/pg << "\t"
						<< coeffUpTable[bsIds][0][msIds][timeSamples-1][rb]
						<< std::endl;
				}
			}
		}
	}
	if(d2dActive){
		// Compute D2D UP RB coefficients
		coeffUpD2DTable.resize(numBs,
				VectorNd<double,4>(numberOfMobileStations));
		for(size_t bsIds=0; bsIds<numBs; ++bsIds){
			numMs = msPositions[bsIds].size();
			for(size_t recMsId=0; recMsId<numberOfMobileStations; ++recMsId){
				coeffUpD2DTable[bsIds][recMsId].resize(numMs,
						vector<vector<double>>(timeSamples,
							vector<double>(upRBs)));
				for(size_t msIds=0; msIds<numMs; ++msIds){
					for(size_t t=0; t<timeSamples; ++t){
						for(size_t rb=0; rb<upRBs; ++rb){
							pg = pathgain(msPositions[bsIds][msIds],msPositions[bsId][recMsId]);
							coeffUpD2DTable[bsIds][recMsId][msIds][t][rb] = pg * exponential(expMean);
						}
					}
				}
			}
		}
		// Compute D2D DOWN RB coefficients
		coeffDownD2DTable.resize(numBs,
				VectorNd<double,4>(numberOfMobileStations));
		for(size_t bsIds=0; bsIds<numBs; ++bsIds){
			numMs = msPositions[bsIds].size();
			for(size_t recMsId=0; recMsId<numberOfMobileStations; ++recMsId){
				coeffDownD2DTable[bsIds][recMsId].resize(numMs,
						vector<vector<double>>(timeSamples,
							vector<double>(upRBs)));
				for(size_t msIds=0; msIds<numMs; ++msIds){
					for(size_t t=0; t<timeSamples; ++t){
						for(size_t rb=0; rb<downRBs; ++rb){
							pg = pathgain(msPositions[bsIds][msIds],msPositions[bsId][recMsId]);
							coeffDownD2DTable[bsIds][recMsId][msIds][t][rb] = pg * exponential(expMean);
						}
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
