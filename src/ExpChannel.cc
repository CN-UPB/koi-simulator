/**
 * @file ExpChannel.cc
 * @brief Implementation of the ExpChannel class.
 */

#include "ExpChannel.h"
#include "includes.h"
#include "Position.h"
#include <vector>

using std::vector;

bool ExpChannel::init(cSimpleModule* module,
		const std::vector<std::vector<Position>>& msPositions, 
		std::map<int,Position>& neighbourPositions){
	// First, execute the parent init method of the Channel class
	Channel::init(module,msPositions,neighbourPositions);
	expMean = module->par("expMean");
	plExp = module->par("plExp");
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
	size_t numBs = msPositions.size();
	// Compute DOWN RB coefficients
	coeffDownTable.resize(numberOfMobileStations,
			vector<vector<vector<double>>>(numBs,
				vector<vector<double>>(timeSamples,
					vector<double>(upRBs))));
	double pl = 0.0;
	for(size_t msIds=0; msIds<numberOfMobileStations; ++msIds){
		for(size_t bsIds=0; bsIds<numBs; ++bsIds){
			for(size_t t=0; t<timeSamples; ++t){
				for(size_t rb=0; rb<downRBs; ++rb){
					pl = pathloss(neighbourPositions[bsIds],msPositions[bsId][msIds]);
					coeffDownTable[msIds][bsIds][t][rb] = pl * exponential(expMean);
				}
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
			for(size_t t=0; t<timeSamples; ++t){
				for(size_t rb=0; rb<upRBs; ++rb){
					pl = pathloss(msPositions[bsIds][msIds],neighbourPositions[bsId]);
					coeffUpTable[bsIds][0][msIds][t][rb] = pl * exponential(expMean);
				}
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
