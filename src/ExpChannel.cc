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
	expMean = par("expMean");
	recomputeCoefficients(msPositions);
	return true;
}

void ExpChannel::recomputeCoefficients(
		const vector<vector<Position>>& msPositions){
	size_t numBs = msPositions.size();
	// Compute DOWN RB coefficients
	coeffDownTable.resize(numberOfMobileStations,
			vector<vector<vector<double>>>(numBs,
				vector<vector<double>>(timeSamples,
					vector<double>(upRBs))));
	for(size_t msIds=0; msIds<numberOfMobileStations; ++msIds){
		for(size_t bsIds=0; bsIds<numBs; ++bsIds){
			for(size_t t=0; t<timeSamples; ++t){
				for(size_t rb=0; rb<downRBs; ++rb){
					coeffDownTable[msIds][bsIds][t][rb] = exponential(expMean);
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
				vector<vector<double>>>(timeSamples,
					vector<double>(upRBs)));
		for(size_t msIds=0; msIds<numMs; ++msIds){
			for(size_t t=0; t<timeSamples; ++t){
				for(size_t rb=0; rb<upRBs){
					coeffUpTable[bsIds][0][msIds][t][rb] = exponential(expMean);
				}
			}
		}
	}
	// Compute D2D UP RB coefficients
	coeffUpD2DTable.resize(numBs,
			vector<vector<vector<vector<double>>>>(numberOfMobileStations));
	size_t numMs;
	for(size_t bsIds=0; bsIds<numBs; ++bsIds){
		numMs = msPositions[bsIds].size();
		for(size_t recMsId=0; recMsId<numberOfMobileStations; ++recMsId){
			coeffUpD2DTable[bsIds][recMsId].resize(numMs,
					vector<vector<double>>>(timeSamples,
						vector<double>(upRBs)));
			for(size_t msIds=0; msIds<numMs; ++msIds){
				for(size_t t=0; t<timeSamples; ++t){
					for(size_t rb=0; rb<upRBs){
						coeffUpD2DTable[bsIds][recMsId][msIds][t][rb] = exponential(expMean);
					}
				}
			}
		}
	}
	// Compute D2D DOWN RB coefficients
	coeffDownD2DTable.resize(numBs,
			vector<vector<vector<vector<double>>>>(numberOfMobileStations));
	size_t numMs;
	for(size_t bsIds=0; bsIds<numBs; ++bsIds){
		numMs = msPositions[bsIds].size();
		for(size_t recMsId=0; recMsId<numberOfMobileStations; ++recMsId){
			coeffDownD2DTable[bsIds][recMsId].resize(numMs,
					vector<vector<double>>>(timeSamples,
						vector<double>(upRBs)));
			for(size_t msIds=0; msIds<numMs; ++msIds){
				for(size_t t=0; t<timeSamples; ++t){
					for(size_t rb=0; rb<downRBs){
						coeffDownD2DTable[bsIds][recMsId][msIds][t][rb] = exponential(expMean);
					}
				}
			}
		}
	}
}
