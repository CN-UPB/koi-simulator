/**
 * @file FactoryChannel.cc
 * @brief Implementation of the FactoryChannel class.
 */

#include "FactoryChannel.h"
#include "includes.h"
#include "Position.h"
#include "util.h"
#include <fstream>
#include <cmath>
#include <vector>

using namespace omnetpp;
using std::vector;

bool FactoryChannel::init(cSimpleModule* module,
		const std::vector<std::vector<Position>>& msPositions, 
		const std::map<int,Position>& neighbourPositions){
	// First, execute the parent init method of the Channel class
	Channel::init(module,msPositions,neighbourPositions);
	plExp = module->par("plExp");
	pl0 = module->par("pl0");
	d0 = module->par("d0");
	kMean = module->par("kMean");
	kSigma = module->par("kSigma");
	shSigma = module->par("shadowSigma");
	expMean = module->par("expMean");
	// We need the gains in linear scale
	bsGain = module->par("bsGain");
	bsGain = std::pow(10,bsGain/10.0);
	msGain = module->par("msGain");
	msGain = std::pow(10,msGain/10.0);
	transPower = module->par("transmissionPower");
	initOffset = module->par("initOffset");
	randEng = boost::random::mt19937(module->getRNG(0)->intRand());
	distExp = boost::random::exponential_distribution<double>(expMean);
	distNorm = boost::random::normal_distribution<double>(0.0,shSigma);

	if(debug){
		std::string fname("coeff_table_down-"+std::to_string(bsId));
		downValues = getResultFile(fname);
		downValues << "TTI\t" << "BS\t" << "MS\t" << "RB\t" << "PL\t" 
			<< "Fade\t" << "Coeff" << "\n"; 
		fname = "coeff_table_up-"+std::to_string(bsId);
		upValues = getResultFile(fname);
		upValues << "TTI\t" << "Cell\t" << "MS\t" << "BS\t" << "RB\t" << "PL\t" 
			<< "Fade\t" << "Coeff" << "\n"; 
	}
	generateShadowing(msPositions);
	recomputeCoefficients(msPositions);
	return true;
}

double FactoryChannel::pathgain(Position sender, Position receiver){
	double distX = pow(sender.x-receiver.x,2);
	double distY = pow(sender.y-receiver.y,2);
	double dist = sqrt(distX+distY);
	double pl = (pl0+10*plExp*log10(dist/d0));
	// Convert pathloss to linear scale and return gain instead of loss
	return std::pow(10,-pl/10);
}

void FactoryChannel::generateShadowing(
				const std::vector<std::vector<Position>>& msPositions){
	int numBs = msPositions.size();
	// Generate DOWN Shadowing
	shDown.resize(numberOfMobileStations,vector<double>(numBs));
	for(int msIds=0; msIds<numberOfMobileStations; ++msIds){
		for(int bsIds=0; bsIds<numBs; ++bsIds){
			shDown[msIds][bsIds] = shadowing();
		}
	}
	// Compute UP shadowing
	shUp.resize(numBs,vector<double>());
	int numMs;
	for(int bsIds=0; bsIds<numBs; ++bsIds){
		numMs = msPositions[bsIds].size();
		shUp[bsIds].resize(numMs);
		for(int msIds=0; msIds<numMs; ++msIds){
			shUp[bsIds][msIds] = shadowing();
		}
	}
	if(d2dActive){
		// Compute D2D UP shadowing
		shD2DUp.resize(numBs,
				VectorNd<double,2>(numberOfMobileStations));
		for(int bsIds=0; bsIds<numBs; ++bsIds){
			numMs = msPositions[bsIds].size();
			for(int recMsId=0; recMsId<numberOfMobileStations; ++recMsId){
				shD2DUp[bsIds][recMsId].resize(numMs);
				for(int msIds=0; msIds<numMs; ++msIds){
					shD2DUp[bsIds][recMsId][msIds] = shadowing();
				}
			}
		}
		// Compute D2D DOWN shadowing
		shD2DDown.resize(numBs,
				VectorNd<double,2>(numberOfMobileStations));
		for(int bsIds=0; bsIds<numBs; ++bsIds){
			numMs = msPositions[bsIds].size();
			for(int recMsId=0; recMsId<numberOfMobileStations; ++recMsId){
				shD2DDown[bsIds][recMsId].resize(numMs);
				for(int msIds=0; msIds<numMs; ++msIds){
					shD2DDown[bsIds][recMsId][msIds] = shadowing();
				}
			}
		}
	}
}

void FactoryChannel::recomputeCoefficients(
		const vector<vector<Position>>& msPositions){
	auto val = ((simTime().dbl()-initOffset)/tti).dbl();
	// +1 because when this method is called, it computes the values for the 
	// NEXT TTI, not the current one.
	int tti = std::floor(val);
	int numBs = msPositions.size();
	// Compute DOWN RB coefficients
	coeffDownTable.resize(numberOfMobileStations,
			VectorNd<double,2>(numBs,
					vector<double>(upRBs)));
	double pg = 0.0;
	double fading = 0.0;
	for(int msIds=0; msIds<numberOfMobileStations; ++msIds){
		for(int bsIds=0; bsIds<numBs; ++bsIds){
			pg = pathgain(neighbourPositions[bsIds],msPositions[bsId][msIds]);
			for(size_t rb=0; rb<downRBs; ++rb){
				fading = fadingExponential();
				coeffDownTable[msIds][bsIds][rb] = pg * fading * shDown[msIds][bsIds];
				if(simTime()>initOffset && debug){
					downValues << tti << "\t" << bsIds << "\t"
						<< msIds << "\t"
						<< rb << "\t"
						<< pg << "\t"
						<< fading << "\t"
						<< coeffDownTable[msIds][bsIds][rb]
						<< std::endl;
				}
			}
		}
	}
	// Compute UP RB coefficients
	coeffUpTable.resize(numBs,
			VectorNd<double,3>(1));
	int numMs;
	for(int bsIds=0; bsIds<numBs; ++bsIds){
		numMs = msPositions[bsIds].size();
		coeffUpTable[bsIds][0].resize(numMs,vector<double>(upRBs));
		for(int msIds=0; msIds<numMs; ++msIds){
			pg = pathgain(msPositions[bsIds][msIds],neighbourPositions[bsId]);
			for(size_t rb=0; rb<upRBs; ++rb){
				fading = fadingExponential();
				coeffUpTable[bsIds][0][msIds][rb] = pg * fading * shUp[bsIds][msIds];
				if(simTime()>initOffset && debug){
					upValues << tti << "\t" << bsIds << "\t"
						<< msIds << "\t"
						<< bsId << "\t"
						<< rb << "\t"
						<< pg << "\t"
						<< fading << "\t"
						<< coeffUpTable[bsIds][0][msIds][rb]
						<< std::endl;
				}
			}
		}
	}
	if(d2dActive){
		// Compute D2D UP RB coefficients
		coeffUpD2DTable.resize(numBs,
				VectorNd<double,3>(numberOfMobileStations));
		for(int bsIds=0; bsIds<numBs; ++bsIds){
			numMs = msPositions[bsIds].size();
			for(int recMsId=0; recMsId<numberOfMobileStations; ++recMsId){
				coeffUpD2DTable[bsIds][recMsId].resize(numMs,vector<double>(upRBs));
				for(int msIds=0; msIds<numMs; ++msIds){
					pg = pathgain(msPositions[bsIds][msIds],msPositions[bsId][recMsId]);
					for(size_t rb=0; rb<upRBs; ++rb){
						fading = fadingExponential();
						coeffUpD2DTable[bsIds][recMsId][msIds][rb] = pg * fading * shD2DUp[bsIds][recMsId][msIds];
					}
				}
			}
		}
		// Compute D2D DOWN RB coefficients
		coeffDownD2DTable.resize(numBs,
				VectorNd<double,3>(numberOfMobileStations));
		for(int bsIds=0; bsIds<numBs; ++bsIds){
			numMs = msPositions[bsIds].size();
			for(int recMsId=0; recMsId<numberOfMobileStations; ++recMsId){
				coeffDownD2DTable[bsIds][recMsId].resize(numMs,vector<double>(upRBs));
				for(int msIds=0; msIds<numMs; ++msIds){
					pg = pathgain(msPositions[bsIds][msIds],msPositions[bsId][recMsId]);
					for(size_t rb=0; rb<downRBs; ++rb){
						fading = fadingExponential();
						coeffDownD2DTable[bsIds][recMsId][msIds][rb] = pg * fading * shD2DDown[bsIds][recMsId][msIds];
					}
				}
			}
		}
	}
}

void FactoryChannel::updateChannel(const vector<vector<Position>>& msPos){
	recomputeCoefficients(msPos);
}

void FactoryChannel::handleMessage(cMessage* msg){
}

void FactoryChannel::recomputePerTTIValues(){
	recomputeCoefficients(msPos);
}

double FactoryChannel::fadingExponential(){
	// Draw 16 values and take the best one
	double choosen = distExp(randEng);
	double curr;
	for(int i = 0; i<15; ++i){
		curr = distExp(randEng);
		if(curr>choosen){
			choosen = curr;
		}
	}
	return choosen;
}

double FactoryChannel::shadowing(){
	double inDB = distNorm(randEng);
	return std::pow(10.0,inDB/10);
}

double FactoryChannel::calcLongtermUpSINR(int rb, int msId, double transPower){
	double received = transPower 
		* pathgain(msPos[bsId][msId],neighbourPositions[bsId]) * shUp[bsId][msId];
	double interference = 0.0;
	interference += getTermalNoise(300,chnBandwidth);
	return 10*log10(received/interference);
}

double FactoryChannel::calcLongtermDownSINR(int rb, int msId, double transPower){
	double received = transPower 
		* pathgain(msPos[bsId][msId],neighbourPositions[bsId]) * shDown[msId][bsId];
	double interference = 0.0;
	interference += getTermalNoise(300,chnBandwidth);
	return 10*log10(received/interference);
}

FactoryChannel::~FactoryChannel(){
	downValues.close();
	upValues.close();
}
