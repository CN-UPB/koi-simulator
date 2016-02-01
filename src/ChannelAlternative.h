/*
 * ChannelAlternative.h
 *
 *  Created on: Jul 22, 2014
 *      Author: Thomas Prinz
 * 
 */
 
#pragma once

#include <omnetpp.h>
#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>
#include "Position.h"
#include "NeighbourIdMatching.h"
#include "VisibilityRegion.h"
#include "cluster.h"
#include "VisibilityRegionMessage_m.h"
#include "ClusterMessage_m.h"
#include "Channel.h"
#include "PointerExchange_m.h"

using namespace std;
using namespace itpp;

class ChannelAlternative : public Channel{
	private:
		double tenlogk;
		double alpha;
		int FADING_PATHS;
		int nrMSs;
		double stationaryDoppler;
		double DELAY_RMS;
		int RANDOM_SEED_ID;
		int myBsID;
		int myID;
		double*** angle_of_arrival;
		double*** delay;
		double* bsXposition;
		double* bsYposition;
		double* msXposition;
		double* msYposition;
		double* shadowLoss;
		vector<double> velocity;
		NeighbourIdMatching *neighbourIdMatching;

	public:
		// Initialize your Channel through ini access via module pointer.
		bool init(cSimpleModule* module, Position** msPositions, std::map <int,Position> neighbourPositions);
		// It may be necessary for the Channel to receive Message from other LPs.
		void handleMessage(cMessage* msg);
		// Computes the pathloss for a given distance using an arbitrary model.
		double calcPathloss(double dist);
		// Computes the Termal Noise.
		double getTermalNoise(double temp, double bandwidth);
		// Calculates the current SINR for given interferers and a given RB.
		double calcSINR(int RB, vector<double> &power, vector<Position> &pos, vector<int> &bsId_, bool up, int msId);
		// Calculates the current SINR for given interferers and a given RB.
		vec calcSINR(vector<double> &power, vector<Position> &pos, vector<int> &bsId_, bool up, int msId);
		// Updates the Channel if necessary for moving MS
		void updateChannel(Position** msPos);
		// Calculates Jackes like Fading
		double calculateFading(int msId, double frequency, double mobile_speed, int intbsID);
		double calculateLoss(int, double, double, int);
		double toDB(double val);
		double toLinear(double val);
		//Destructor
		~ChannelAlternative();
};
