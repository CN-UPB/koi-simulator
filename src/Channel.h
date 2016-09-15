/*
 * Channel.h
 *
 *  Created on: Jul 15, 2014
 *      Author: Thomas Prinz
 * 
 */
 
#pragma once

#include "includes.h"
#include <itpp/itbase.h>
#include "Position.h"
#include "TransInfo_m.h"
#include <forward_list>

using namespace std;
using namespace itpp;

class Channel{
	protected:
		const static double speedOfLightVac;

	public:
		// Initialize your Channel through ini access via module pointer. The MS/BS Positions may not be needed for every channel.
		virtual bool init(cSimpleModule* module, const vector<vector<Position>>& msPositions, std::map <int,Position> neighbourPositions)=0;
		// It may be necessary for the Channel to receive Message from other LPs.
		// Note: If you send a msg with kind x, it is received with kind (x+1) by neighbour cells
		// When it has the name "CHANNEL_INFO"
		virtual void handleMessage(cMessage* msg) = 0;
		// Computes the pathloss for a given distance using an arbitrary model.
		virtual double calcPathloss(double dist) = 0;
		// Computes the Termal Noise.
		virtual double getTermalNoise(double temp, double bandwidth) = 0;
		// Updates the Channel if necessary for moving MS
		virtual void updateChannel(const vector<vector<Position>>& msPos) = 0;
		        
		virtual double calcUpSINR(int RB, 
				std::forward_list<TransInfo*> &interferers,
				int msId,
				double transPower) = 0;

		virtual double calcDownSINR(int RB, 
				std::forward_list<TransInfo*> &interferers,
				int msId,
				double transPower) = 0;

		virtual double calcD2DSINR(int RB, 
				std::forward_list<TransInfo*> &interferers,
				int sendMsID,
				int receiveMsId,
				MessageDirection direction,
				double transPower) = 0;

		virtual double calcAvgUpSINR(int RB, 
				std::forward_list<TransInfo*> &interferers,
				int msId,
				double transPower) = 0;

		virtual double calcAvgDownSINR(int RB, 
				std::forward_list<TransInfo*> &interferers,
				double transPower) = 0;

		virtual double calcAvgD2DDownSINR(int RB, 
				std::forward_list<TransInfo*> &interferers,
				int msId,
				double transPower) = 0;

		virtual ~Channel(){}
};
