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
		
		// Power values and Positions.
	    vector<Position> senderPosition;
        vector<Position> targetPosition;
        vector<vector<Position>> interfererPositions;
        vector<vector<double>> interfererPower;
        vector<double> senderPower;

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
		// Calculates the current SINR for given interferers and a given RB.
		virtual double calcSINR(int RB, vector<double> &power, vector<Position> &pos, vector<int> &bsId_, bool up, int msId) = 0;
		// Calculates the current SINR for given interferers and a given RB.
		virtual vec calcSINR(vector<double> &power, vector<Position> &pos, vector<int> &bsId_, bool up, int msId) = 0;
		// Updates the Channel if necessary for moving MS
		virtual void updateChannel(const vector<vector<Position>>& msPos) = 0;
		        
		// Store positions of sender/receiver/interferes (Id because pointer is shared)
		void setSenderPosition(Position p, double power, int Id) { senderPosition[Id] = p; senderPower[Id] = power; }

		virtual double calcUpSINR(int RB, 
				std::forward_list<TransInfo*> &interferers,
				int msId,
				double transPower){return 0.0;}

		virtual double calcDownSINR(int RB, 
				std::forward_list<TransInfo*> &interferers,
				int msId,
				double transPower){return 0.0;}

		virtual double calcD2DSINR(int RB, 
				std::forward_list<TransInfo*> &interferers,
				int sendMsID,
				int receiveMsId,
				MessageDirection direction,
				double transPower){return 0.0;}

		Position getSenderPosition(int Id) { return senderPosition.at(Id); }
		void setTargetPosition(Position p, int Id) { targetPosition[Id] = p; }
		Position getTargetPosition(int Id) { return targetPosition.at(Id); }
		void clearInterfererPostitions(int Id) { interfererPositions.at(Id).clear(); interfererPower.at(Id).clear(); }
		void addInterfererPosition(Position p, double power, int Id) { interfererPositions.at(Id).push_back(p); interfererPower.at(Id).push_back(power); }

		virtual ~Channel(){}
};
