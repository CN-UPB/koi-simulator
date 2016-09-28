/*
 * Channel.h
 *
 *  Created on: Jul 15, 2014
 *      Author: Thomas Prinz
 * 
 */
 
#pragma once

#include "includes.h"
#include "NeighbourIdMatching.h"
#include "Position.h"
#include "TransInfo_m.h"
#include <itpp/itbase.h>
#include <forward_list>
#include <utility>
#include <vector>

class Channel{
	protected:
		/**
		 * Unique ID of the BS in the cell this channel serves
		 */
		int bsId;
		/**
		 * Number of down resource blocks
		 */
		int downRBs;
		/**
		 * True iff METISChannel::init has been called
		 */
		bool initialized;
		/**
		 * Pointer to OMNeT module for intermodule communication
		 */
		cSimpleModule *initModule;
		/**
		 * Max Number of Neighbour BS
		 */
		int maxNumberOfNeighbours;
		/**
		 * Information about neighbouring cells
		 */
		NeighbourIdMatching *neighbourIdMatching;
		/**
		 * Positions of Neighbour BS
		 */
		std::map<int,Position> neighbourPositions;
		/**
		 * Number of MSs within the local cell
		 */
		int numberOfMobileStations;
		/**
		 * Size of playground in the X dimension
		 */
		double sizeX;
		/**
		 * Size of playground in the Y dimension
		 */
		double sizeY;
		/**
		 * The speed of light in m/s
		 */
		const static double speedOfLight;
		/**
		 * Interference information from local and neighbouring senders
		 */
		std::pair<std::vector<std::forward_list<TransInfo*>>,
			std::vector<std::forward_list<TransInfo*>>> transInfo;
		/**
		 * Transmission Time Interval
		 */
		double tti;
		/**
		 * Number of up resource blocks
		 */
		int upRBs;

	public:
		/**
		 * @brief Default channel constructor
		 */
		Channel(){
			bsId = -1;
			initialized = false;
		}
		// Initialize your Channel through ini access via module pointer.
		// The MS/BS Positions may not be needed for every channel.
		virtual bool init(cSimpleModule* module,
				const std::vector<std::vector<Position>>& msPositions, 
				std::map<int,Position>& neighbourPositions)=0;
		// It may be necessary for the Channel to receive Message from other LPs.
		// Note: If you send a msg with kind x, it is received with kind (x+1) by neighbour cells
		// When it has the name "CHANNEL_INFO"
		virtual void handleMessage(cMessage* msg) = 0;
		// Computes the pathloss for a given distance using an arbitrary model.
		virtual double calcPathloss(double dist) = 0;
		// Computes the Termal Noise.
		virtual double getTermalNoise(double temp, double bandwidth) = 0;
		// Updates the Channel if necessary for moving MS
		virtual void updateChannel(const std::vector<std::vector<Position>>& msPos)
			= 0;
		        
		virtual double calcUpSINR(int RB, 
				int msId,
				double transPower) = 0;

		virtual double calcDownSINR(int RB, 
				int msId,
				double transPower) = 0;

		virtual double calcD2DSINR(int RB, 
				int sendMsID,
				int receiveMsId,
				MessageDirection direction,
				double transPower) = 0;

		virtual double calcAvgUpSINR(int RB, 
				int msId,
				double transPower) = 0;

		virtual double calcAvgDownSINR(int RB, 
				double transPower) = 0;

		virtual double calcAvgD2DDownSINR(int RB, 
				int msId,
				double transPower) = 0;

		virtual void addTransInfo(TransInfo* trans);

		virtual ~Channel(){}
};
