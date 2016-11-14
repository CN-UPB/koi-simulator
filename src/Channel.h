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
#include <forward_list>
#include <utility>
#include <vector>

using std::vector;

class Channel{
	protected:
		/**
		 * Channel bandwidth in Hz
		 */
		double chnBandwidth;
		/**
		 * Unique ID of the BS in the cell this channel serves
		 */
		int bsId;
		/**
		 * @brief Table to save downlink coefficients
		 *
		 * [receivingMsId][sendingBsId][counter][RB]
		 */
		vector<vector<vector<vector<double>>>> coeffDownTable;
		/**
		 * @brief Table to save uplink coefficients
		 *
		 * [sendingMsCellID][0][sendingMsId][counter][RB]
		 */
		vector<vector<vector<vector<vector<double>>>>> coeffUpTable;
		/**
		 * @brief Table to save D2D DOWN Rb coefficients
		 *
		 * [sendingMsCellID][receivingMsId][sendingMsId][counter][RB]
		 */
		vector<vector<vector<vector<vector<double>>>>> coeffDownD2DTable;
		/**
		 * @brief Table to save D2D UP Rb coefficients
		 *
		 * [sendingMsCellID][receivingMsId][sendingMsId][counter][RB]
		 */
		vector<vector<vector<vector<vector<double>>>>> coeffUpD2DTable;
		/**
		 * Should interference be considered or not
		 */
		bool considerInterference;
		/**
		 * Be more verbose.
		 *
		 * WARNING! This may potentially produce Gigabytes of additional data!
		 */
		bool debug;
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
		 * Bandwidth in Hz per resource block
		 */
		double rbBandwidth;
		/**
		 * If position resend intervall > 1, it counts the current TTI
		 */
		int SINRcounter;
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
		 * Number of TTIs until Position is updated 
		 *
		 * (Number of Time Samples for Channel Model)
		 */
		int timeSamples;
		/**
		 * Interference information from local and neighbouring senders
		 */
		std::pair<vector<std::forward_list<TransInfo*>>,
			vector<std::forward_list<TransInfo*>>> transInfo;
		/**
		 * Transmission Time Interval
		 */
		double tti;
		/**
		 * Number of up resource blocks
		 */
		int upRBs;
		/**
		 * @brief Calculate interference for transmission
		 */
		virtual double calcInterference(std::forward_list<TransInfo*>& interferers,
				int rb,
				int receiverId,
				int SINRCounter,
				MessageDirection dir);
		/**
		 * @brief Computes the Termal Noise
		 *
		 * The default implementation is the Johnson-Nyquist noise
		 */
		double getTermalNoise(double temp, double bandwidth);

	public:
		/**
		 * @brief Default channel constructor
		 */
		Channel(){
			bsId = -1;
			initialized = false;
		}

		virtual void addTransInfo(TransInfo* trans);

		virtual double calcUpSINR(int RB, 
				int msId,
				double transPower);

		virtual double calcDownSINR(int RB, 
				int msId,
				double transPower);

		virtual double calcD2DSINR(int RB, 
				int sendMsID,
				int receiveMsId,
				MessageDirection direction,
				double transPower);

		virtual double calcAvgUpSINR(int RB, 
				int msId,
				double transPower);

		virtual double calcAvgDownSINR(int RB, 
				double transPower);

		virtual double calcAvgD2DDownSINR(int RB, 
				int msId,
				double transPower);

		virtual void clearTransInfo();

		/**
		 * @brief It may be necessary for the Channel to receive Messages
		 * 
		 */
		virtual void handleMessage(cMessage* msg) = 0;

		/**
		 * @brief Initialize the channel model with the given positions
		 */
		virtual bool init(cSimpleModule* module,
				const vector<vector<Position>>& msPositions, 
				std::map<int,Position>& neighbourPositions)=0;

		/**
		 * @brief Output the Up coefficient table to out stream
		 */
		virtual std::ostream& printCoeffUpTables(std::ostream& out);

		/**
		 * @brief Output the Down coefficient table to out stream
		 */
		virtual std::ostream& printCoeffDownTables(std::ostream& out);

		/**
		 * @brief Updates the Channel if necessary for moving MS
		 */
		virtual void updateChannel(const vector<vector<Position>>& msPos) = 0;
		        
		virtual ~Channel();
};
