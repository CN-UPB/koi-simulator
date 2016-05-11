/*
 * MsChannel.h
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#pragma once

#include "includes.h"
#include "Position.h"
#include "NeighbourIdMatching.h"
#include "Channel.h"
#include "TransInfoBs_m.h"
#include <vector>
#include <forward_list>


class MsChannel : public cSimpleModule  {
    private:
        int maxNumberOfNeighbours;
        Position *bsPositions;
        bool useSimpleChannelCalc;
        int simpleChannelCalcNops;
        int currentChannel;
        int downResourceBlocks;
        simtime_t epsilon;
        simtime_t tti;
        double packetLoss;
        Position msPosition;
        int bsId;
        int msId;
        vec eesm_beta_values;
        NeighbourIdMatching *neighbourIdMatching; //map the bsId to the pos in the data structures
        Channel* channel;
	/**
	 * Holds transmission information for the current TTI
	 *
	 * The vector holds lists of TransInfoBs messages, sorted by the resource 
	 * block the transmission uses.
	 */
	std::vector<std::forward_list<TransInfoBs*>> transInfos;

    protected:
        virtual void initialize();
        virtual void handleMessage(cMessage *msg);
        virtual simtime_t getProcessingDelay(cMessage *msg);

    public:
        ~MsChannel();
};
