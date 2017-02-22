/*
 * BsChannel.h
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 * Last edited: Jun 25, 2014
 *      Author Thomas Prinz
 */

#pragma once

#include "BsMsPositions_m.h"
#include "Channel.h"
#include "Coding.h"
#include "includes.h"
#include "NeighbourIdMatching.h"
#include "Position.h"
#include "TransInfo_m.h"
#include <itpp/itbase.h>
#include <forward_list>
#include <ostream>
#include <vector>

class BsChannel : public omnetpp::cSimpleModule  {
    private:
				Coding coding;
        int maxNumberOfNeighbours;
        int numberOfMobileStations;
				int numMSAntenna;
				int numBSAntenna;
        int upResBlocks;
        int downResBlocks;
        int bsId;
        int init_counter;
        Position bsPosition;
				omnetpp::simtime_t tti;
				omnetpp::simtime_t epsilon;
				omnetpp::simtime_t initOffset;
        std::vector<std::vector<Position>> msPositions;
        std::map <int,Position> neighbourPositions;
        Channel* channel;

    protected:
        void initialize() override;
        void handleMessage(omnetpp::cMessage *msg) override;

    public:
        ~BsChannel() override;
};
