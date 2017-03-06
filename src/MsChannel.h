/*
 * MsChannel.h
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#pragma once

#include "Channel.h"
#include "Coding.h"
#include "includes.h"
#include "NeighbourIdMatching.h"
#include "Position.h"
#include "TransInfo_m.h"
#include <itpp/itbase.h>
#include <fstream>
#include <vector>
#include <forward_list>

class MsChannel : public omnetpp::cSimpleModule  {
    private:
				Coding coding;
        int maxNumberOfNeighbours;
        Position *bsPositions;
        int downResourceBlocks;
        int upResourceBlocks;
				int numBSAntenna;
				int numMSAntenna;
				std::ofstream sinrFile;
				omnetpp::simtime_t epsilon;
				omnetpp::simtime_t tti;
				omnetpp::simtime_t initOffset;
        int bsId;
        int msId;
				itpp::vec eesm_beta_values;
        Channel* channel;
				bool debug;
				bool d2dActive;
				std::ofstream* mcs_file;

    protected:
        virtual void initialize();
				virtual void finish();
        virtual void handleMessage(omnetpp::cMessage *msg);

    public:
        ~MsChannel();
};
