/*
 * MsChannel.h
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#pragma once

#include "Channel.h"
#include "includes.h"
#include "NeighbourIdMatching.h"
#include "Position.h"
#include "TransInfo_m.h"
#include <itpp/itbase.h>
#include <fstream>
#include <vector>
#include <forward_list>

class MsChannel : public cSimpleModule  {
    private:
        int maxNumberOfNeighbours;
        Position *bsPositions;
        int downResourceBlocks;
        int upResourceBlocks;
				std::ofstream sinrFile;
        simtime_t epsilon;
        simtime_t tti;
        simtime_t initOffset;
        int bsId;
        int msId;
				itpp::vec eesm_beta_values;
        Channel* channel;

    protected:
        virtual void initialize();
				virtual void finish();
        virtual void handleMessage(cMessage *msg);

    public:
        ~MsChannel();
};
