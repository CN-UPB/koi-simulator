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
#include <vector>
#include <forward_list>
#include <ostream>

class BsChannel : public cSimpleModule  {
    private:
				Coding coding;
        int maxNumberOfNeighbours;
        int numberOfMobileStations;
				int numMSAntenna;
				int numBSAntenna;
        int upResBlocks;
        int downResBlocks;
        int **schedules;
        int *scheduleDirection;
        double **schedulePower;
        double *maxPower;
        int bsId;
        size_t init_counter;
        int numberOfVR;
        Position bsPosition;
        bool scheduleCatch;
        simtime_t tti;
        simtime_t epsilon;
				simtime_t initOffset;
        std::vector<std::vector<Position>> msPositions;
        NeighbourIdMatching *neighbourIdMatching;          // map the bsId to the pos in the data structures
        std::map <int,Position> neighbourPositions;
        Channel* channel;
				itpp::vec eesm_beta_values;

    protected:
        virtual void initialize();
        virtual void handleMessage(cMessage *msg);
        NeighbourIdMatching* getNeighbourMatching(){return neighbourIdMatching;}

    public:
        virtual ~BsChannel();
};
