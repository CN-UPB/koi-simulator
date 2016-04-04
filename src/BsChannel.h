/*
 * BsChannel.h
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 * Last edited: Jun 25, 2014
 *      Author Thomas Prinz
 */

#pragma once

#include "includes.h"
#include "Position.h"
#include "BsMsPositions_m.h"
#include "NeighbourIdMatching.h"
#include "VisibilityRegion.h"
#include "cluster.h"
#include "Channel.h"
#include <vector>

class BsChannel : public cSimpleModule  {
    private:
        int maxNumberOfNeighbours;
        int numberOfMobileStations;
        int upResBlocks;
        int **schedules;
        int *scheduleDirection;
        double **schedulePower;
        double *maxPower;
        bool useSimpleChannelCalc;
        int simpleChannelCalcNops;
        int bsId;
        int init_counter;
        int numberOfVR;
        Position bsPosition;
        double packetLoss;
        bool scheduleCatch;
        simtime_t tti;
        simtime_t epsilon;
        Position **msPositions;
        NeighbourIdMatching *neighbourIdMatching;          // map the bsId to the pos in the data structures
        std::map <int,Position> neighbourPositions;
        std::vector<VisibilityRegion> VR;
        VisibilityRegion LOS_VR;
        std::vector<VisibilityRegion> ForeignVR;
        Cluster localcluster_bs;
        Cluster localcluster_ms;
        std::vector<Cluster> remoteCluster;
        Channel* channel;
        vec eesm_beta_values;

    protected:
        virtual void initialize();
        virtual void handleMessage(cMessage *msg);
        virtual simtime_t getProcessingDelay(cMessage *msg);
        NeighbourIdMatching* getNeighbourMatching(){return neighbourIdMatching;}

    public:
        virtual ~BsChannel();
};
