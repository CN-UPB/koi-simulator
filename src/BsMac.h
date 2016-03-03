/*
 * BsMac.h
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#pragma once

#include <itpp/itbase.h>
#include "includes.h"
#include "Position.h"
#include "util.h"
#include <vector>
#include "NeighbourIdMatching.h"
#include "randomScheduler.h"
#include "roundRobinScheduler.h"
#include "bestCQIScheduler.h"
#include "proportionalFairScheduler.h"

using namespace std;
using namespace itpp;

enum DATADIRECTION {UP, DOWN, SPECIAL};

class BsMac : public cSimpleModule  {
    private:
        int numberOfMobileStations;
        int maxNumberOfNeighbours;
        int bsId;
        int currentChannel;
        int upResBlocks;
        int downResBlocks;
        int sinr_est;
        vector<int> *downSchedule;
        Position pos;
        double initOffset;
        vector<bool> *msPosUpdateArrived;
        vector<int> *transmitRequests;
        Position *msPositions;
        mat SINR_;
        NeighbourIdMatching *neighbourIdMatching; //matching for bsId <-> dataStrPos
        int ownDataStrId; //pos in the dataStr for bsId
        int downToUpPeriodicity; // The periodicity for change in Up and Downlink
        int currentPeriodicity; // The current periodicity
        Scheduler *scheduler;
        vec eesm_beta_values;
        mat blerTable;

        simtime_t tti;
        simtime_t epsilon;
        cQueue *packetQueue;

    protected:
        virtual void initialize();
        virtual void handleMessage(cMessage *msg);
        virtual simtime_t getProcessingDelay(cMessage *msg);
        virtual void calcNewSchedule();
        virtual DATADIRECTION getDataDirection();
        inline virtual simtime_t scheduleTime();
        
        double getEffectiveSINR(vector<double> sinrValues);
        double getBler(int cqi, double sinr);

        //important: these methods doesn't delete the msg!
        virtual void sendToNeighbourCells(cMessage *msg);
        virtual void sendDelayedToNeighbourCells(cMessage *msg, simtime_t delay);

    public:
        BsMac();
        ~BsMac();
};
