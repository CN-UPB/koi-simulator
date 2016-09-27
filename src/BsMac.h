/*
 * BsMac.h
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#pragma once

#include <itpp/itbase.h>
#include "includes.h"
#include "KoiData_m.h"
#include "Position.h"
#include "util.h"
#include "NeighbourIdMatching.h"
#include <vector>
#include <unordered_map>
#include <list>

using namespace std;
using namespace itpp;

class BsMac : public cSimpleModule  {
    private:
        int numberOfMobileStations;
        int maxNumberOfNeighbours;
        int bsId;
        int currentChannel;
				int resourceBlocks;
        int sinr_est;
        Position pos;
        double initOffset;
        vector<bool> *msPosUpdateArrived;
        Position *msPositions;
        vector<double> sinrDown;
        NeighbourIdMatching *neighbourIdMatching; //matching for bsId <-> dataStrPos
        int ownDataStrId; //pos in the dataStr for bsId
        vec eesm_beta_values;
        mat blerTable;
				double transmissionPower;
        simtime_t tti;
        simtime_t epsilon;
				simsignal_t avgRatePerStation;
				unordered_map<unsigned long,std::list<KoiData*>> streamQueues;
				void writePositions();

    protected:
        virtual void initialize();
        virtual void handleMessage(cMessage *msg);
        
        double getEffectiveSINR(vector<double> sinrValues);
        double getBler(int cqi, double sinr);

        //important: these methods doesn't delete the msg!
        virtual void sendToNeighbourCells(cMessage *msg);
        virtual void sendDelayedToNeighbourCells(cMessage *msg, simtime_t delay);

    public:
        ~BsMac();
};
