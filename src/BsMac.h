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

#include <fstream>
#include <functional>
#include <list>
#include <unordered_map>
#include <vector>

using namespace std;
using namespace itpp;

class BsMac : public cSimpleModule  {
    private:
				std::function<bool(const KoiData*, const KoiData*)> comparator;
        int numberOfMobileStations;
        int maxNumberOfNeighbours;
        int bsId;
        int sinr_est;
				int sinrEstCount;
        Position pos;
        double initOffset;
        vector<bool> *msPosUpdateArrived;
        Position *msPositions;
        NeighbourIdMatching *neighbourIdMatching; //matching for bsId <-> dataStrPos
        int ownDataStrId; //pos in the dataStr for bsId
        vec eesm_beta_values;
        mat blerTable;
				double transmissionPower;
        simtime_t tti;
        simtime_t epsilon;
				unordered_map<unsigned long,std::list<KoiData*>> streamQueues;
				ofstream delays_file;
				ofstream rate_file;
				void writePositions();

    protected:
        virtual void initialize();
				virtual void finish();
        virtual void handleMessage(cMessage *msg);
        
        double getEffectiveSINR(vector<double> sinrValues);
        double getBler(int cqi, double sinr);

        //important: these methods doesn't delete the msg!
        virtual void sendToNeighbourCells(cMessage *msg);
        virtual void sendDelayedToNeighbourCells(cMessage *msg, simtime_t delay);

    public:
        ~BsMac();
};
