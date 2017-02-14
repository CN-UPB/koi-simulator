/*
 * BsMac.h
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#pragma once

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

class BsMac : public omnetpp::cSimpleModule  {
    private:
				std::function<bool(const KoiData*, const KoiData*)> comparator;
        int numberOfMobileStations;
        int maxNumberOfNeighbours;
        int bsId;
				int sinrEstCount;
        Position pos;
        double initOffset;
        vector<bool> msPosUpdateArrived;
				std::vector<Position> msPositions;
        NeighbourIdMatching *neighbourIdMatching; //matching for bsId <-> dataStrPos
				double transmissionPower;
				omnetpp::simtime_t tti;
				omnetpp::simtime_t epsilon;
				unordered_map<unsigned long,std::list<KoiData*>> streamQueues;
				ofstream delays_file;
				ofstream rate_file;
				void writePositions();

    protected:
        void initialize() override;
				void finish() override;
        void handleMessage(omnetpp::cMessage *msg) override;
        
        //important: these methods doesn't delete the msg!
        virtual void sendToNeighbourCells(omnetpp::cMessage *msg);
        virtual void sendDelayedToNeighbourCells(omnetpp::cMessage *msg, 
						const omnetpp::simtime_t& delay);

    public:
        ~BsMac() override;
};
