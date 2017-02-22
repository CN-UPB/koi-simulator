/*
 * MsMac.h
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#pragma once

#include "Channel.h"
#include "includes.h"
#include "KoiData_m.h"
#include "Position.h"
#include "SINR_m.h"
#include "util.h"

#include <itpp/itbase.h>

#include <algorithm>
#include <fstream>
#include <functional>
#include <list>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace itpp;
using namespace std;
using TTISchedule = std::pair<int,std::vector<int>>;
using ScheduleList = std::vector<TTISchedule>;

class MsMac : public omnetpp::cSimpleModule  {
    private:
			unordered_map<unsigned long,list<KoiData*>> streamQueues;
			std::function<bool(const KoiData*, const KoiData*)> comparator;
			ScheduleList staticSchedule;
			ScheduleList::iterator staticIter;
			Channel *chn;
			SINR *longTermEst;
			int staticSchedLength;
			Position msPosition;
			std::ofstream* rateFile;
			int msId;
			int bsId;
			int numberOfMobileStations;
			int positionResendInterval;
			int packetLength;
			omnetpp::simtime_t initOffset;
			omnetpp::simtime_t epsilon;
			omnetpp::simtime_t tti;
			double radius;
			vector<double> velocity;
			Position initBsPos; //just used for position centering in Tkenv ms pos calc in init
			double transmissionPower;
			inline omnetpp::simtime_t positionResendTime();
			Position initMsPositionLinear();
			Position initMsPositionRand();
			Position initMsPositionLine();
			void scheduleStatic(omnetpp::cMessage* msg);
			/**
			 * @enum Placement
			 * All possible methods to determine initial Mobile Station placement
			 */
			enum Placement: int{uniformRand,params,linear,line};

			/**
			 * @var MsMac::Placement MsMac::uniformRand
			 * Place mobile stations uniformly at random in the cell.
			 */
			/**
			 * @var MsMac::Placement MsMac::params
			 * Place mobile stations as set in the initMsXPos,initMsYPos params.
			 */
			/**
			 * @var MsMac::Placement MsMac::linear
			 * Place mobile stations linearly along a "road".
			 */
			/**
			 * @var MsMac::Placement MsMac::line
			 * Place mobile stations in a straight line to the right of their BS.
			 *
			 * The grater the MS ID, the greater the distance.
			 */
				

    protected:
        virtual void initialize();
        virtual void handleMessage(omnetpp::cMessage *msg);
				virtual void finish();

    public:
        ~MsMac();
};
