/*
 * bestCQIScheduler.h
 *
 *  Created on: Jul 10, 2014
 *      Author: Thomas Prinz
 */
 
#pragma once
 
#include <itpp/itbase.h>
#include <omnetpp.h>
#include <vector>
#include "Schedule_m.h"
#include "Scheduler.h"

using namespace itpp;

 class bestCQIScheduler : public Scheduler{
	 protected:
		int numberOfRessourceBlocks;
		int numberOfMobileStations;
		mat sinr_;
		
	 public:
		void init(cModule *module);
		void setSINR(mat sinr);
		Schedule* calcUpSchedule(cQueue *packetQueue, std::vector<int> *transmitRequests);
		std::vector<int> calcDownSchedule(cQueue *packetQueue);
		~bestCQIScheduler();
 };
