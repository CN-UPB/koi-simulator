/*
 * Scheduler.h
 *
 *  Created on: Jul 10, 2014
 *      Author: Thomas Prinz
 */
 
#pragma once
 
#include <itpp/itbase.h>
#include "includes.h"
#include <vector>
#include "Schedule_m.h"

using namespace itpp;

class Scheduler{
 protected:
 
 public:
	virtual void init(cModule *module) = 0;
	virtual void setSINR(mat sinr) = 0;
	virtual void updateSINR() = 0;
	virtual Schedule* calcUpSchedule(cQueue *packetQueue, std::vector<int> *transmitRequests) = 0;
	virtual std::vector<int> calcDownSchedule(cQueue *packetQueue) = 0;
	virtual ~Scheduler() = 0;
};
