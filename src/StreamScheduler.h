/**
 * \class StreamScheduler StreamScheduler.h
 *
 * The Stream Scheduler for assigning packet streams to ressource blocks.
 *
 * This module is responsible for assigning streams to ressource blocks based 
 * on the specific streams' packet interarrival time. This scheduling is 
 * repeated periodically.
 */

#pragma once

#include "includes.h"
#include "RBScheduler.h"
#include "StreamInfo_m.h"
#include "StreamTransReq_m.h"

#include <unordered_map>
#include <vector>

class StreamScheduler: public cSimpleModule{

	private:
		simtime_t initOffset;
		int numberOfMs;
		int resourceBlocks;
		simtime_t streamSchedPeriod;
		simtime_t tti;
		std::vector<StreamInfo*> infos;		
		std::unordered_map<int,std::vector<StreamTransReq*>> requests;
		std::unordered_map<int,std::unordered_map<int,int>> rbAssignments;
		virtual void scheduleStreams();

	protected:
		virtual void initialize();
		virtual void handleMessage(cMessage *msg);
	
	public:
		~StreamScheduler() = default;
};
