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
#include "SINR_m.h"
#include "StreamInfo_m.h"
#include "StreamTransReq_m.h"
#include "MessageTypes.h"

#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

class StreamScheduler: public cSimpleModule{

	private:
		using ResAssign = std::pair<MessageDirection,int>;
		std::unordered_map<int,std::unordered_map<int,std::vector<StreamTransReq*>>> requests;
		std::set<int> scheduledStations;
		virtual void scheduleStreams();

	protected:
		simtime_t initOffset;
		simtime_t epsilon;
		int numberOfMs;
		int upRB;
		int downRB;
		simtime_t streamSchedPeriod;
		simtime_t tti;
		std::vector<StreamInfo*> infos;		
		std::vector<SINR*> sinrEstimate;
		SINR *estimateBS;
		std::unordered_map<unsigned long,std::unordered_map<int,ResAssign>> rbAssignments;
		virtual void initialize();
		virtual void handleMessage(cMessage *msg);
		virtual void handleSINREstimate(SINR *msg);
		virtual void printAssignment();
	
	public:
		~StreamScheduler() = default;
};
