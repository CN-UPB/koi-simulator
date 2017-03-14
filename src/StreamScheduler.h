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
#include "StaticSchedule_m.h"
#include "StreamInfo_m.h"
#include "StreamTransReq_m.h"
#include "MessageTypes.h"

#include <functional> 
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

class StreamScheduler: public omnetpp::cSimpleModule{

	private:
		using ResAssign = std::pair<MessageDirection,int>;
		std::unordered_map<int,std::unordered_map<int,std::vector<StreamTransReq*>>> requests;
		std::function<bool(const KoiData*,const KoiData*)> defaultPacketSorter;

	protected:
		omnetpp::simtime_t initOffset;
		omnetpp::simtime_t epsilon;
		bool debug;
		int numberOfMs;
		int upRB;
		int downRB;
		std::vector<int> assignedUpRB;
		std::vector<int> assignedDownRB;
		int staticSchedLength;
		bool upStatic;
		bool downStatic;
		omnetpp::simtime_t streamSchedPeriod;
		std::set<int> scheduledStations;
		omnetpp::simtime_t tti;
		std::vector<StreamInfo*> infos;		
		std::vector<SINR*> sinrEstimate;
		SINR *estimateBS;
		std::vector<SINR*> longtermSinrEstimate;
		SINR *longtermEstimateBS;
		std::unordered_map<unsigned long,std::unordered_map<int,ResAssign>> rbAssignments;
		void initialize() override;
		void handleMessage(omnetpp::cMessage *msg) override;
		virtual void handleSINREstimate(SINR *msg);
		virtual void printAssignment();
		virtual void scheduleDynStreams();
		virtual void distributeStaticSchedules();
		virtual std::unordered_map<int,ScheduleList> scheduleStatStreams();
	
	public:
		~StreamScheduler() override;
};
