/**
 * \class RBScheduler RBScheduler.h
 *
 * The Packet Scheduler for a specific Resource Block.
 *
 * This scheduler class is responsible for scheduling packet transmissions on 
 * a single Resource Block. It does not have a direct connection to Mobile and 
 * Base stations. Instead, it gets all requests it should schedule from the 
 * StreamScheduler, which is responsible for assigning streams to Resource 
 * Blocks. 
 *
 * To schedule packets, the RBScheduler receives all packets for all streams 
 * it is responsible for, and then decides which of these packages gets to be 
 * send in the current TTI.
 */

#pragma once

#include "includes.h"
#include "MessageTypes.h"
#include "SINR_m.h"
#include "StreamTransSched_m.h"
#include "StreamTransReq_m.h"
#include "KoiData_m.h"

#include <unordered_map>
#include <vector>

class RBScheduler: public cSimpleModule{
	private:
		int rbNumber;
		virtual StreamTransSched *getSchedule(
				std::vector<StreamTransReq*>& reqs,
				int direction,
				const std::unordered_map<int,SINR*>* estimates);

	protected:
		virtual void initialize();
		virtual void handleMessage(cMessage *msg);
	
	public:
		virtual bool comparator(const KoiData *left,const KoiData *right) const;
		~RBScheduler() = default;
};
