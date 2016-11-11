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

#include <functional>
#include <memory>
#include <unordered_map>
#include <vector>

class RBScheduler: public cSimpleModule{
	private:
		int rbNumber;
		int numSubcarriers;
		const std::function<bool(const KoiData*,const KoiData*)> comparator = 
			[](const KoiData *left,const KoiData *right) -> bool{
				return left->getCreationTime()<right->getCreationTime();
			};

		virtual StreamTransSched *getSchedule(
				std::vector<StreamTransReq*>& reqs,
				int direction,
				const std::shared_ptr<std::unordered_map<int,SINR*>> estimates);

	protected:
		virtual void initialize();
		virtual void handleMessage(cMessage *msg);
	
	public:
		~RBScheduler() = default;
};
