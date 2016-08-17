/**
 * \class EDFRBScheduler EDFRBScheduler.h
 *
 * The EDF Packet Scheduler for a specific Resource Block.
 *
 * This scheduler class is responsible for scheduling packet transmissions on 
 * a single Resource Block. In each TTI, the packet with the earliest deadline 
 * among all packets awaiting transmission from all streams assigned to this 
 * resource block will be scheduled for transmission.
 */

#pragma once

#include "includes.h"
#include "RBScheduler.h"
#include "StreamTransSched_m.h"
#include "StreamTransReq_m.h"
#include "KoiData_m.h"

#include <vector>

class EDFRBSched: public RBScheduler{
	public:
		virtual bool comparator(const KoiData *left, 
				const KoiData *right) const;
		~EDFRBSched() = default;
};
