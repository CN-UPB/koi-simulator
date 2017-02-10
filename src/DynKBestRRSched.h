/**
 * @class DynKBestRRSched DynKBestRRSched.h
 *
 * Assign as many RB to each sender via Round Robin as needed for one packet
 *
 * Instead of assigning streams to resource blocks until all resource blocks 
 * are in use, this round robin variant assigns k resource blocks to each sender
 * until all resource blocks are assigned and continues with assignment in the 
 * next TTI where it stopped in the previous one.
 *
 * In this variant, k is dynamic. This means that each sender gets as many RB
 * as he needs to send at least one packet.
 *
 * Thus, not all senders get to send in the current TTI, but those that do 
 * can utilize multiple resource blocks.
 */

#pragma once

#include "KBestRRStreamScheduler.h"
#include "MessageTypes.h"
#include "StreamTransReq_m.h"

#include <fstream>
#include <set>
#include <unordered_map>
#include <vector>

class DynKBestRRSched: public KBestRRStreamScheduler{
	private:
		int packetLength;

	protected:
		std::set<int>::iterator scheduleKBest(std::set<int>::iterator iter,
				std::vector<int>& blocks,MessageDirection dir,unsigned k) override;
		std::set<int>::iterator scheduleKBestStatic(
				std::set<int>::iterator iter,
				std::vector<int> blocks,MessageDirection dir,unsigned k,
				std::unordered_map<int,std::vector<int>>& schedule) override;
    void initialize() override;
};
