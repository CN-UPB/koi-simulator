/**
 * @class KBestStreamScheduler KBestStreamScheduler.h
 *
 * This scheduler assigns k RB to each sender via Round Robin
 *
 * Instead of assigning streams to resource blocks until all resource blocks 
 * are in use, this round robin variant assigns k resource blocks to each sender
 * until all resource blocks are assigned and continues with assignment in the 
 * next TTI where it stopped in the previous one.
 *
 * Thus, not all senders get to send in the current TTI, but those that do 
 * can utilize multiple resource blocks.
 */

#pragma once

#include "StreamScheduler.h"
#include "MessageTypes.h"
#include "StreamTransReq_m.h"

#include <set>
#include <unordered_map>
#include <vector>

class KBestStreamScheduler: public StreamScheduler{
	private:
		std::set<int>::iterator originUpIter;
		std::set<int>::iterator originDownIter;
		std::set<int> allOrigins;
		std::set<int> currOrigins;
		std::unordered_map<int,std::unordered_map<int,std::vector<StreamTransReq*>>> requests;
		std::unordered_map<int,std::unordered_map<int,std::vector<int>>> originAssignments;
		int upK;
		int downK;
		std::set<int>::iterator scheduleKBest(std::set<int>::iterator iter,
				std::vector<int>& blocks,MessageDirection dir,int k);
    virtual void scheduleStreams();
	
	protected:
    virtual void initialize();
		virtual void handleMessage(cMessage *msg);
};
