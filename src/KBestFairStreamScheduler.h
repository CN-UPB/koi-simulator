/**
 * @class KBestFairStreamScheduler KBestFairStreamScheduler.h
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

#include <fstream>
#include <set>
#include <unordered_map>
#include <vector>

class KBestFairStreamScheduler: public StreamScheduler{
	private:
		std::vector<int> allOrigins;
		std::set<int> currOrigins;
		std::unordered_map<int,std::unordered_map<int,std::vector<StreamTransReq*>>> requests;
		std::unordered_map<int,std::unordered_map<int,std::vector<int>>> originAssignments;
		std::unordered_map<int,std::unordered_map<int,double>> packCount;
		int upK;
		int downK;
		int bsId;
		int packetLength;
		int numSubcarriers;
		std::ofstream upSchedule;
		std::ofstream downSchedule;
		void scheduleKBest(std::vector<int>& blocks,MessageDirection dir,int k);
    virtual void scheduleStreams();
	
	protected:
    virtual void initialize();
		virtual void finish();
		virtual void handleMessage(cMessage *msg);
};
