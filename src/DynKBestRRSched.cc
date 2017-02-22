/**
 * @file DynKBestRRSched.cc
 * Implementation of the DynKBestRRSched class.
 */

#include "DynKBestRRSched.h"
#include "includes.h"
#include "MessageTypes.h"
#include "SINR_m.h"
#include "TransReqList_m.h"
#include "util.h"

#include <algorithm>
#include <cmath>
#include <forward_list>
#include <fstream>
#include <iterator>
#include <numeric>
#include <vector>

using namespace omnetpp;
using std::ofstream;
using std::set;
using std::unordered_map;
using std::vector;

Define_Module(DynKBestRRSched);

void DynKBestRRSched::initialize(){
	KBestRRStreamScheduler::initialize();
	packetLength = par("packetLength");
}

std::set<int>::iterator DynKBestRRSched::scheduleKBest(
		std::set<int>::iterator iter,std::vector<int>& blocks,
		MessageDirection dir,unsigned k){
	bool assigned = true;
	while(!blocks.empty()){
		if(iter==allOrigins.end()){
			// We're at the end of the set of senders, start at the beginning
			iter = allOrigins.begin();
			// If resource blocks were assigned in the previous run through the list,
			// set assigned to false and start another run. If there were no 
			// assignments in the previous run, break the loop. No sender needs 
			// resource blocks.
			if(assigned){
				assigned = false;
			}
			else{
				break;
			}
		}
		if(currOrigins.find(*iter)==currOrigins.end() || requests[*iter][dir].empty()){
			// Current origin does not have a request for the current direction 
			// and thus will be skipped
			++iter;
			continue;
		}
		if(*iter==-1 && dir==MessageDirection::up){
			// Current origin is the base station and we are assigning resource blocks
			// in the UP frequencies, which the BS does not need.
			++iter;
			continue;
		}
		SINR *estimate;
		int id = *iter;
		if(id==-1){
			estimate = estimateBS;
		}
		else{
			estimate = sinrEstimate[id];
		}
		auto blockComp = [&](int& first, int& second) -> bool {
			if(dir==MessageDirection::up){
				return estimate->getUp(first) > estimate->getUp(second);
			}
			else{
				return estimate->getDown(first) > estimate->getDown(second);
			}
		};
		// Sort the list of resource blocks in decreasing order, so that 
		// the first k elements are also the first k best resource blocks with the 
		// highest SINR.
		std::sort(blocks.begin(),blocks.end(),blockComp);
		auto bIter = blocks.begin();
		double sum = 0.0;
		while(bIter!= blocks.end() && sum<packetLength){
			if(dir==MessageDirection::up){
				sum += estimate->getRUp(*bIter);
			}
			else{
				sum += estimate->getRDown(*bIter);
			}
			++bIter;
		}
		// Assign the k best resource blocks to station id and remove them from 
		// the list of available resource blocks.
		std::move(blocks.begin(),bIter,std::back_inserter(originAssignments[id][dir]));
		blocks.erase(blocks.begin(),bIter);
		++iter;
		assigned = true;
	}
	return iter;
}

std::set<int>::iterator DynKBestRRSched::scheduleKBestStatic(
		std::set<int>::iterator iter,
		std::vector<int> blocks,MessageDirection dir,unsigned k,
		std::unordered_map<int,std::vector<int>>& schedule){
	while(!blocks.empty()){
		if(iter==allOrigins.end()){
			// We're at the end of the set of senders, start at the beginning
			iter = allOrigins.begin();
		}
		if(*iter==-1 && dir==MessageDirection::up){
			// Current origin is the base station and we are assigning resource blocks
			// in the UP frequencies, which the BS does not need.
			++iter;
			continue;
		}
		SINR *estimate;
		int id = *iter;
		if(id==-1){
			estimate = longtermEstimateBS;
		}
		else{
			estimate = longtermSinrEstimate[id];
		}
		auto blockComp = [&](int& first, int& second) -> bool {
			if(dir==MessageDirection::up){
				return estimate->getUp(first) > estimate->getUp(second);
			}
			else{
				return estimate->getDown(first) > estimate->getDown(second);
			}
		};
		// Sort the list of resource blocks in decreasing order, so that 
		// the first k elements are also the first k best resource blocks with the 
		// highest SINR.
		std::sort(blocks.begin(),blocks.end(),blockComp);
		auto bIter = blocks.begin();
		double sum = 0.0;
		while(bIter!= blocks.end() && sum<packetLength){
			if(dir==MessageDirection::up){
				sum += estimate->getRUp(*bIter);
			}
			else{
				sum += estimate->getRDown(*bIter);
			}
			++bIter;
		}
		// Assign the k best resource blocks to station id and remove them from 
		// the list of available resource blocks.
		std::move(blocks.begin(),bIter,std::back_inserter(schedule[id]));
		blocks.erase(blocks.begin(),bIter);
		++iter;
	}
	return iter;
}
