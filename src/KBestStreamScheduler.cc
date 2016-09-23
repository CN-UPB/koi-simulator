/**
 * @file KBestStreamScheduler.cc
 * Implementation of the KBestStreamScheduler class.
 */

#include "KBestStreamScheduler.h"
#include "includes.h"
#include "MessageTypes.h"
#include "SINR_m.h"
#include "TransReqList_m.h"

#include <algorithm>
#include <numeric>
#include <vector>

using std::set;
using std::unordered_map;
using std::vector;

Define_Module(KBestStreamScheduler);

void KBestStreamScheduler::initialize(){
	StreamScheduler::initialize();
	this->upK = par("upK");
	this->downK = par("downK");
	// Fill the origins set with all possible origins for transmission requests,
	// a.k.a all local mobile stations plus the local base station.
	// First add id -1, representing the local base station
	allOrigins.insert(-1);
	// Then insert all mobile station ids
	for(int i=0; i<numberOfMs; i++){
		allOrigins.insert(i);
	}
	originUpIter = allOrigins.begin();
	originDownIter = allOrigins.begin();
}

std::set<int>::iterator KBestStreamScheduler::scheduleKBest(
		std::set<int>::iterator iter,std::vector<int>& blocks,
		MessageDirection dir,int k){
	while(blocks.size()>=k){
		if(iter==allOrigins.end()){
			// We're at the end of the set of senders, start at the beginning
			iter = allOrigins.begin();
		}
		if(currOrigins.find(*iter)==currOrigins.end()){
			// Current origin does not have a request and thus will be skipped
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
		// Assign the k best resource blocks to station id and remove them from 
		// the list of available resource blocks.
		originAssignments[id][dir].resize(k);
		std::move(blocks.begin(),blocks.begin()+k,originAssignments[id][dir].begin());
		blocks.erase(blocks.begin(),blocks.begin()+k);
		++iter;
	}
	return iter;
}

void KBestStreamScheduler::scheduleStreams(){
	// Clear out the current assignments
	originAssignments.clear();
	// Build lists of resource block for UP/DOWN bands
	vector<int> upBlocks(upRB);
	vector<int> downBlocks(downRB);
	// Fill lists with resource block numbers
	std::iota(upBlocks.begin(),upBlocks.end(),0);
	std::iota(downBlocks.begin(),downBlocks.end(),0);
	// Schedule UP blocks
	originUpIter = scheduleKBest(originUpIter,upBlocks,MessageDirection::up,upK);
	// Schedule DOWN blocks
	originDownIter = scheduleKBest(originDownIter,downBlocks,
			MessageDirection::down,downK);
	currOrigins.clear();
}

void KBestStreamScheduler::handleMessage(cMessage *msg){
	switch(msg->getKind()){
    case MessageType::scheduleStreams:{
			if(currOrigins.size()>0){
				// Only compute a schedule if there actually are any requests which 
				// would make use of the schedule.
				this->scheduleStreams();
			}
			scheduleAt(simTime()+this->streamSchedPeriod,msg);
		} break;
		case MessageType::streamTransReq:{
			StreamTransReq *req = dynamic_cast<StreamTransReq*>(msg);
			currOrigins.insert(req->getRequestOrigin());
			if(req->getMessageDirection()==MessageDirection::d2d){
				requests[req->getRequestOrigin()][MessageDirection::down].push_back(req);
				requests[req->getRequestOrigin()][MessageDirection::up].push_back(req);
			}
			else{
				requests[req->getRequestOrigin()][req->getMessageDirection()].push_back(req);
			}
		}	break;
		case MessageType::scheduleRBs:
			scheduledStations.clear();
			// Iterate over all origins (sending stations)
			for(auto iterOrig=this->requests.begin();
				 iterOrig!=requests.end();
				 ++iterOrig){
				// Iterate over both UP/DOWN transmission bands
				for(auto iterDir=iterOrig->second.begin();
					 iterDir!=iterOrig->second.end();
					 ++iterDir){
					TransReqList *lst = new TransReqList();
					lst->setRequests(iterDir->second);
					// Gather SINR estimates for all 
					// MS with streams in this request
					unordered_map<int,SINR*>* estimates = new unordered_map<int,SINR*>();
					if(iterOrig->first==-1){
						// Request from BS
						(*estimates)[-1] = estimateBS;
					}
					else{
						// Request from MS
						(*estimates)[iterOrig->first] = sinrEstimate[iterOrig->first];
					}
					lst->setEstimates(estimates);
					lst->setMessageDirection(iterDir->first);
					for(int& rb:originAssignments[iterOrig->first][iterDir->first]){
						switch(iterDir->first){
							case MessageDirection::up:
								send(lst,"upRB$o",rb);
								break;
							case MessageDirection::down:
								send(lst,"downRB$o",rb);
								break;
						}
					}
				}
			}
			scheduleAt(simTime()+tti,msg);
			scheduleAt(simTime()+epsilon,new cMessage("",MessageType::sendSchedules));
			break;
		case MessageType::sendSchedules:{
			StreamTransSched* sched(nullptr);
			for(const int& id:scheduledStations){
				sched = new StreamTransSched();
				sched->setSrc(id);
				if(id==-1){
					// -1 is the local base station
					send(sched,"toBs");
				}
				else{
					// Send to the indicated mobile station
					send(sched,"toMs",id);
				}
			}
			// At this point, all requests for the current TTI have been handled, 
			// and the messages can be deleted.
			for(auto iterOrig=this->requests.begin();
				 iterOrig!=requests.end();
				 ++iterOrig){
				// Iterate over both UP/DOWN transmission bands
				for(auto iterDir=iterOrig->second.begin();
					 iterDir!=iterOrig->second.end();
					 ++iterDir){
					for(StreamTransReq* req:iterDir->second){
						delete req;
					}
					iterDir->second.clear();
				}
				iterOrig->second.clear();
			}
			requests.clear();
		} break;
		default:
			// Message type is handled by StreamScheduler class
			StreamScheduler::handleMessage(msg);
	}
}
