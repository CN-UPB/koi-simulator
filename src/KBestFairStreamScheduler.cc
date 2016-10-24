/**
 * @file KBestFairStreamScheduler.cc
 * Implementation of the KBestFairStreamScheduler class.
 */

#include "KBestFairStreamScheduler.h"
#include "includes.h"
#include "MessageTypes.h"
#include "SINR_m.h"
#include "TransReqList_m.h"
#include "util.h"

#include <algorithm>
#include <cmath>
#include <forward_list>
#include <fstream>
#include <numeric>
#include <vector>

using std::forward_list;
using std::ofstream;
using std::set;
using std::unordered_map;
using std::vector;

Define_Module(KBestFairStreamScheduler);

void KBestFairStreamScheduler::initialize(){
	StreamScheduler::initialize();
	this->upK = par("upK");
	this->downK = par("downK");
	this->bsId = par("bsId");
	this->packetLength = par("packetLength");
	this->numSubcarriers = par("numSubcarriers");
	// Fill the origins set with all possible origins for transmission requests,
	// a.k.a all local mobile stations plus the local base station.
	// First add id -1, representing the local base station
	allOrigins.push_back(-1);
	// Then insert all mobile station ids
	for(int i=0; i<numberOfMs; i++){
		allOrigins.push_back(i);
	}
	// Set initial number of packets send last TTI to 0
	for(int& id:allOrigins){
		packCount[id][MessageDirection::up] = 0.0;
		packCount[id][MessageDirection::down] = 0.0;
	}
	// Prepare result files for schedules
	std::string fname("schedule-up-cell-"+std::to_string(bsId));
	upSchedule = std::move(getResultFile(fname));
	upSchedule << "TTI\t" << "Cell\t" 
		<< "MS\t" << "RB\t" << "SINR" << std::endl;
	fname = "schedule-down-cell-"+std::to_string(bsId);
	downSchedule = std::move(getResultFile(fname));
	downSchedule << "TTI\t" << "Cell\t"
		<< "MS\t" << "RB\t" << "SINR" << std::endl;
}

void KBestFairStreamScheduler::finish(){
	upSchedule.close();
	downSchedule.close();
}

void KBestFairStreamScheduler::scheduleKBest(
		std::vector<int>& blocks, MessageDirection dir,int k){
	bool assigned = true;
	auto iter = allOrigins.begin();
	while(blocks.size()>=k){
		if(iter==allOrigins.end()){
			// We're at the end of the set of senders, start at the beginning
			iter = allOrigins.begin();
			// If resource blocks were assigned in the previous run through the list,
			// set assigned to false and start another run. If there were no 
			// assignments in the previous run, break the loop. No sender needs 
			// resource blocks.
			if(assigned==true){
				assigned = false;
			}
			else{
				break;
			}
		}
		if(currOrigins.find(*iter)==currOrigins.end()){
			// Current origin does not have a request and thus will be skipped
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
		// Assign the k best resource blocks to station id and remove them from 
		// the list of available resource blocks.
		originAssignments[id][dir].resize(k);
		std::move(blocks.begin(),blocks.begin()+k,originAssignments[id][dir].begin());
		blocks.erase(blocks.begin(),blocks.begin()+k);
		if(simTime()>initOffset){
			auto val = (simTime()-initOffset)/tti;
			int tti = std::floor(val);
			if(dir==MessageDirection::up){
				for(auto assRB:originAssignments[id][dir]){
					upSchedule 
						<< tti << "\t" 
						<< bsId << "\t" 
						<< id << "\t"
						<< assRB << "\t" 
						<< estimate->getUp(assRB)
						<< std::endl;
				}
			}
			else if(dir==MessageDirection::down){
				downSchedule << tti << bsId << id;
				for(auto assRB:originAssignments[id][dir]){
					downSchedule 
						<< tti << "\t" 
						<< bsId << "\t" 
						<< id << "\t"
						<< assRB << "\t" 
						<< estimate->getDown(assRB)
						<< std::endl;
				}
			}
		}
		// Set send packets count
		double cap = 0.0;
		vector<double> vals(1);
		for(auto assRB:originAssignments[id][dir]){
			if(dir==MessageDirection::down){
				vals[0] = estimate->getDown(assRB);
			}
			else{
				vals[0] = estimate->getUp(assRB);
			}
			cap += getChannelCapacity(vals,numSubcarriers);
		}
		packCount[id][dir] = cap/packetLength;
		++iter;
		assigned = true;
	}
}

void KBestFairStreamScheduler::scheduleStreams(){
	// Clear out the current assignments
	originAssignments.clear();
	// Build lists of resource block for UP/DOWN bands
	vector<int> upBlocks(upRB);
	vector<int> downBlocks(downRB);
	// Fill lists with resource block numbers
	std::iota(upBlocks.begin(),upBlocks.end(),0);
	std::iota(downBlocks.begin(),downBlocks.end(),0);
	// Sort allOrigins ascending by number of packets send last tti
	MessageDirection dir = MessageDirection::up;
	auto compSendPackUp = [&](int& first, int& second) -> bool {
		return packCount[first][dir] < packCount[second][dir];
	};
	std::sort(allOrigins.begin(),allOrigins.end(),compSendPackUp);
	for(int& id:allOrigins){
		packCount[id][MessageDirection::up] = 0.0;
	}
	// Schedule UP blocks	
	scheduleKBest(upBlocks,MessageDirection::up,upK);
	dir = MessageDirection::down;
	auto compSendPackDown = [&](int& first, int& second) -> bool {
		return packCount[first][dir] < packCount[second][dir];
	};
	std::sort(allOrigins.begin(),allOrigins.end(),compSendPackDown);
	for(int& id:allOrigins){
		packCount[id][MessageDirection::down] = 0.0;
	}
	// Schedule DOWN blocks
	scheduleKBest(downBlocks,MessageDirection::down,downK);
	currOrigins.clear();
}

void KBestFairStreamScheduler::handleMessage(cMessage *msg){
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
					std::shared_ptr<unordered_map<int,SINR*>> estimates = std::make_shared<unordered_map<int,SINR*>>();
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
								send(lst->dup(),"upRB$o",rb);
								break;
							case MessageDirection::down:
								send(lst->dup(),"downRB$o",rb);
								break;
						}
					}
					delete lst;
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
			forward_list<StreamTransReq*> d2dReqs;
			for(auto iterOrig=this->requests.begin();
				 iterOrig!=requests.end();
				 ++iterOrig){
				// Iterate over both UP/DOWN transmission bands
				for(auto iterDir=iterOrig->second.begin();
					 iterDir!=iterOrig->second.end();
					 ++iterDir){
					for(StreamTransReq* req:iterDir->second){
						if(iterDir->first==MessageDirection::up 
								&& req->getMessageDirection()==MessageDirection::d2d){
							// Special handling for d2d requests, because they are added to
							// the request lists for both transmission bands. But delete 
							// must only be called on a pointer once.
							// Thus, D2D requests only get deleted when found in the UP
							// list of requests.
							d2dReqs.push_front(req);
						}
						else if(req->getMessageDirection()!=MessageDirection::d2d){
							// Delete any non-D2D request
							delete req;
						}
					}
					iterDir->second.clear();
				}
				iterOrig->second.clear();
			}
			requests.clear();
			// Now delete all the D2D requests
			for(auto ptr:d2dReqs){
				delete ptr;
			}
		} break;
		default:
			// Message type is handled by StreamScheduler class
			StreamScheduler::handleMessage(msg);
	}
}
