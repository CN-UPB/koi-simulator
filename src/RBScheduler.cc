/**
 * \file RBScheduler.cc
 *
 * This file contains the default implementation of the RBScheduler interface.
 *
 * This resource block scheduler schedules packets depending on their creation
 * time. The oldest packet gets to be send first.
 */

#include "RBScheduler.h"
#include "MessageTypes.h"
#include "ScheduleInfo_m.h"
#include "TransReqList_m.h"
#include "util.h"
#include <algorithm>
#include <list>

using std::vector;
using std::list;
using std::unordered_map;

Define_Module(RBScheduler);

void RBScheduler::initialize(){
	this->rbNumber = par("rbNumber");
	this->numSubcarriers = par("numSubcarriers");
	ScheduleInfo *s = new ScheduleInfo();
	s->setSortfn(comparator);
	scheduleAt(simTime(),s);
}

StreamTransSched *RBScheduler::getSchedule(
		std::vector<StreamTransReq*>& reqs,
		int direction,
		const std::shared_ptr<std::unordered_map<int,SINR*>> estimates){
	if(!reqs.empty()){
		StreamTransReq *bestReq = nullptr;
		KoiData *bestPacket = nullptr;
		// The following code iterates over all transmission requests and finds 
		// the rquest with the best packet, according to the scheduler's "compare" 
		// member.
		StreamTransReq *currReq;
		for(auto iter=reqs.begin(); iter!=reqs.end(); iter++){
			currReq = *iter;
			// Iterate over all packets in the queue for the current stream
			// and find the which has the lowest value according to the 
			// scheduler's compare method.
			KoiData *currPacket;
			for(auto iterQueue = (currReq->getPackets())->begin();
					iterQueue!=(currReq->getPackets())->end();
					++iterQueue
				 ){
				currPacket = *iterQueue;
				if(!currPacket->getScheduled() 
						&& (bestPacket==nullptr || comparator(currPacket,bestPacket))){
					// The current packet is better than 
					// the previous best packet, so it 
					// becomes the new package to be 
					// scheduled.
					bestReq = currReq;
					bestPacket = currPacket;
					// Packets are sorted according to comparator, no further searching
					// needed as there are not going to be better packets deeper into 
					// the queue.
					break;
				}
			}
		}
		// Check whether we found a best request. That search might fail when 
		// all packets have already been scheduled.
		if(bestReq==nullptr){
			return nullptr;
		}
		// At this point, we know the request with the best packet.
		// The BS/MS where that request originated now gets as many 
		// packets scheduled as the expected transmission rate allows.
		// First, remove all requests from other BS/MS
		auto delIter(reqs.begin());
		while(delIter!=reqs.end()){
			if((*delIter)->getRequestOrigin()!=bestReq->getRequestOrigin()
					|| (*delIter)->getDest()!=bestReq->getDest()){
				delIter = reqs.erase(delIter);
			}
			else{
				++delIter;
			}
		}
		// Determine potential channel capacity.
		double channelCap;
		if(direction==MessageDirection::up){
			channelCap = estimates->at(bestReq->getRequestOrigin())->getRUp(rbNumber);
		}
		else{
			channelCap = estimates->at(bestReq->getRequestOrigin())->getRDown(rbNumber);
		}
		// Schedule the best packets according to the compare method
		// until the channel capacity is used up.
		while(channelCap>0){
			StreamTransReq *currReq = nullptr;
			bestPacket = nullptr;
			for(auto iter=reqs.begin(); iter!=reqs.end(); iter++){
				currReq = *iter;
				// Iterate over all packets in the queue for the current stream
				// and find the which has the lowest value according to the 
				// scheduler's compare method.
				KoiData *currPacket;
				for(auto iterQueue = (currReq->getPackets())->begin();
						iterQueue!=(currReq->getPackets())->end();
						++iterQueue
					 ){
					currPacket = *iterQueue;
					if(!currPacket->getScheduled() 
							&& (bestPacket==nullptr || comparator(currPacket,bestPacket))){
						// The current packet is better than 
						// the previous best packet, so it 
						// becomes the new package to be 
						// scheduled.
						bestPacket = currPacket;
						break;
					}
				}
			}
			if(bestPacket!=nullptr){
				bestPacket->setResourceBlock(rbNumber);
				int dir;
				if(bestPacket->getD2d()){
					if(direction==MessageDirection::down){
						dir = MessageDirection::d2dDown;
					}
					else{
						dir = MessageDirection::d2dUp;
					}
				}
				else{
					dir = direction;
				}
				bestPacket->setMessageDirection(dir);
				if(channelCap-bestPacket->getBitLength()>0){
					channelCap -= bestPacket->getBitLength();
					bestPacket->setScheduled(true);
				}
				else{
					bestPacket->setBitLength(bestPacket->getBitLength()-channelCap);
					KoiData *leftover = new KoiData(*bestPacket);
					leftover->setBitLength(channelCap);
					leftover->setScheduled(true);
					list<KoiData*>* squeue(currReq->getPackets());
					auto p = std::lower_bound(squeue->begin(),squeue->end(),leftover,
							comparator);
					squeue->insert(p,leftover);
					channelCap = 0;
				}
			}
			else{
				// No new best packet has been found, meaning that all packets 
				// of the choosen station have been scheduled. The remaining capacity
				// goes unused.
				break;
			}
		}
		StreamTransSched *schedule = new StreamTransSched();
		schedule->setSrc(bestReq->getRequestOrigin());
		schedule->setMessageDirection(direction);
		return schedule;
	}
	else{
		return nullptr;
	}
}

void RBScheduler::handleMessage(cMessage *msg){
	if(msg->getKind()==MessageType::transReqList){
		TransReqList *req = dynamic_cast<TransReqList*>(msg);
		StreamTransSched *sched = this->getSchedule(req->getRequests(),
				req->getMessageDirection(),
				req->getEstimates());
		if(sched!=nullptr){
			if(sched->getMessageDirection()==MessageDirection::d2d){
				switch(req->getMessageDirection()){
					case MessageDirection::down:
						sched->setMessageDirection(MessageDirection::d2dDown);
						break;
					case MessageDirection::up:
						sched->setMessageDirection(MessageDirection::d2dUp);
						break;
				}
			}
			send(sched,"scheduler$o");
		}
		delete req;
	}
	else if(msg->getKind()==MessageType::scheduleInfo){
		send(msg,"scheduler$o");
	}
}
