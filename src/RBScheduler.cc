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
#include "TransReqList_m.h"

using std::vector;

Define_Module(RBScheduler);

void RBScheduler::initialize(){
	this->rbNumber = par("rbNumber");
}

StreamTransSched *RBScheduler::getSchedule(
		const std::vector<StreamTransReq*>& reqs){
	if(!reqs.empty()){
		KoiData *bestPacket = nullptr;
		int bestPacketIndex;
		int bestSrc;
		int bestDest;
		bool bestBs;
		// The following code iterates over all transmission requests and finds 
		// the next packet to be send, according to the scheduler's "compare" 
		// method.
		StreamTransReq *currReq;
		for(auto iter=reqs.begin(); iter!=reqs.end(); iter++){
			currReq = *iter;
			// Iterate over all packets in the queue for the current stream
			// and find the which has the lowest value according to the 
			// scheduler's compare method.
			int index = 0;
			KoiData *currPacket;
			for(cQueue::Iterator iterQueue(
						*(currReq->getPackets()));
					!iterQueue.end();
					){
				currPacket = dynamic_cast<KoiData*>(iterQueue++);
				if(bestPacket==nullptr 
						|| comparator(currPacket,bestPacket)){
					// The current packet is better than 
					// the previous best packet, so it 
					// becomes the new package to be 
					// scheduled.
					bestPacket = currPacket;
					bestPacketIndex = index;
					bestSrc = bestPacket->getSrc();
					bestDest = bestPacket->getDest();
					bestBs = currReq->getBs();					
				}
				index++;
			}
			delete currReq;
		}
		StreamTransSched *schedule = new StreamTransSched();
		schedule->setSrc(bestSrc);
		schedule->setDest(bestDest);
		schedule->setRb(this->rbNumber);
		schedule->setPacketIndex(bestPacketIndex);
		schedule->setBs(bestBs);
		return schedule;
	}
	else{
		return nullptr;
	}
}

void RBScheduler::handleMessage(cMessage *msg){
	if(msg->getKind()==MessageType::transReqList){
		TransReqList *req = dynamic_cast<TransReqList*>(msg);
		StreamTransSched *sched = this->getSchedule(req->getRequests());
		if(sched!=nullptr){
			send(sched,"toScheduler");
		}
	}
}

bool RBScheduler::comparator(const KoiData *left, const KoiData *right) const{
	return left->getCreationTime()<right->getCreationTime();
}
