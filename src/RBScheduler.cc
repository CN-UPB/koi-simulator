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
#include <list>

using std::vector;
using std::list;

Define_Module(RBScheduler);

void RBScheduler::initialize(){
	this->rbNumber = par("rbNumber");
}

StreamTransSched *RBScheduler::getSchedule(
		const std::vector<StreamTransReq*>& reqs){
	if(!reqs.empty()){
                list<KoiData*>::const_iterator bestPos;
                bool init = false;
		int bestDir;
		// The following code iterates over all transmission requests and finds 
		// the next packet to be send, according to the scheduler's "compare" 
		// method.
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
				if(!init || comparator(currPacket,*bestPos)){
					// The current packet is better than 
					// the previous best packet, so it 
					// becomes the new package to be 
					// scheduled.
                                        bestPos = iterQueue;
					bestDir = currReq->getMessageDirection();
                                        // We need to indicate that the
                                        // iterator has been initialized, 
                                        // because a default constructed iterator
                                        // cannot be checked
                                        init = true;
				}
			}
			delete currReq;
		}
		StreamTransSched *schedule = new StreamTransSched();
		schedule->setSrc((*bestPos)->getSrc());
		schedule->setDest((*bestPos)->getDest());
		schedule->setStreamId((*bestPos)->getStreamId());
		schedule->setRb(this->rbNumber);
		schedule->setPacketPos(bestPos);
		schedule->setMessageDirection(bestDir);
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
}

bool RBScheduler::comparator(const KoiData *left, const KoiData *right) const{
	return left->getCreationTime()<right->getCreationTime();
}
