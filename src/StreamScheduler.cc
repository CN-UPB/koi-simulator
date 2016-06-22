/**
 * \file StreamScheduler.cc
 *
 * This file contains the default implementation of the StreamScheduler interface.
 *
 * This implementation of the StreamScheduler implements simple round robin 
 * assignment of streams to resource blocks.
 */

#include "StreamScheduler.h"
#include "TransReqList_m.h"
#include "MessageTypes.h"
#include <iostream>

using std::vector;
using std::unordered_map;

Define_Module(StreamScheduler);

void StreamScheduler::initialize(){
	this->initOffset = par("initOffset");
	this->downRB = par("downResourceBlocks");
	this->upRB = par("upResourceBlocks");
	this->numberOfMs = par("numberOfMobileStations");
	this->infos = vector<StreamInfo*>();
	this->streamSchedPeriod = par("streamSchedPeriod");
	this->tti = par("tti");
	// Produce the first schedule right at the init offset
	// We will need to make certain that all mobile stations have reported 
	// their streams at the time the first schedule is computed.
	scheduleAt(simTime()+initOffset,
			new cMessage("",MessageType::scheduleStreams));
	scheduleAt(simTime()+initOffset,
			new cMessage("",MessageType::scheduleRBs));
}

void StreamScheduler::scheduleStreams(){
	if(!this->infos.empty()){
		// First, clear the current assignment
		this->rbAssignments.clear();
		// This simple algorithm assigns streams to resource blocks 
		// by a round robin.
		int rbUp = 0;
		int rbDown = 0;
		for(auto info:this->infos){
			if(info->getD2d()){
				// If the stream is a D2D link, the scheduler 
				// does not only need to decide which ressource 
				// block the stream should use, but also 
				// which frequency band half, the one assigned 
				// to down or up traffic.
				//
				// We don't need assignments for up/down here 
				// as D2D streams only have one direction, 
				// directly from the sending MS to the receiving
				// MS.
				if(rbUp/(double)upRB<=rbDown/(double)downRB){
					this->rbAssignments[info->getStreamId()][MessageDirection::d2d] 
						= std::make_pair<MessageDirection,int>(MessageDirection::up,rbUp%upRB);
					rbUp++;
				}
				else{
					this->rbAssignments[info->getStreamId()][MessageDirection::d2d] 
						= std::make_pair<MessageDirection,int>(MessageDirection::down,rbDown%downRB);
					rbDown++;
				}
			}
			else{
				this->rbAssignments[info->getStreamId()][MessageDirection::up] 
					= std::make_pair<MessageDirection,int>(MessageDirection::up,rbUp%upRB);
				this->rbAssignments[info->getStreamId()][MessageDirection::down] 
					= std::make_pair<MessageDirection,int>(MessageDirection::down,rbDown%downRB);
				rbUp++;
				rbDown++;
			}
			delete info;
		}
		// Remove all stream infos
		this->infos.clear();
	}
}

void StreamScheduler::handleMessage(cMessage *msg){
	switch(msg->getKind()){
		case MessageType::streamInfo:{
				StreamInfo *tmp = dynamic_cast<StreamInfo*>(msg);
				this->infos.push_back(tmp);
			} break;
		case MessageType::streamTransReq:{
				StreamTransReq *req = dynamic_cast<StreamTransReq*>(msg);
				ResAssign& assignment = rbAssignments[req->getStreamId()][req->getMessageDirection()];
				this->requests[assignment.first][assignment.second].push_back(req);
			} break;
		case MessageType::scheduleStreams:{
				this->scheduleStreams();
				scheduleAt(simTime()+this->streamSchedPeriod,msg);
			} break;
		case MessageType::scheduleRBs:
			// Iterate over all message directions (up/down)
			for(auto iterDir=this->requests.begin();
					iterDir!=requests.end();
					++iterDir){
				// Iterate over all Resource blocks in the current 
				// transmission direction
				for(auto iterRb=iterDir->second.begin();
						iterRb!=iterDir->second.end();
						++iterRb){
					TransReqList *lst = new TransReqList();
					lst->setRequests(iterRb->second);
					// Clear all transmission requests for this 
					// stream from the stream schedulers queue.
					iterRb->second.clear();
					switch(iterDir->first){
						case MessageDirection::up:
							send(lst,"upRB$o",iterRb->first);
							break;
						case MessageDirection::down:
							send(lst,"downRB$o",iterRb->first);
							break;
					}
				}
			}
			scheduleAt(simTime()+tti,msg);
			break;
		case MessageType::streamSched:{
				StreamTransSched *sched = dynamic_cast<StreamTransSched*>(msg);
				switch(sched->getMessageDirection()){
					case MessageDirection::up:
						send(sched,"toMs",sched->getSrc());
						break;
					case MessageDirection::down:
						send(sched,"toBs");
						break;
				}
			} break;
		default:
			std::cerr << "Received invalid Message in "
				  << "StreamScheduler handleMessage"
				  << std::endl;
	}
}
