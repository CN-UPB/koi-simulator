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
	this->infos = vector<StreamInfo*>(this->numberOfMs);
	this->initOffset = par("initOffset");
	this->resourceBlocks = par("resourceBlocks");
	this->numberOfMs = par("numberOfMobileStations");
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
		int rb = 0;
		for(auto info:this->infos){
			this->rbAssignments[info->getSrc()][info->getDest()] = rb%resourceBlocks;
			rb++;
			std::cout << "Did RB Assignment: "
				<< info->getSrc()
				<< "->"
				<< info->getDest()
				<< " to RB "
				<< this->rbAssignments[info->getSrc()][info->getDest()]
				<< std::endl;
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
				int assignedRB = rbAssignments[req->getSrc()][req->getDest()];
				this->requests[assignedRB].push_back(req);
			} break;
		case MessageType::scheduleStreams:{
				this->scheduleStreams();
				scheduleAt(simTime()+this->streamSchedPeriod,msg);
			} break;
		case MessageType::scheduleRBs:
			for(auto iter=this->requests.begin();
					iter!=requests.end();
					++iter){
				TransReqList *lst = new TransReqList();
				lst->setRequests(iter->second);
				// Clear all transmission requests for this 
				// stream from the stream schedulers queue.
				iter->second.clear();
				send(lst,"toRB",iter->first);
			}
			scheduleAt(simTime()+tti,msg);
			break;
		case MessageType::streamSched:{
				StreamTransSched *sched = dynamic_cast<StreamTransSched*>(msg);
				if(sched->getBs()){
					send(sched,"toBs");
				}
				else{
					send(sched,"toMs",sched->getSrc());
				}
			} break;
		default:
			std::cerr << "Received invalid Message in "
				  << "StreamScheduler handleMessage"
				  << std::endl;
	}
}
