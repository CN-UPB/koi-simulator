/**
 * \file TrafficGen.cc
 * Implementation of the KOI traffic generator
 */

#include "TrafficGen.h"
#include "KoiData_m.h"
#include "Traffic_m.h"
#include "StreamInfo_m.h"
#include "MessageTypes.h"

#include <iostream>

using std::vector;

Define_Module(TrafficGen);

void TrafficGen::initialize(){
	this->bsId = par("bsId");
	this->deadline = par("deadline");
	this->initOffset = par("initOffset");
	this->msId = par("msId");
	this->packetLength = par("packetLength");
	this->period = par("period");
	this->periodicTraffic = par("periodicTraffic");

	// convert communication partners string to vector of ints
	const char *partners = par("commPartners").stringValue();
	this->commPartners = cStringTokenizer(partners).asIntVector();
	// schedule initial events randomly in [current time,initOffset]
	simtime_t currTime = simTime();
	for(int i:this->commPartners){
		if(this->periodicTraffic){
			Traffic *msg = new Traffic();
			msg->setPartner(i);
			msg->setTrafficType(TrafficType::periodic);
			scheduleAt(currTime+uniform(currTime,currTime+this->initOffset),
					msg);
		}
		// Send stream notification message, which will inform the 
		// mobile station's MAC, the base station's MAC and the 
		// stream scheduler about the communication stream between 
		// this MS and it's comm partner.
		StreamInfo *info = new StreamInfo();
		info->setSrc(msId);
		info->setDest(i);
		info->setInterarrival(period);
		send(info,"toMac");
		std::cout << "Sending successfull" << std::endl;
	}
}

void TrafficGen::handleMessage(cMessage *msg){
	simtime_t currTime = simTime();
	if(msg->isSelfMessage()) {
		switch(msg->getKind()){
			case MessageType::traffic:
				Traffic *genMsg = dynamic_cast<Traffic*>(msg);
				KoiData *pack = new KoiData();
				pack->setBsId(this->bsId);
				pack->setDeadline(currTime+this->deadline);
				pack->setSrc(this->msId);
				pack->setDest(genMsg->getPartner());
				pack->setTrafficType(TrafficType::periodic);
				pack->setBitLength(this->packetLength);
				pack->setInterarrival(this->period);
				send(pack,"toMac");
				scheduleAt(currTime+this->period,msg);
				std::cout << "MS " << this->msId << " send msg " << pack->getId() << std::endl;
				break;
			
		}
	}
	else if(msg->arrivedOn("fromMac")){
		KoiData *pack = dynamic_cast<KoiData*>(msg);
		std::cout << "MS " << this->msId 
			<< " got message " << pack->getId() 
			<< " from " << pack->getSrc() 
			<< std::endl;
		delete pack;
	}
}
