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
#include <vector>

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
	vector<int> partnerIds = cStringTokenizer(partners,",").asIntVector();
	// schedule initial events randomly in [current time,initOffset]
	simtime_t currTime = simTime();
	unsigned long sId;
	StreamDef stream{0};
	for(int i:partnerIds){
		sId = simulation.getUniqueNumber();
		stream = StreamDef{sId,bsId,i,period,false};
		streams[sId] = stream;
		if(this->periodicTraffic){
			Traffic *msg = new Traffic();
			msg->setPartner(stream.destMsId);
			msg->setTrafficType(TrafficType::periodic);
			msg->setStreamId(stream.streamId);
			scheduleAt(currTime+uniform(currTime,currTime+this->initOffset),
					msg);
		}
		// Send stream notification message, which will inform the 
		// mobile station's MAC, the base station's MAC and the 
		// stream scheduler about the communication stream between 
		// this MS and it's comm partner.
		StreamInfo *info = new StreamInfo();
		info->setSrc(msId);
		info->setDest(stream.destMsId);
		info->setInterarrival(stream.period);
		info->setStreamId(stream.streamId);
		send(info,"toMac");
	}
}

void TrafficGen::handleMessage(cMessage *msg){
	simtime_t currTime = simTime();
	if(msg->isSelfMessage()) {
		switch(msg->getKind()){
			case MessageType::traffic:
				Traffic *genMsg = dynamic_cast<Traffic*>(msg);
				KoiData *pack = new KoiData();
				StreamDef stream = streams[genMsg->getStreamId()];
				pack->setBsId(this->bsId);
				pack->setDeadline(currTime+this->deadline);
				pack->setSrc(this->msId);
				pack->setDest(stream.destMsId);
				pack->setTrafficType(TrafficType::periodic);
				pack->setBitLength(this->packetLength);
				pack->setInterarrival(stream.period);
				pack->setStreamId(stream.streamId);
				send(pack,"toMac");
				scheduleAt(currTime+this->period,msg);
				std::cout << "Stream " << stream.streamId << ": " 
					<< "MS " << this->msId 
					<< " Generated Msg " << pack->getId() 
					<< " to " << pack->getDest() 
					<< std::endl;
				break;
			
		}
	}
	else if(msg->arrivedOn("fromMac")){
		KoiData *pack = dynamic_cast<KoiData*>(msg);
		std::cout << "Stream " << pack->getStreamId() << ": " 
			<< "MS " << this->msId 
			<< " got message " << pack->getId() 
			<< " from " << pack->getSrc() 
			<< " to " << pack->getDest() 
			<< std::endl;
		delete pack;
	}
}
