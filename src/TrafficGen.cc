/**
 * \file TrafficGen.cc
 * Implementation of the KOI traffic generator
 */

#include "TrafficGen.h"
#include "KoiData_m.h"
#include "Traffic_m.h"
#include "StreamInfo_m.h"
#include "MessageTypes.h"
#include "omnetpp.h"

#include <iostream>
#include <vector>

using std::vector;
using std::string;

Define_Module(TrafficGen);

void TrafficGen::initialize(){
	this->bsId = par("bsId");
	this->initOffset = par("initOffset");
	this->msId = par("msId");
	this->packetLength = par("packetLength");
	this->periodicTraffic = par("periodicTraffic");
        this->genPackets = registerSignal("genPackets");
        this->recPackets = registerSignal("recPackets");
        this->missedDL = registerSignal("missedDL");

	string xmlPath = par("commTable");
	vector<StreamDef> parsedStreams(parseCommTable(xmlPath,bsId,msId));
	simtime_t currTime = simTime();
	for(StreamDef& stream:parsedStreams){
		this->streams[stream.streamId] = stream;
		if(this->periodicTraffic){
			Traffic *msg = new Traffic();
                        msg->setPartner(stream.destMsId);
			msg->setTrafficType(TrafficType::periodic);
                        msg->setStreamId(stream.streamId);
                        // schedule initial events randomly in [current time,initOffset]
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
		info->setD2d(stream.d2d);
                info->setDeadline(stream.deadline);
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
				StreamDef& stream = streams[genMsg->getStreamId()];
				pack->setBsId(this->bsId);
				pack->setDeadline(currTime+stream.deadline);
				pack->setSrc(this->msId);
				pack->setDest(stream.destMsId);
				pack->setTrafficType(TrafficType::periodic);
				pack->setBitLength(this->packetLength);
				pack->setInterarrival(stream.period);
				pack->setStreamId(stream.streamId);
				pack->setD2d(stream.d2d);
				send(pack,"toMac");
				scheduleAt(currTime+stream.period,msg);
				std::cout << "Stream " << stream.streamId << ": " 
					<< "MS " << this->msId 
					<< " Generated Msg " << pack->getId() 
					<< " to " << pack->getDest() 
					<< std::endl;
                                // Emit signal for statistics gathering
                                emit(genPackets,true);
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
                emit(recPackets,true);
                if(pack->getDeadline()>simTime()){
                  emit(missedDL,pack->getDeadline()-simTime());
                }
		delete pack;
	}
}

vector<TrafficGen::StreamDef> TrafficGen::parseCommTable(const string& commTable,
		int bsId,
		int msId){
	vector<TrafficGen::StreamDef> parsedStreams;
	string xpath("/root/cell[@id='"
			+std::to_string(bsId)
			+"']/ms[@id='"+std::to_string(msId)+"']");
	cXMLElement *msNode = ev.getXMLDocument(commTable.c_str(),
			xpath.c_str());
	unsigned long sId;
	for(cXMLElement *curr=msNode->getFirstChild();curr!=nullptr;
			curr=curr->getNextSibling()){
		sId = simulation.getUniqueNumber();
		parsedStreams.emplace_back(
			sId,
			std::stoi(curr->getAttribute("cell")),
			std::stoi(curr->getAttribute("destid")),
			std::stod(curr->getAttribute("period")),
			std::stoi(curr->getAttribute("d2d")),
			std::stod(curr->getAttribute("deadline")));
	}
	return parsedStreams;
}
