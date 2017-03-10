/**
 * \file TrafficGen.cc
 * Implementation of the KOI traffic generator
 */

#include "TrafficGen.h"
#include "KoiData_m.h"
#include "Traffic_m.h"
#include "ResultFileExchange_m.h"
#include "StreamInfo_m.h"
#include "MessageTypes.h"
#include "omnetpp.h"
#include "util.h"

#include <fstream>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <vector>

using namespace omnetpp;
using std::vector;
using std::string;

Define_Module(TrafficGen);

std::once_flag TrafficGen::tFlag;
cXMLElement* TrafficGen::commTable;

void TrafficGen::initialize(){
	this->bsId = par("bsId");
	this->initOffset = par("initOffset");
	this->msId = par("msId");
	this->packetLength = par("packetLength");
	this->fixedPacketSizes = par("fixedPacketSizes");
	this->periodicTraffic = par("periodicTraffic");
	this->tti = par("tti");
	this->d2dActive = par("d2dActive");
	randEng = boost::random::mt19937(getRNG(0)->intRand());
	if(!fixedPacketSizes){
		// Initialize the random distribution for packet sizes
		packetSizeDist = boost::random::normal_distribution<double>(packetLength,25);
	}

	string xmlPath = par("commTable");
	// Only load the table once!
	std::call_once(tFlag,loadComTable,xmlPath);

	vector<StreamDef> parsedStreams(parseCommTable(bsId,msId));
	for(StreamDef& stream:parsedStreams){
		this->streams[stream.streamId] = stream;
		Traffic *msg = new Traffic();
		msg->setPartner(stream.destMsId);
		msg->setStreamId(stream.streamId);
		if(this->periodicTraffic){
			msg->setTrafficType(TrafficType::periodic);
			scheduleAt(initOffset-(2*tti),msg);
		}
		else{
			msg->setTrafficType(TrafficType::gaussian);
			periodDist = boost::random::normal_distribution<double>(stream.period,0.001);
			scheduleAt(initOffset-(2*tti)+periodDist(randEng),msg);
		}
		if(stream.d2d && !d2dActive){
			// D2D needs to be active if we have D2D streams!
			throw std::runtime_error("The Option d2dActive needs to be true when there are D2D streams!");
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
				if(fixedPacketSizes){
					pack->setBitLength(this->packetLength);
				}
				else{
					double rand;
					do{
						rand = packetSizeDist(randEng);
					}
					while(rand<1);
					pack->setBitLength(floor(rand));
				}
				pack->setInterarrival(stream.period);
				pack->setStreamId(stream.streamId);
				pack->setD2d(stream.d2d);
				simtime_t sendTime;
				if(genMsg->getTrafficType()==TrafficType::periodic){
					pack->setTrafficType(TrafficType::periodic);
					sendTime = currTime+stream.period;
				}
				else{
					pack->setTrafficType(TrafficType::gaussian);
					do{
						double rand = periodDist(randEng);
						sendTime = currTime+rand;
					}
					while(sendTime<currTime);
				}
				send(pack,"toMac");
				scheduleAt(sendTime,msg);
				break;
			
		}
	}
	else if(msg->getKind()==MessageType::koidata){
		KoiData *pack = dynamic_cast<KoiData*>(msg);
                /**
		std::cout << "Stream " << pack->getStreamId() << ": " 
			<< "MS " << this->msId 
			<< " got message " << pack->getId() 
			<< " from " << pack->getSrc() 
			<< " to " << pack->getDest() 
			<< std::endl;
                */
		// Write total transmission time in ms to file
		*delays << msId << "\t" << pack->getTotalQueueDelay()*1000 << std::endl;
		delete pack;
	}
	else if(msg->isName("DELAYS_FILE")){
		// Store pointer to delay results file
		ResultFileExchange* ex = dynamic_cast<ResultFileExchange*>(msg);
		delays = ex->getPtr();
		delete ex;
	}
}

vector<TrafficGen::StreamDef> TrafficGen::parseCommTable(int bsId, int msId){
	vector<TrafficGen::StreamDef> parsedStreams;
	string xpath("/root/cell[@id='"
			+std::to_string(bsId)
			+"']/ms[@id='"+std::to_string(msId)+"']");
	cXMLElement *msNode = commTable->getElementByPath(xpath.c_str(),commTable);
	unsigned long sId;
	for(cXMLElement *curr=msNode->getFirstChild();curr!=nullptr;
			curr=curr->getNextSibling()){
		sId = omnetpp::getSimulation()->getUniqueNumber();
		parsedStreams.emplace_back(
			sId,
			std::stoi(curr->getAttribute("cell")),
			std::stoi(curr->getAttribute("destid")),
			std::stod(curr->getAttribute("period")),
			std::stod(curr->getAttribute("deadline")),
			std::stoi(curr->getAttribute("d2d")));
	}
	return parsedStreams;
}

void TrafficGen::loadComTable(const std::string& fpath){
	commTable = omnetpp::getEnvir()->getXMLDocument(fpath.c_str(),nullptr);
}

TrafficGen::~TrafficGen(){
}
