/*
 * MsMac.cc
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#include "MsMac.h"
#include "PointerExchange_m.h"
#include "PositionExchange_m.h"
#include "KoiData_m.h"
#include "ScheduleInfo_m.h"
#include "TransmitRequest_m.h"
#include "ResultFileExchange_m.h"
#include "Schedule_m.h"
#include "StaticSchedule_m.h"
#include "StreamInfo_m.h"
#include "StreamTransReq_m.h"
#include "StreamTransSched_m.h"
#include "SINR_m.h"
#include "TransInfo_m.h"
#include "util.h"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <random>
#include <set>
#include <utility>

using namespace omnetpp;
using std::set;

Define_Module(MsMac);

inline simtime_t MsMac::positionResendTime()  {
    if(simTime() == initOffset)  {
        return simTime() + (positionResendInterval * tti) - (2 * epsilon); //send with the schedules; one sync point
    }
    else  {
        return simTime() + positionResendInterval * tti;
    }
}

Position MsMac::initMsPositionLinear()  { //for MS Position along a straight road
	Position msPos{};
	ofstream P_MS;
	P_MS.open ("MS_Positions.txt", fstream::app);

        velocity.push_back(30); //velocity in x-direction in m/s
	velocity.push_back(0); //velocity in y-direction in m/s

	msPos.x = 6;
	msPos.y = radius - 9; //for movement of the MS on a straight road
	msPos.z = 0; // MS on road are all the same height
	P_MS << msPos.x << " " << msPos.y << "\n" << std::endl;
	P_MS.close();
	return msPos;
}

Position MsMac::initMsPositionRand(){
	double angle = uniform(0,360);
	double distToBs = uniform(1.1,radius);
	// Convert angle to radians for math library functions
	angle = angle * (M_PI/180);
	Position initPos{};
	initPos.x = initBsPos.x+distToBs*std::cos(angle);
	initPos.y = initBsPos.y+distToBs*std::sin(angle);
	initPos.z = initBsPos.z;
	return initPos;
}

Position MsMac::initMsPositionLine(){
	Position initPos{};
	double msDists = (radius-1.01)/numberOfMobileStations;
	initPos.x = initBsPos.x+1.01+(msId*msDists);
	initPos.y = initBsPos.y;
	initPos.z = initBsPos.z;
	return initPos;
}

void MsMac::initialize()  {
	positionResendInterval = par("positionResendInterval");
	msId = par("msId");
	bsId = par("bsId");
	numberOfMobileStations = par("numberOfMobileStations");
	epsilon = par("epsilon");
	radius = par("radius");
	initBsPos.x = par("initBsXPos");
	initBsPos.y = par("initBsYPos");
	initBsPos.z = par("initBsZPos");
	initOffset = par("initOffset");
	tti = par("tti");
	transmissionPower = par("transmissionPower");
	longTermEst = nullptr;
	chn = nullptr;

	switch((int)par("positioning")){
		case MsMac::Placement::params:
			msPosition.x = par("initMsXPos");
			msPosition.y = par("initMsYPos");
			msPosition.z = par("initMsZPos");
			break;
		case MsMac::Placement::linear:
			msPosition = initMsPositionLinear(); //for MS position along a straight road
			break;
		case MsMac::Placement::uniformRand:
			msPosition = initMsPositionRand();
			break;
		case MsMac::Placement::line:
			msPosition = initMsPositionLine();
			break;
		default:
			// TODO Notify the user that the value for MS positioning
			// is invalid.
			std::cout << "Invalid Ms placement algorithm " << (int) par("positioning") << std::endl;
	}

	// Send MS Position once at the very beginning for cluster generation
	PositionExchange *posEx = new PositionExchange("MS_POS_UPDATE");
	posEx->setId(msId);
	posEx->setPosition(msPosition);
	send(posEx, "toBsMac");
}

void MsMac::finish(){
}

void MsMac::scheduleStatic(cMessage* msg){
	TTISchedule& curr(*staticIter);
	++staticIter;
	if(staticIter==staticSchedule.end()){
		staticIter = staticSchedule.begin();
	}
	TTISchedule& next(*staticIter);

	std::pair<set<int>,set<int>> infos = std::make_pair(set<int>(),set<int>());
	int rate = 0;
	int rbRate = 0;
	simtime_t delay = 0.0;
	list<KoiData*>* bestStream;
	bool assigned = true;
	for(int rb:curr.second){
		if(chn->senseUpSINR(rb,msId,transmissionPower)>=longTermEst->getUp(rb)){
			// Only use this RB if the estimated SINR is above the long term SINR,
			// which allows the initially predicted MCS to be used.
			rbRate = longTermEst->getRUp(rb);
			while(rbRate>0 && assigned){
				// find the best packet
				bestStream = nullptr;
				assigned = false;
				for(auto& streamIter: streamQueues){
					list<KoiData*>& currList = streamIter.second;
					if(!currList.empty() && (bestStream==nullptr 
								|| comparator(currList.front(),bestStream->front()))){
						bestStream = &currList;
					}
				}
				if(bestStream!=nullptr){
					KoiData* packet(bestStream->front());
					bestStream->pop_front();
					packet->setResourceBlock(rb);
					// Send out the best packet of the best stream
					if(rbRate>packet->getBitLength()){
						rate += packet->getBitLength();
						rbRate -= packet->getBitLength();
					}
					else{
						packet->setBitLength(packet->getBitLength()-rbRate);
						KoiData *leftover = new KoiData(*packet);
						leftover->setBitLength(rbRate);
						bestStream->push_front(leftover);
						rate += packet->getBitLength();
						rbRate = 0;
					}
					packet->setTransPower(transmissionPower);
					delay = packet->getTotalQueueDelay() + (simTime() - packet->getLatestQueueEntry());
					packet->setTotalQueueDelay(delay);
					sendDelayed(packet, epsilon, "toPhy");
					if(packet->getMessageDirection()==MessageDirection::up
							|| packet->getMessageDirection()==MessageDirection::d2dUp){
						infos.first.insert(packet->getResourceBlock());
					}
					else{
						infos.second.insert(packet->getResourceBlock());
					}
					assigned = true;
				}
			}
		}
	}
	// Store rate for later evaluation
	*rateFile << msId << "\t" << rate << std::endl;
	// Send out transmission info
	for(auto& rb:infos.first){
		TransInfo *info = new TransInfo();
		info->setBsId(bsId);
		info->setPower(transmissionPower);
		info->setRb(rb);
		info->setMsId(msId);
		info->setMessageDirection(MessageDirection::up);
		send(info,"toBsMac");
	}
	for(auto& rb:infos.second){
		TransInfo *info = new TransInfo();
		info->setBsId(bsId);
		info->setPower(transmissionPower);
		info->setRb(rb);
		info->setMsId(msId);
		info->setMessageDirection(MessageDirection::d2dDown);
		send(info,"toBsMac");
	}
	// Schedule the next transmission opportunity
	int nextTTI = 0;
	if(curr.first==next.first){
		// There is only one TTI in which this MS may transmit, so the next 
		// opportunity is in exactly staticSchedLength many tti
		nextTTI = staticSchedLength;
	}
	else{
		if(next.first < curr.first){
			nextTTI = staticSchedLength - curr.first+next.first;
		}
		else{
			nextTTI = next.first - curr.first;
		}
	}
	scheduleAt(simTime()+(nextTTI*tti),msg);
}

void MsMac::handleMessage(cMessage *msg)  {
	if(msg->getKind()==MessageType::streamSched)  {
		StreamTransSched *schedule = dynamic_cast<StreamTransSched*>(msg);
		std::pair<set<int>,set<int>> infos = std::make_pair(set<int>(),set<int>());

		if(schedule->getSrc() == msId)  {
			KoiData *currPacket = nullptr;
			int rate = 0;
			simtime_t delay = 0.0;
			for(auto& streamIter: streamQueues){
				list<KoiData*>& currList = streamIter.second;
				for(auto packetIter = currList.begin(); 
						packetIter!=currList.end();){
					currPacket = *packetIter;
					if(currPacket->getScheduled()){
						packetIter = currList.erase(packetIter);
						currPacket->setTransPower(transmissionPower);
						// Set CQI for a fixed value until we decide on how to 
						// compute it
						currPacket->setCqi(15);
						delay = currPacket->getTotalQueueDelay() + (simTime() - currPacket->getLatestQueueEntry());
						currPacket->setTotalQueueDelay(delay);
						this->take(currPacket);
						sendDelayed(currPacket, epsilon, "toPhy");
						// Store RB in frequency half dependent sets, so that we send 
						// at most 1 TransInfo for any given resource block, even 
						// if we transmit multiple packets using that block.
						if(currPacket->getMessageDirection()==MessageDirection::up
								|| currPacket->getMessageDirection()==MessageDirection::d2dUp){
							infos.first.insert(currPacket->getResourceBlock());
						}
						else{
							infos.second.insert(currPacket->getResourceBlock());
						}
						rate += currPacket->getBitLength();
					}
					else{
						++packetIter;
					}
				}
			}
		*rateFile << msId << "\t" << rate << std::endl;
		}
		for(auto& rb:infos.first){
			TransInfo *info = new TransInfo();
			info->setBsId(bsId);
			info->setPower(transmissionPower);
			info->setRb(rb);
			info->setMsId(msId);
			info->setMessageDirection(MessageDirection::up);
			send(info,"toBsMac");
		}
		for(auto& rb:infos.second){
			TransInfo *info = new TransInfo();
			info->setBsId(bsId);
			info->setPower(transmissionPower);
			info->setRb(rb);
			info->setMsId(msId);
			info->setMessageDirection(MessageDirection::d2dDown);
			send(info,"toBsMac");
		}
		delete schedule;
	}
	else if(msg->isName("POINTER_EXCHANGE2")){
		PointerExchange *PtrMessage = dynamic_cast<PointerExchange*>(msg);
		chn = PtrMessage->getPtr();
		delete msg;
	}
	else if(msg->isName("SCHEDULE_STATIC")){
		scheduleStatic(msg);
	}
	else if(msg->getKind()==MessageType::transInfo){
		send(msg,"toPhy");
	}
	else if(msg->getKind()==MessageType::scheduleInfo){
		ScheduleInfo *s = dynamic_cast<ScheduleInfo*>(msg);
		comparator = s->getSortfn();
		if(s->getUpStatic()){
			SINR *longTermSINR = new SINR();
			longTermSINR->setKind(MessageType::longTermSinrEst);
			send(longTermSINR,"toPhy");
		}
		else{
			// If the UP direction is not static, generate transmit requests at 
			// at the end of each TTI for the next TTI.
			scheduleAt(initOffset-epsilon, new cMessage("GEN_TRANSMIT_REQUEST"));
		}
		delete msg;
	}
	else if(msg->getKind()==MessageType::staticSchedule){
		StaticSchedule* s = dynamic_cast<StaticSchedule*>(msg);
		staticSchedule = s->getSchedule();
		staticIter = staticSchedule.begin();
		staticSchedLength = s->getScheduleLength();
		TTISchedule& first(staticSchedule.front());
		scheduleAt(initOffset+epsilon+(first.first*tti),new cMessage("SCHEDULE_STATIC"));
		delete s;
	}
	else if(msg->isName("RATES_FILE")){
		// Store pointer to the rate results file
		ResultFileExchange* ex = dynamic_cast<ResultFileExchange*>(msg);
		rateFile = ex->getPtr();
		delete ex;
	}
	else if(msg->isName("DELAYS_FILE")){
		// Forward pointer to delays results file to the traffic generator
		send(msg,"toApp");
	}
	else if(msg->isName("GEN_TRANSMIT_REQUEST"))  {
		// Send requests for each stream originating from this MS to the 
		// scheduler if that stream has packets.
		for(auto iter=this->streamQueues.begin(); iter!=streamQueues.end();
				++iter){
			if(!iter->second.empty()){
				StreamTransReq *req = new StreamTransReq();
				KoiData *queueHead = dynamic_cast<KoiData*>(iter->second.front());
				req->setSrc(this->msId);
				req->setDest(queueHead->getDest());
				req->setStreamId(iter->first);
				req->setPeriod(queueHead->getInterarrival());
				req->setPackets(&(iter->second));
				req->setRequestOrigin(msId);
				if(queueHead->getD2d()){
					req->setMessageDirection(MessageDirection::d2d);
				}
				else{
					req->setMessageDirection(MessageDirection::up);
				}
				send(req,"toScheduler");
			}
		}
		scheduleAt(simTime() + tti, msg);
	}
	else if(msg->getKind()==MessageType::sinrEst)  {
		SINR* est = dynamic_cast<SINR*>(msg);
		// Forward estimates to BS Mac
		send(est,"toBsMac");
	}
	else if(msg->getKind()==MessageType::longTermSinrEst)  {
		// Store long term SINR estimate
		SINR* est = dynamic_cast<SINR*>(msg);
		if(longTermEst!=nullptr){
			delete longTermEst;
		}
		longTermEst = est;
		// Forward long term estimate to scheduler
		send(est->dup(),"toScheduler");
	}
	else if(msg->isName("MCS_FILE")){
		send(msg,"toPhy");
	}
	else if(msg->arrivedOn("fromApp"))  {
		// Packet arrived for sending from traffic generator
		switch(msg->getKind()){
			case MessageType::streamInfo:{
				// Add queue for the new stream
				StreamInfo *tmp = dynamic_cast<StreamInfo*>(msg);
				this->streamQueues[tmp->getStreamId()];
				send(tmp->dup(),"toScheduler");
				send(tmp->dup(),"toBsMac");
				delete tmp;
			} break;
			case MessageType::koidata:{
				KoiData *data = dynamic_cast<KoiData*>(msg);
				data->setLatestQueueEntry(simTime());
				list<KoiData*>& squeue(streamQueues[data->getStreamId()]);
				auto p = std::lower_bound(squeue.begin(),squeue.end(),data,comparator);
				this->streamQueues[data->getStreamId()].insert(p,data);
			} break;
		}
	}
	else if(msg->arrivedOn("fromPhy"))  {
		// Unpack the data bundle and forward data packets to app
		if(msg->getKind()==MessageType::koidata){
			send(msg,"toApp");
		}
	}
	else{
		std::cout << "Received unhandled message: " << msg->getKind() << std::endl;
	}
}
