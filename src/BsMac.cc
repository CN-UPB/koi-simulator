/*
 * BsMac.cc
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 * Last edited: May 20, 2014
 *      Author: Thomas Prinz
 */

#include "BsMac.h"
#include "SINR_m.h"
#include "PositionExchange_m.h"
#include "ChannelExchange_m.h"
#include "BsMsPositions_m.h"
#include "VisibilityRegionMessage_m.h"
#include "PointerExchange_m.h"
#include "ClusterMessage_m.h"
#include "StreamInfo_m.h"
#include "StreamTransReq_m.h"
#include "StreamTransSched_m.h"
#include "KoiData_m.h"
#include "TransInfo_m.h"
#include <algorithm>
#include <fstream>
#include <set>

using namespace itpp;
using std::set;

Define_Module(BsMac);

void BsMac::initialize()  {
    pos.x = par("xPos");
    pos.y = par("yPos");
    bsId = par("bsId");
    currentChannel = par("currentChannel");
    maxNumberOfNeighbours = par("maxNumberOfNeighbours");
    resourceBlocks = par("downResourceBlocks");
    numberOfMobileStations = par("numberOfMobileStations");
    sinr_est = 0;
    transmissionPower = par("transmissionPower");
    initOffset = par("initOffset");
    tti = par("tti");
    epsilon = par("epsilon");
		avgRatePerStation = registerSignal("avgRatePerStation");

		sinrEstCount = 0;
    
    //find the neighbours and store the pair (bsId, position in data structures) in a map
    cModule *cell = getParentModule()->getParentModule();
    neighbourIdMatching = new NeighbourIdMatching(bsId, maxNumberOfNeighbours, cell);

    /*ev << "Neigbour map---------------" << endl;
    NeighbourIdMap *map = neighbourIdMatching->getNeighbourIdMap();
    for(NeighbourIdMap::iterator i = map->begin(); i != map->end(); i++)  {
        ev << "BS " << i->first << ": Real pos " << i->second << endl;
    }
    ev << "-------------------------------------------------------" << endl;*/

    ownDataStrId = neighbourIdMatching->getDataStrId(bsId);

    //stores the pos of the own ms
    msPositions = new Position[numberOfMobileStations]; //TODO change to dynamic size per bs

    //only send the ms positions to the other cells when all have arrived
    msPosUpdateArrived = new vector<bool>(numberOfMobileStations);
    for(int i = 0; i < numberOfMobileStations; ++i)
        msPosUpdateArrived->at(i) = false;

    //exchange the bs positions one time at the sim begin
    scheduleAt(simTime(), new cMessage("BS_POSITION_INIT"));

    // At the end of each tti, send transmit request to scheduler for next 
    // tti.
    scheduleAt(simTime() + initOffset-epsilon, new cMessage("GEN_TRANSMIT_REQUEST"));

    #ifndef NDEBUG
    scheduleAt(simTime()+initOffset,new cMessage("DEBUG"));
    #endif
}

/* important: this method does not delete the msg! */
void BsMac::sendToNeighbourCells(cMessage *msg)  {
    sendDelayedToNeighbourCells(msg, SIMTIME_ZERO);
}

/* important: this method does not delete the msg! */
void BsMac::sendDelayedToNeighbourCells(cMessage *msg, simtime_t delay)  {
    //send the msg to the connected neighbour cells
    NeighbourMap *map = neighbourIdMatching->getNeighbourMap();
    for(NeighbourMap::iterator i = map->begin(); i != map->end(); i++)  {
        if(i->first != bsId)  { //skip the own cell
            sendDelayed(msg->dup(), delay, "toCell", (i->second).second);
        }
    }
}

void BsMac::handleMessage(cMessage *msg)  {
	if(msg->isName("CHANNEL_INFO")){
		if(msg->getKind() % 2 == 1){
			msg->setKind(msg->getKind() + 1);
			sendToNeighbourCells(msg->dup());
			send(msg->dup(), "toBsChannel", 0);
			delete msg;
		}else if(msg->getKind() % 2 == 0){
			send(msg->dup(), "toBsChannel", 0);
			delete msg;
		}
	}
	else if(msg->getKind()==MessageType::streamInfo){
		StreamInfo *info = dynamic_cast<StreamInfo*>(msg);
		// Add a packet queue for this stream, using associative containers 
		// automatic insertion of new entries when subscripting.
		streamQueues[info->getStreamId()];
		delete info;
	}
	else if(msg->getKind()==MessageType::transInfo){
		// Route transmission information according to origin
		TransInfo *trans = dynamic_cast<TransInfo*>(msg);
		if(trans->arrivedOn("fromMsMac")){
			// The message is transmission information from 
			// one of the local mobile stations. Forward to 
			// neighbours.
			sendToNeighbourCells(trans);
		}
		else if(trans->arrivedOn("fromCell")){
			// Message arrived from neighbouring cell
			send(trans->dup(),"toBsChannel",0);
		}
		delete msg;
	}
	else if(msg->isName("POINTER_EXCHANGE2")){
		for(int i = 1; i < numberOfMobileStations; ++i)  {
			send(msg->dup(), "toBsChannel", i);
		}
		delete msg;
	}
	else if(msg->getKind()==MessageType::sinrEst){
		SINR *sinrMessage = (SINR *) msg;
		sinrEstCount++;
		if(sinrEstCount==numberOfMobileStations+1){
			// All local stations have completed their SINR estimates, so the 
			// last TTI's transmission infos can now be cleared from the local 
			// channel.
			send(new cMessage(nullptr,MessageType::clearTransInfo),"toBsChannel",
					0);
			sinrEstCount = 0;
		}
		// Provide the estimates to the scheduler, too
		send(msg,"toScheduler");
	}
	else if(msg->getKind()==MessageType::streamSched)  {
		StreamTransSched *sched = dynamic_cast<StreamTransSched*>(msg);
		KoiData *currPacket = nullptr;
		set<int> infos;
		int rate = 0;
		for(auto streamIter = streamQueues.begin();
				streamIter!=streamQueues.end(); ++streamIter){
			list<KoiData*>& currList = streamIter->second;
			for(auto packetIter = currList.begin(); packetIter!=currList.end();){
				currPacket = *packetIter;
				if(currPacket->getScheduled()){
					packetIter = currList.erase(packetIter);
					currPacket->setTransPower(transmissionPower);
					// Set CQI for a fixed value until we decide on how to 
					// compute it
					currPacket->setCqi(15);
					sendDelayed(currPacket, epsilon, "toPhy");

					// Store used resource block in set for later generation of
					// TransInfo messages. That way, we can make sure that only one 
					// TransInfo is send out per used RB.
					infos.insert(currPacket->getResourceBlock());
					rate += currPacket->getBitLength();
				}
				else{
					++packetIter;
				}
			}
		}
		emit(avgRatePerStation,rate);
		// Send out exactly one TransInfo per used resource block
		for(auto& rb:infos){
			TransInfo *info = new TransInfo();
			info->setBsId(bsId);
			info->setPower(transmissionPower);
			info->setRb(rb);
			// It is the BS itself sending, not a MS, which we indicate
			// with an index of -1 for the MS
			info->setMsId(-1);
			info->setMessageDirection(MessageDirection::down);
			sendToNeighbourCells(info);
			delete info;
		}
		delete sched;
	}
	else if(msg->isName("BS_MS_POSITIONS"))  {
		//send the ms positions to the own bs channels
		//ev << "Forwarding the ms positions to the bs channels" << endl;
		for(int i = 0; i < numberOfMobileStations; ++i)  {
			send(msg->dup(), "toBsChannel", i);
		}
		delete msg;
	}
	else if(msg->isName("MS_POS_UPDATE"))  {
		//save the position update of the mobile stations
		PositionExchange *posEx = (PositionExchange *) msg;
		msPositions[posEx->getId()] = posEx->getPosition();
		msPosUpdateArrived->at(posEx->getId()) = true;

		//if all ms positions arrived; send the positions to the other cells
		bool allArrived = true;
		for(int i = 0; i < numberOfMobileStations; i++)  {
			if(!msPosUpdateArrived->at(i))
				allArrived = false;
		}
		if(allArrived)  {
			//ev << "BS " << bsId << " all ms positions arrived!" << endl;

			BsMsPositions *msPos = new BsMsPositions("BS_MS_POSITIONS");
			msPos->setBsId(bsId);
			msPos->setPositionsArraySize(numberOfMobileStations);
			for(int i = 0; i < numberOfMobileStations; ++i)  {
				msPos->setPositions(i, msPositions[i]);
			}
			//send the ms positions to the other connected cells
			sendToNeighbourCells(msPos);
			//send the ms positions to the own bs channels
			for(int i = 0; i < numberOfMobileStations; ++i)  {
				send(msPos->dup(), "toBsChannel", i);
			}
			delete msPos;

			//reset the arrived vector
			for(int i = 0; i < numberOfMobileStations; ++i)
				msPosUpdateArrived->at(i) = false;
			writePositions();
		}

		delete msg;
	}
	else if(msg->isSelfMessage())  {
		if(msg->isName("BS_POSITION_INIT"))  {
			//exchange the bs positions one time at the sim begin
			PositionExchange *bsPos = new PositionExchange("BS_POSITION_MSG");
			bsPos->setId(bsId);
			bsPos->setPosition(pos);
			//send the bs position to the other cells
			send(bsPos->dup(), "toBsChannel", 0);
			sendToNeighbourCells(bsPos);
			/* send the position also to the own ms
			 * also possible to do with the config
			 * but i think this is more flexible */
			send(bsPos, "toPhy");
			delete msg;
		}
		else if(msg->isName("GEN_TRANSMIT_REQUEST"))  {
			// Send requests for each stream to the 
			// scheduler if that stream has packets.
			for(auto iter=this->streamQueues.begin(); iter!=streamQueues.end();
					++iter){
				if(!iter->second.empty()){
					StreamTransReq *req = new StreamTransReq();
					KoiData *queueHead = dynamic_cast<KoiData*>(iter->second.front());
					req->setSrc(queueHead->getSrc());
					req->setDest(queueHead->getDest());
					req->setStreamId(iter->first);
					req->setPeriod(queueHead->getInterarrival());
					req->setPackets(&(iter->second));
					req->setMessageDirection(MessageDirection::down);
					req->setRequestOrigin(-1);
					send(req,"toScheduler");
				}
			}
			scheduleAt(simTime() + tti-epsilon, msg);
		}
		else if(msg->isName("DEBUG")){
			// forward debug messages to BS channel
			send(msg->dup(),"toBsChannel",0);
			delete msg;
		}
	}
	else if(msg->isName("BS_POSITION_MSG"))  {

		//DEBUG: send the BS position messages also to BS Channel 0.
		send(msg->dup(), "toBsChannel", 0);

		//forward bs position msg from the other cells to the ms
		send(msg, "toPhy");
	}
	//data packet
	else if(msg->arrivedOn("fromPhy"))  {
		KoiData *packet = (KoiData *) msg;
		streamQueues[packet->getStreamId()].push_back(packet);
	}
}

void BsMac::writePositions(){
	int run = std::stoi(ev.getConfig()->substituteVariables("${runnumber}"));
	std::fstream fout;
	fout.open("./results/run_"+std::to_string(run)+"_pos_cell_"+std::to_string(bsId)+".dat",std::fstream::out);
	fout << "[" << "BS_" << bsId << "]" << std::endl
		<< "posX=" << pos.x << std::endl
		<< "posY=" << pos.y << std::endl;
	for(int i = 0; i<numberOfMobileStations; i++){
		fout << "[MS_" << i << "]" << std::endl
			<< "posX=" << msPositions[i].x << std::endl
			<< "posY=" << msPositions[i].y << std::endl
			<< "bs=" << bsId << std::endl;
	}
	fout.close();
}

BsMac::~BsMac()  {
    delete[] msPositions;
    delete msPosUpdateArrived;
    delete neighbourIdMatching;
}
