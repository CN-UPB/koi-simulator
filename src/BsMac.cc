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
#include "BsMsPositions_m.h"
#include "ScheduleInfo_m.h"
#include "StreamInfo_m.h"
#include "StreamTransReq_m.h"
#include "StreamTransSched_m.h"
#include "KoiData_m.h"
#include "ResultFileExchange_m.h"
#include "TransInfo_m.h"
#include "util.h"

#include <algorithm>
#include <fstream>
#include <set>

using std::set;

Define_Module(BsMac);

void BsMac::initialize()  {
    pos.x = par("xPos");
    pos.y = par("yPos");
    bsId = par("bsId");
    maxNumberOfNeighbours = par("maxNumberOfNeighbours");
    numberOfMobileStations = par("numberOfMobileStations");
    transmissionPower = par("transmissionPower");
    initOffset = par("initOffset");
    tti = par("tti");
    epsilon = par("epsilon");

		sinrEstCount = 0;
    
    //find the neighbours and store the pair (bsId, position in data structures) in a map
    cModule *cell = getParentModule()->getParentModule();
    neighbourIdMatching = new NeighbourIdMatching(bsId, maxNumberOfNeighbours, cell);

    //stores the pos of the own ms
		msPositions.reserve(numberOfMobileStations);

    //only send the ms positions to the other cells when all have arrived
		msPosUpdateArrived.resize(numberOfMobileStations,false);

    //exchange the bs positions one time at the sim begin
    scheduleAt(simTime(), new cMessage("BS_POSITION_INIT"));

    // At the end of each tti, send transmit request to scheduler for next 
    // tti.
    scheduleAt(simTime() + initOffset-epsilon, new cMessage("GEN_TRANSMIT_REQUEST"));
		// Result file streams for this cell need to be opened and pointers send
		// out to all local mobile stations.
		scheduleAt(simTime()+epsilon, new cMessage("RESULT_FILES"));

    #ifndef NDEBUG
    scheduleAt(simTime()+initOffset,new cMessage("DEBUG"));
    #endif
}

void BsMac::finish(){
	delays_file.close();
	rate_file.close();
}

/* important: this method does not delete the msg! */
void BsMac::sendToNeighbourCells(cMessage *msg)  {
    sendDelayedToNeighbourCells(msg, SIMTIME_ZERO);
}

/* important: this method does not delete the msg! */
void BsMac::sendDelayedToNeighbourCells(cMessage *msg, const simtime_t& delay){
    //send the msg to the connected neighbour cells
    NeighbourMap *map = neighbourIdMatching->getNeighbourMap();
    for(auto& i:*map)  {
        if(i.first != bsId)  { //skip the own cell
            sendDelayed(msg->dup(), delay, "toCell", (i.second).second);
        }
    }
}

void BsMac::handleMessage(cMessage *msg)  {
	if(msg->getKind()==MessageType::streamInfo){
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
		send(msg->dup(),"toMsMac",0);
		for(int i = 1; i < numberOfMobileStations; ++i)  {
			send(msg->dup(), "toBsChannel", i);
			send(msg->dup(),"toMsMac",i);
		}
		delete msg;
	}
	else if(msg->getKind()==MessageType::sinrEst){
		sinrEstCount++;
		if(sinrEstCount==numberOfMobileStations+1){
			// All local stations have completed their SINR estimates, so the 
			// last TTI's transmission infos can now be cleared from the local 
			// channel.
			send(new cMessage(nullptr,MessageType::cleanupTTI),"toBsChannel",
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
		simtime_t delay = 0.0;
		for(auto& stream:streamQueues){
			list<KoiData*>& currList = stream.second;
			for(auto packetIter = currList.begin(); packetIter!=currList.end();){
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
		PositionExchange *posEx = dynamic_cast<PositionExchange*>(msg);
		msPositions[posEx->getId()] = posEx->getPosition();
		msPosUpdateArrived[posEx->getId()] = true;

		//if all ms positions arrived; send the positions to the other cells
		bool allArrived = std::all_of(msPosUpdateArrived.cbegin(),
				msPosUpdateArrived.cend(),[](bool val) -> bool{return val;});
		if(allArrived)  {
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
			std::fill(msPosUpdateArrived.begin(),msPosUpdateArrived.end(),false);
			writePositions();
		}

		delete msg;
	}
	else if(msg->getKind()==MessageType::scheduleInfo){
		ScheduleInfo *s = dynamic_cast<ScheduleInfo*>(msg);
		comparator = s->getSortfn();
		if(s->getDownStatic()){
			SINR *longTermSINR = new SINR();
			longTermSINR->setKind(MessageType::longTermSinrEst);
			send(longTermSINR,"toBsChannel");
		}
		delete msg;
	}
	else if(msg->getKind()==MessageType::longTermSinrEst){
		send(msg,"toScheduler");
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
			scheduleAt(simTime() + tti, msg);
		}
		else if(msg->isName("DEBUG")){
			// forward debug messages to BS channel
			send(msg->dup(),"toBsChannel",0);
			delete msg;
		}
		else if(msg->isName("RESULT_FILES")){
			// When the number of mobile stations and/or cells increases, opening
			// seperate result files for each MS would hit the maximum number of open
			// files allowed per process. Thus, all per-MS data will now be 
			// written to a single file per cell. As a file can only be opened once,
			// we need to open it centrally here, and then send pointers to the 
			// stream to all local MS.
			// Not all too nice, but the ony possibility right now.
			
			ResultFileExchange *delays = new ResultFileExchange("DELAYS_FILE");
			ResultFileExchange *rates = new ResultFileExchange("RATES_FILE");
			std::string fname("delays-cell-"+std::to_string(bsId));
			delays_file = std::move(getResultFile(fname));
			delays_file << "MS\t" << "Delay" << std::endl;
			delays->setPtr(&delays_file);
			fname = "rates-cell-"+std::to_string(bsId);
			rate_file = std::move(getResultFile(fname));
			rate_file << "MS\t" << "Rate" << std::endl;
			rates->setPtr(&rate_file);
			for(int i = 0; i<numberOfMobileStations; i++){
				send(delays->dup(),"toMsMac",i);
				send(rates->dup(),"toMsMac",i);
			}
			delete msg;
		}
	}
	else if(msg->isName("BS_POSITION_MSG"))  {
		// send the BS position messages BS Channel 0, where it is needed for 
		// channel model initialization.
		send(msg->dup(), "toBsChannel", 0);
	}
	//data packet
	else if(msg->arrivedOn("fromPhy"))  {
		KoiData *data = dynamic_cast<KoiData*>(msg);
		data->setLatestQueueEntry(simTime());
		list<KoiData*>& squeue(streamQueues[data->getStreamId()]);
		auto p = std::lower_bound(squeue.begin(),squeue.end(),data,comparator);
		this->streamQueues[data->getStreamId()].insert(p,data);
	}
}

void BsMac::writePositions(){
	std::string fname("pos-cell-"+std::to_string(bsId));
	std::ofstream fout(getResultFile(fname));
	fout << "Station\t" << "PosX\t" << "PosY" << std::endl; 
	fout << "-1" << "\t" << pos.x << "\t" << pos.y << std::endl;
	for(int i = 0; i<numberOfMobileStations; i++){
		fout << i << "\t" << msPositions[i].x << "\t" << msPositions[i].y << std::endl;
	}
	fout.close();
}

BsMac::~BsMac()  {
    delete neighbourIdMatching;
}
