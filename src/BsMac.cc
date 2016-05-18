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
#include "DataPacketBundle_m.h"
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
#include "TransInfoBs_m.h"
#include <algorithm>

using namespace itpp;

Define_Module(BsMac);

void BsMac::initialize()  {
    pos.x = par("xPos");
    pos.y = par("yPos");
    bsId = par("bsId");
    currentChannel = par("currentChannel");
    maxNumberOfNeighbours = par("maxNumberOfNeighbours");
    resourceBlocks = par("resourceBlocks");
    numberOfMobileStations = par("numberOfMobileStations");
    sinr_est = 0;
    transmissionPower = par("transmissionPower");

    initOffset = par("initOffset");
    tti = par("tti");
    epsilon = par("epsilon");
    
    // We need an default SINR to start with, for the first schedule, but 
    // cannot calculate it without the previous schedule, which is non existent.
    // Currently: Random init SINR.
    SINR_ = zeros(numberOfMobileStations,resourceBlocks);
    
    // EESM Beta values for effective SINR
    string eesm_beta = par("eesm_beta");
    eesm_beta_values = vec(eesm_beta);

    // Read Block Error Rate Table
    string bler = par("bler_table");
    blerTable = mat(bler);

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

    //every tti send transmit requests to stream scheduler
    scheduleAt(simTime() + initOffset + 2*tti-epsilon, new cMessage("GEN_TRANSMIT_REQUEST"));
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
		streamQueues[info->getSrc()][info->getDest()];
		delete info;
	}
	else if(msg->getKind()==MessageType::transInfoBs){
		for(int i = 0; i < numberOfMobileStations; ++i)  {
			send(msg->dup(), "toMs", i);
		}
		delete msg;
	}
	else if(msg->isName("POINTER_EXCHANGE2")){
		for(int i = 1; i < numberOfMobileStations; ++i)  {
			send(msg->dup(), "toBsChannel", i);
		}
		delete msg;
	}
	else if(msg->isName("SINR_ESTIMATION")){
		SINR *sinrMessage = (SINR *) msg;
		vec sinr_new;
		int msId = sinrMessage->getMsId();
		for(int i = 0;i < resourceBlocks; i++){
			sinr_new.ins(i,sinrMessage->getSINR(i));
		}
		SINR_.set_row(msId,sinr_new);
		delete msg;
	}
    else if(msg->getKind()==MessageType::streamSched)  {
	    StreamTransSched *sched = dynamic_cast<StreamTransSched*>(msg);
	    DataPacketBundle *packetBundle = new DataPacketBundle("DATA_BUNDLE");

	    int srcMs = sched->getSrc();
	    int destMs = sched->getDest();
            vector<double> sinr_values;
	    sinr_values.push_back(SINR_(destMs,sched->getRb()));
            
            double channel_capacity = getChannelCapacity(sinr_values);
	    /**
            int cqi;
            if(sinr_values.size() > 0){
				cqi = SINR_to_CQI(*(std::min_element(sinr_values.begin(), sinr_values.end())));
			}else{
				cqi = 1;
			}
            **/
	    // For now, only 1 packet will be send per RB in each TTI
            if(channel_capacity > 0)  {
                packetBundle->setPacketsArraySize(1);
		KoiData *packet = dynamic_cast<KoiData*>(
				streamQueues[srcMs][destMs].get(
					sched->getPacketIndex()));
		streamQueues[srcMs][destMs].remove(packet);
		packetBundle->setPackets(0, *packet);
		packetBundle->setRBsArraySize(1);
		packetBundle->setRBs(0,sched->getRb());
		packetBundle->setTransPower(transmissionPower);
		// Set CQI for a fixed value until we decide on how to 
		// compute it
		//packetBundle->setCqi(cqi);
		packetBundle->setCqi(15);
		packetBundle->setMsId(destMs);
		packetBundle->setBsId(packet->getBsId());
		delete packet;

		TransInfoBs *info = new TransInfoBs();
		info->setBsId(bsId);
		info->setPower(transmissionPower);
		info->setRb(sched->getRb());
            
		sendToNeighbourCells(info);
		delete info;
                sendDelayed(packetBundle, epsilon, "toPhy");
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
		for(auto iterSrc=this->streamQueues.begin(); iterSrc!=streamQueues.end();
				++iterSrc){
			if(!iterSrc->second.empty()){
				for(auto iterDest=iterSrc->second.begin(); 
						iterDest!=iterSrc->second.end();
						++iterDest){
					if(!iterDest->second.empty()){
						StreamTransReq *req = new StreamTransReq();
						KoiData *queueHead = dynamic_cast<KoiData*>(iterDest->second.front());
						req->setSrc(iterSrc->first);
						req->setDest(iterDest->first);
						req->setPeriod(queueHead->getInterarrival());
						req->setPackets(&(iterDest->second));
						req->setBs(true);
						send(req,"toScheduler");
					}
				}
			}
		}
		scheduleAt(simTime() + tti-epsilon, msg);
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
        DataPacketBundle *bundle = (DataPacketBundle *) msg;
	KoiData packet;
	for(unsigned int i=0; i<bundle->getPacketsArraySize(); i++){
		packet = bundle->getPackets(i);
		streamQueues[packet.getSrc()][packet.getDest()].insert(packet.dup());
	}
	delete bundle;
    }
    else if(msg->arrivedOn("fromMsMac")){
	if(msg->getKind()==MessageType::transInfoMs){
		// Forward transmission info from local MS to neighbour cells
		sendToNeighbourCells(msg);
		delete msg;
	}
    }
    else if(msg->arrivedOn("fromCell")){
	if(msg->getKind()==MessageType::transInfoMs){
		// Forward transmission info from neighbouring MS to BS 
		// channels via the BS PHY
		send(msg,"toPhy");
	}
    }
}

BsMac::~BsMac()  {
    delete[] msPositions;
    delete msPosUpdateArrived;
    delete neighbourIdMatching;
}
