/*
 * BsMac.cc
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 * Last edited: May 20, 2014
 *      Author: Thomas Prinz
 */

#include "SINR_m.h"
#include "BsMac.h"
#include "DataPacket_m.h"
#include "DataPacketBundle_m.h"
#include "PositionExchange_m.h"
#include "ChannelExchange_m.h"
#include "BsMsPositions_m.h"
#include "TransmitRequest_m.h"
#include "VisibilityRegionMessage_m.h"
#include "PointerExchange_m.h"
#include "ClusterMessage_m.h"
#include "Schedule_m.h"
#include <algorithm>

using namespace itpp;

Define_Module(BsMac);

BsMac::BsMac()  {
    packetQueue = NULL;
}

void BsMac::initialize()  {
    pos.x = par("xPos");
    pos.y = par("yPos");
    bsId = par("bsId");
    currentChannel = par("currentChannel");
    maxNumberOfNeighbours = par("maxNumberOfNeighbours");
    upResBlocks = par("upResourceBlocks");
    downResBlocks = par("downResourceBlocks");
    numberOfMobileStations = par("numberOfMobileStations");
    sinr_est = 0;

    initOffset = par("initOffset");
    tti = par("tti");
    epsilon = par("epsilon");
    
    scheduler = new proportionalFairScheduler();
    scheduler->init(this);
    
    // We need an default SINR to start with, for the first schedule, but 
    // cannot calculate it without the previous schedule, which is non existent.
    // Currently: Random init SINR.
    SINR_ = zeros(numberOfMobileStations,upResBlocks);
    
    // EESM Beta values for effective SINR
	string eesm_beta = par("eesm_beta");
	eesm_beta_values = vec(eesm_beta);
	
	// Read Block Error Rate Table
	string bler = par("bler_table");
	blerTable = mat(bler);
	
	// Down/Up Configuration 0
	downToUpPeriodicity = par("downToUpPeriodicity");
	currentPeriodicity = par("currentPeriodicity");

    //Every App gets an own sending queue
    packetQueue = new cQueue[numberOfMobileStations];

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

    //for the transmit request of the ms
    transmitRequests = new vector<int>(numberOfMobileStations);
    for(int i = 0; i < numberOfMobileStations; ++i)
        transmitRequests->at(i) = 0;

    //only send the ms positions to the other cells when all have arrived
    msPosUpdateArrived = new vector<bool>(numberOfMobileStations);
    for(int i = 0; i < numberOfMobileStations; ++i)
        msPosUpdateArrived->at(i) = false;

    //every tti estimate a new schedule and send it to the MSs and the dataChn of the BS
    cMessage *scheduleMsg = new cMessage("CALC_NEW_SCHEDULE");
    scheduleMsg->setEventDuration(tti); //allows horizon to parallelize other events
    scheduleAt(simTime() + initOffset, scheduleMsg);

    //init the down schedule from bs to ms
    downSchedule = new vector<int>(downResBlocks);

    //exchange the bs positions one time at the sim begin
    scheduleAt(simTime(), new cMessage("BS_POSITION_INIT"));

    //exchange the bs channels one time at the sim begin
    scheduleAt(simTime(), new cMessage("BS_CHANNEL_INIT"));
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

inline simtime_t BsMac::scheduleTime()  {
    return simTime() + tti;
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
		for(int i = 0;i < downResBlocks; i++){
			//std::cout << "SINR for RB " << i << " from MS " << msId << " " << sinrMessage->getSINR(i) << std::endl;
			sinr_new.ins(i,sinrMessage->getSINR(i));
		}
		SINR_.set_row(msId,sinr_new);
		if((sinr_est++) == numberOfMobileStations){
			sinr_est = 0;
			scheduler->setSINR(SINR_);
			scheduler->updateSINR();
			//std::cout << "New SINR from Message: " << SINR_ << std::endl;
		}
		delete msg;
	}
    else if(msg->isName("SCHEDULE"))  {
        Schedule *schedule = (Schedule *) msg;
        //only forward schedules with the same channels
        if(schedule->getChannel() == currentChannel)  {
            //Forward schedule to the own bs channels
            for(int i = 0; i < numberOfMobileStations; ++i)  {
                send(msg->dup(), "toBsChannel", i);
            }
        }
        delete schedule;
    }
    else if(msg->isName("EXEC_DOWN_SCHEDULE"))  {
		// NEW:
		// Get the Number of RBs, that is assigned to this MS and the minimal SINR
        int numberOfRB[numberOfMobileStations];
        vector<vector<double>> sinrValues(numberOfMobileStations,vector<double>());
        vector<vector<double>> RBs(numberOfMobileStations,vector<double>());
        for(int i = 0; i < numberOfMobileStations; ++i){
			numberOfRB[i] = 0;
		}
        for(int i = 0; i < downResBlocks; ++i)  {
            int toMsId = downSchedule->at(i);
            if(toMsId == -1)
                break; //only empty slots can follow
            numberOfRB[toMsId]++;
            sinrValues[toMsId].push_back(SINR_(toMsId,i));
            RBs[toMsId].push_back(i);
        }
        // Get "Sending Capacity" for each MS
        int sendingCapacity[numberOfMobileStations];
        int cqi[numberOfMobileStations];
        
        for(int i = 0; i < numberOfMobileStations; i++){
			sendingCapacity[i] = 0;
			sendingCapacity[i] = getChannelCapacity(sinrValues[i]);
			if(sinrValues[i].size() > 0){
				cqi[i] = SINR_to_CQI(*(std::min_element(sinrValues[i].begin(), sinrValues[i].end())));
			}else{
				cqi[i] = 1;
			}
		}
		
	for(int i = 0; i < numberOfMobileStations; i++){
			//std::cout << "Channel Capacity in bits/TTI for MS " << i << " is " << sendingCapacity[i] << std::endl;
		}
		
		//how much packets to each ms; for the bundle packets
        int packetsToSend[numberOfMobileStations];
        for(int i = 0; i < numberOfMobileStations; ++i)
            packetsToSend[i] = 0;
		
		// Get more Packets, as long as Capacity is left:
		for(int i = 0; i < numberOfMobileStations; i++){
			int toMsId = i;
			while(sendingCapacity[i] > 0){
				//cout << "current capacity: " << sendingCapacity[i] << endl;
				DataPacket *packet = (DataPacket *) packetQueue[toMsId].get(packetsToSend[toMsId]);
				if(packet == NULL){
					break;
				}
				//packet->setResourceBlock(i); //set the resource block of the packet
				if((sendingCapacity[i] - packet->getBitLength()) > 0){
					sendingCapacity[i] = sendingCapacity[i] - packet->getBitLength();
				}else{
					//std::cout << "Last one cut off at: " << packet->getBitLength() - sendingCapacity[i] << " Bits" << std::endl;
					DataPacket *packet_new = new DataPacket("DATA");
					packet_new->setBitLength(packet->getBitLength() - sendingCapacity[i]);
					packet_new->setBsId(packet->getBsId());
					packet_new->setMsId(packet->getMsId());
					packetQueue[toMsId].insert(packet_new);
					packet->setBitLength(sendingCapacity[i]);
					sendingCapacity[i] = 0;
				}

				packetsToSend[toMsId]++;
			}
		}
		
		for(int i = 0; i < numberOfMobileStations; ++i)  {
            		if(packetsToSend[i] <= 0)
                		continue;

            		//ev << "BS " << bsId << " sending " << packetsToSend[i] << " packets to ms " << i << endl;
            		DataPacketBundle *bundle = new DataPacketBundle("DATA_BUNDLE");
            		bundle->setMsId(i);
	            	bundle->setBsId(bsId);
            		bundle->setCqi(cqi[i]);
            		//cout << "Setting cqi for this bundle to: " << cqi[i] << endl;
            		bundle->setPacketsArraySize(packetsToSend[i]);
            		bundle->setRBsArraySize(RBs[i].size());
            		for(int j = 0; j < packetsToSend[i]; ++j)  {
                		DataPacket *packet = (DataPacket *) packetQueue[i].pop();
                		bundle->setPackets(j, *packet);
                		delete packet;
            		}
            		for(uint j = 0; j < RBs[i].size(); j++){
				bundle->setRBs(j,RBs[i][j]);
			}
            	sendDelayed(bundle, epsilon, "toPhy");
        	}
		
		/*
        //how much packets to each ms; for the bundle packets
        int packetsToSend[numberOfMobileStations];
        for(int i = 0; i < numberOfMobileStations; ++i)
            packetsToSend[i] = 0;

        //send the packets from the down schedule to the mobile stations
        for(int i = 0; i < downResBlocks; ++i)  {
            int toMsId = downSchedule->at(i);
            if(toMsId == -1)
                break; //only empty slots can follow

            DataPacket *packet = (DataPacket *) packetQueue[toMsId].get(packetsToSend[toMsId]);
            packet->setResourceBlock(i); //set the resource block of the packet

            packetsToSend[toMsId]++;
        }
		*/
		/*
        for(int i = 0; i < numberOfMobileStations; ++i)  {
            if(packetsToSend[i] <= 0)
                continue;

            //ev << "BS " << bsId << " sending " << packetsToSend[i] << " packets to ms " << i << endl;
            DataPacketBundle *bundle = new DataPacketBundle("DATA_BUNDLE");
            bundle->setMsId(i);
            bundle->setBsId(i);
            bundle->setPacketsArraySize(packetsToSend[i]);
            for(int j = 0; j < packetsToSend[i]; ++j)  {
                DataPacket *packet = (DataPacket *) packetQueue[i].pop();
                bundle->setPackets(j, *packet);
                delete packet;
            }
            sendDelayed(bundle, epsilon, "toPhy");
        }
        */
        delete msg;
    }
    else if(msg->isName("TRANSMIT_REQUEST"))  {
        //save the transmission request
        TransmitRequest *transReq = (TransmitRequest *) msg;
        int fromMsId = transReq->getId();
        transmitRequests->at(fromMsId) = transReq->getPacketCount();
        //ev << "BS Mac" << bsId << ": Received transmit reqeuest from MS " << fromMsId << endl;
        delete msg;
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
        }
        if(msg->isName("BS_CHANNEL_INIT"))  {
            ChannelExchange *chnEx = new ChannelExchange("BS_CHANNEL_UPDATE");
            chnEx->setId(bsId);
            chnEx->setChannel(currentChannel);
            //send the bs channel to the other cells
            sendToNeighbourCells(chnEx);
            // send the channel also to the own ms
            send(chnEx, "toPhy");
        }
        else if(msg->isName("CALC_NEW_SCHEDULE"))  {
            //calculate the schedule for the next TTI
            calcNewSchedule();
            scheduleAt(scheduleTime(), msg); //next schedule calculation
        }
    }
    else if(msg->isName("BS_POSITION_MSG"))  {
        
        //DEBUG: send the BS position messages also to BS Channel 0.
        send(msg->dup(), "toBsChannel", 0);
        
        //forward bs position msg from the other cells to the ms
        send(msg, "toPhy");
    }
    else if(msg->isName("BS_CHANNEL_UPDATE"))  {
        //forward bs channel update from the other cells to the ms
        send(msg, "toPhy");
    }
    //data packet
    else if(msg->arrivedOn("fromApp"))  {
        //insert msg from app into sending queue
        DataPacket *packet = (DataPacket *) msg;
        int msId = packet->getMsId();
        //ev << "Message arrived from App. Stored in Queue[" << msId << "]" << endl;
        if(packetQueue[msId].getLength() < 25){
			packetQueue[msId].insert(packet);
		}else{
			delete msg;
		}
        //ev << "#Queue[" << msId << "]: " << packetQueue[msId].length() << endl;
    }
    //data packet
    else if(msg->arrivedOn("fromPhy"))  {
        DataPacketBundle *bundle = (DataPacketBundle *) msg;
        int msId = bundle->getMsId();
        //ev << "Forwarding Packet to App " << msId << endl;
        send(msg, "toApp", msId);
    }
}

simtime_t BsMac::getProcessingDelay(cMessage *msg)  {
    if(msg->isName("CALC_NEW_SCHEDULE"))
        return tti - 3 * epsilon;
    else
        return 0;
}

// Returns the Data Direction for the current TTI
// Due to finite speed of light and antenna switching time it is necessary
// to have a special guard TTI when switching from Downlink to Uplink
DATADIRECTION BsMac::getDataDirection(){
	switch(downToUpPeriodicity){
		case 0:
			switch(currentPeriodicity){
				case 0:
					return DOWN;
				case 1:
					return SPECIAL;
				case 2:
					return UP;
				case 3:
					return UP;
				case 4:
					return UP;
				case 5:
					return DOWN;
				case 6:
					return SPECIAL;
				case 7:
					return UP;
				case 8:
					return UP;
				case 9:
					return UP;
			}
			break;
		case 1:
			switch(currentPeriodicity){
				case 0:
					return DOWN;
				case 1:
					return SPECIAL;
				case 2:
					return UP;
				case 3:
					return UP;
				case 4:
					return DOWN;
				case 5:
					return DOWN;
				case 6:
					return SPECIAL;
				case 7:
					return UP;
				case 8:
					return UP;
				case 9:
					return DOWN;
			}
			break;
		case 2:
			switch(currentPeriodicity){
				case 0:
					return DOWN;
				case 1:
					return SPECIAL;
				case 2:
					return UP;
				case 3:
					return DOWN;
				case 4:
					return DOWN;
				case 5:
					return DOWN;
				case 6:
					return SPECIAL;
				case 7:
					return UP;
				case 8:
					return DOWN;
				case 9:
					return DOWN;
			}
			break;
		case 3:
			switch(currentPeriodicity){
				case 0:
					return DOWN;
				case 1:
					return SPECIAL;
				case 2:
					return UP;
				case 3:
					return UP;
				case 4:
					return UP;
				case 5:
					return DOWN;
				case 6:
					return DOWN;
				case 7:
					return DOWN;
				case 8:
					return DOWN;
				case 9:
					return DOWN;
			}
			break;
		case 4:
			switch(currentPeriodicity){
				case 0:
					return DOWN;
				case 1:
					return SPECIAL;
				case 2:
					return UP;
				case 3:
					return UP;
				case 4:
					return DOWN;
				case 5:
					return DOWN;
				case 6:
					return DOWN;
				case 7:
					return DOWN;
				case 8:
					return DOWN;
				case 9:
					return DOWN;
			}
			break;
		case 5:
			switch(currentPeriodicity){
				case 0:
					return DOWN;
				case 1:
					return SPECIAL;
				case 2:
					return UP;
				case 3:
					return DOWN;
				case 4:
					return DOWN;
				case 5:
					return DOWN;
				case 6:
					return DOWN;
				case 7:
					return DOWN;
				case 8:
					return DOWN;
				case 9:
					return DOWN;
			}
			break;
		case 6:
			switch(currentPeriodicity){
				case 0:
					return DOWN;
				case 1:
					return SPECIAL;
				case 2:
					return UP;
				case 3:
					return UP;
				case 4:
					return UP;
				case 5:
					return DOWN;
				case 6:
					return SPECIAL;
				case 7:
					return UP;
				case 8:
					return UP;
				case 9:
					return DOWN;
			}
			break;
		case 7:
			// For Testing Purpose
			switch(currentPeriodicity){
				case 0:
					return DOWN;
				case 1:
					return DOWN;
				case 2:
					return DOWN;
				case 3:
					return DOWN;
				case 4:
					return DOWN;
				case 5:
					return DOWN;
				case 6:
					return DOWN;
				case 7:
					return DOWN;
				case 8:
					return DOWN;
				case 9:
					return DOWN;
			}
			break;
	}
	// Default.
	return SPECIAL;
}

void BsMac::calcNewSchedule()  {
	//std::cout << "Calculation Schedule with Periodicity: " << currentPeriodicity << std::endl;
    //from the bs to ms
    vector<int> packetCount(numberOfMobileStations);
    //vector that stores the ids of the sending apps; dynamic size numberOfSendingApps
    vector<int> sendingApps(numberOfMobileStations);
    
    int numberOfPacketsToSend = 0;
    int numberOfSendingApps = 0;
    //save the length of the queues; and the total number of packets
    for(int i = 0; i < numberOfMobileStations; ++i)  {
        int length = packetQueue[i].length();
        numberOfPacketsToSend += length;
        packetCount.at(i) = length;

        if(length != 0)  {
            sendingApps.at(numberOfSendingApps) = i;
            numberOfSendingApps++;
        }
    }
    
	//schedule msg for the ms and the bs channels
	Schedule *schedule = new Schedule("SCHEDULE");
	schedule->setBsId(bsId);
	schedule->setChannel(currentChannel);
    
    if(getDataDirection() == DOWN){
		//cout << "Calculation of Schedule. We are currently in Down Schedule." << endl;
		schedule->setScheduleDirection(1);
		
		*downSchedule = scheduler->calcDownSchedule(packetQueue);

		// Just mark all slots as free. We are in Down Schedule, so 
		// there is no Uplink scheduling.
		schedule->setUpScheduleArraySize(upResBlocks);
		schedule->setPowerAdaptationArraySize(upResBlocks);
		for(int i = 0; i < upResBlocks; ++i)  {
			schedule->setUpSchedule(i, -1); //init the slot as free
			schedule->setPowerAdaptation(i, 1 / upResBlocks); //power equal for all RBs
		}
	}else if(getDataDirection() == UP){
		//cout << "Calculation of Schedule. We are currently in Up Schedule." << endl;
		schedule->setScheduleDirection(0);
		//std::cout << "UP Scheduling! " << bsId << std::endl;
		// Just mark all slots as free. We are in Up Schedule, so 
		// there is no Downlink scheduling.
		for(int i = 0; i < downResBlocks; ++i)  {
			downSchedule->at(i) = -1; //mark the slot as free
		}

		//from the ms to bs
		int numberOfSendingMs = 0;
		int numberOfPacketsToReceive = 0;
		//vector that stores the ids of the sending ms; dynamic size numberOfSendingMs
		vector<int> sendingMs(numberOfMobileStations);
		for(int i = 0; i < numberOfMobileStations; ++i)  {
			numberOfPacketsToReceive += transmitRequests->at(i);
			if(transmitRequests->at(i) != 0)  {
				sendingMs.at(numberOfSendingMs) = i;
				numberOfSendingMs++;
			}
		}
		
		schedule = scheduler->calcUpSchedule(packetQueue,transmitRequests);
		schedule->setBsId(bsId);
		schedule->setChannel(currentChannel);
	}else{
		//cout << "Calculation of Schedule. We are currently in Guard Frame." << endl;
		// Just mark all slots as free.
		//std::cout << "Guard Schedule Computation (No Scheduling): " << std::endl;
		schedule->setScheduleDirection(2);
		for(int i = 0; i < downResBlocks; ++i)  {
			downSchedule->at(i) = -1; //mark the slot as free
		}
		
		schedule->setUpScheduleArraySize(upResBlocks);
		schedule->setPowerAdaptationArraySize(upResBlocks);
		for(int i = 0; i < upResBlocks; ++i)  {
			schedule->setUpSchedule(i, -1); //init the slot as free
			schedule->setPowerAdaptation(i, 1.0 / upResBlocks); //power equal for all RBs
		}
	}
		
	//send the schedule to the ms
	for(int i = 0; i < numberOfMobileStations; ++i)  {
		sendDelayed(schedule->dup(), tti - 2 * epsilon, "toMsMac", i);
		//send(schedule->dup(), "toMsMac", i);
	}

	//send the schedule to the bs channels
	for(int i = 0; i < numberOfMobileStations; ++i)  {
		sendDelayed(schedule->dup(), tti - 2 * epsilon, "toBsChannel", i);
		//send(schedule->dup(), "toBsChannel", i);
	}

	//send the schedule to the other cells
	sendDelayedToNeighbourCells(schedule, tti - 2 * epsilon);

	//schedule the down schedule
	cMessage *msg = new cMessage("EXEC_DOWN_SCHEDULE");
	scheduleAt(scheduleTime() - 2 * epsilon, msg);

	delete schedule;
	currentPeriodicity = (currentPeriodicity + 1) % 10;

    /*ev << "BsId: " << bsId << " Down Schedule--------------------" << endl;
    for(int i = 0; i < downResBlocks; ++i)  {
        ev << "Slot " << i <<": " << downSchedule->at(i) << endl;
    }
    ev << "Up Schedule-----------------------" << endl;
    for(int i = 0; i < upResBlocks; ++i)  {
        ev << "Slot " << i <<": " << schedule->getUpSchedule(i) << endl;
    }
    ev << "----------------------------------" << endl;*/
}

BsMac::~BsMac()  {
    /*for(int i = 0; i < numberOfMobileStations; ++i)  {
        ev << "BS Queue" << bsId << "to MS " << i << ": Queue " << packetQueue[i].length() << endl;
    }*/

    if(packetQueue != NULL)
    delete[] packetQueue;
    delete msPositions;
    delete downSchedule;
    delete transmitRequests;
    delete msPosUpdateArrived;
    delete neighbourIdMatching;
    
    // Destruction of scheduler triggers writing datarate information to file
    delete scheduler;
}
