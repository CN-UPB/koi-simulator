/*
 * BsChannel.cc
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 * Last edited: Jun 25, 2014
 *      Author: Thomas Prinz
 * 
 * All book references target the book "Pervasive Mobile and Ambient Wireless Communications"
 */

#include "cluster.h"
#include "BsChannel.h"
#include "DataPacket_m.h"
#include "DataPacketBundle_m.h"
#include "SINR_m.h"
#include "PositionExchange_m.h"
#include "PointerExchange_m.h"
#include "VisibilityRegionMessage_m.h"
#include "ClusterMessage_m.h"
#include "Schedule_m.h"
#include "SimpleChannelCalc.h"
#include "Cost2100Channel.h"
#include "ChannelAlternative.h"
#include "util.h"
#include "METISChannel.h"
#include <cmath>
#include <algorithm>

Define_Module(BsChannel);

void BsChannel::initialize()  {
    maxNumberOfNeighbours = par("maxNumberOfNeighbours");
    upResBlocks = par("upResourceBlocks");
    numberOfMobileStations = par("numberOfMobileStations");
    useSimpleChannelCalc = par("useSimpleChannelCalc");
    simpleChannelCalcNops = par("simpleChannelCalcNops");
    packetLoss = par("packetLoss");
    tti = par("tti");
    epsilon = par("epsilon");
    bsId = par("bsId");
    
	//find the neighbours and store the pair (bsId, position in data structures) in a map
    cModule *cell = getParentModule()->getParentModule();
    neighbourIdMatching = new NeighbourIdMatching(bsId, maxNumberOfNeighbours, cell);
    
    // EESM Beta values for effective SINR
	string eesm_beta = par("eesm_beta");
	eesm_beta_values = vec(eesm_beta);
    
    // This counter counts how many neighbour have already transmitted all necessary VR information.
    // Initialized with one because BS knows its own VR.
    init_counter = 1;
    
    // Counts the received schedules from the own mac.
    scheduleCatch = false;
    
    int ch = par("channelModel");
    
    if(ch == 0){
		channel = new METISChannel();
	}else if(ch == 1){
		channel = new Cost2100Channel();
	}else{
		channel = new ChannelAlternative();
	}
      
    // Save the CURRENT Direction of the schedule for next TTI for all Neighbours
    scheduleDirection = new int[neighbourIdMatching->numberOfNeighbours()];
    maxPower = new double[neighbourIdMatching->numberOfNeighbours()];
    schedulePower = new double*[neighbourIdMatching->numberOfNeighbours()];
    for(int i = 0; i < neighbourIdMatching->numberOfNeighbours(); i++){
		schedulePower[i] = new double[numberOfMobileStations];
	}

    //stores the pos of all ms from all bs
    msPositions = new Position*[neighbourIdMatching->numberOfNeighbours()];
    if(bsId == 3){
		std::cout << "Num of Neighbours: " << neighbourIdMatching->numberOfNeighbours() << std::endl;
	}
    for(int i = 0; i < neighbourIdMatching->numberOfNeighbours(); i++)  {
        msPositions[i] = new Position[numberOfMobileStations]; //TODO change to dynamic size per bs; when uneven partition is allowed
    }

    //the position of the base stations
    bsPosition.x = par("xPos");
    bsPosition.y = par("yPos");

    //stores the schedule of all neighbours
    schedules = new int*[neighbourIdMatching->numberOfNeighbours()];
    for(int i = 0; i < neighbourIdMatching->numberOfNeighbours(); ++i)  {
        schedules[i] = new int[upResBlocks];
        for(int j = 0; j < upResBlocks; j++)
            schedules[i][j] = -1;
    }
    if(this->getIndex() == 0){
		scheduleAt(simTime() + 750*tti - epsilon, new cMessage("TEST_MATRIX")); //set to 750*tti
		
		PtrExchange Pointer;
		Pointer.ptr = (uintptr_t) channel;
		PointerExchange *PtrMessage = new PointerExchange("POINTER_EXCHANGE2");
		PtrMessage->setPtr(Pointer);
		sendDelayed(PtrMessage, 999*tti , "toPhy"); //set to 999*tti originally
		
	}
}

simtime_t BsChannel::getProcessingDelay(cMessage *msg)  {
    if(msg->isName("DATA_BUNDLE"))
        return tti - 2 * epsilon;
    else
        return 0;
}

void BsChannel::handleMessage(cMessage *msg)  {
	if(msg->isName("CHANNEL_INFO"))  {
		channel->handleMessage(msg);
	}
	else if(msg->isName("POINTER_EXCHANGE2"))  {
		//cout << "Received Channel Pointer!" << endl;
		PointerExchange *PtrMessage = (PointerExchange*) msg;
		Channel *old = channel;
		channel = (Channel*) (PtrMessage->getPtr()).ptr;
		if(old!=channel){
			// when this BS channel makes use of the agreed-upon
			// channel instance, we need to delete the instance 
			// created in it's own initialize method.
			delete old;
		}
		delete msg;
	}
	else if(msg->isName("TEST_MATRIX")){
		channel->init(this, msPositions, neighbourPositions);
		//scheduleAt(simTime() + tti, msg);
		delete msg;
    }
	else if(msg->isName("BS_POSITION_MSG")) {
        PositionExchange *bsPos = (PositionExchange *) msg;
        neighbourPositions[bsPos->getId()] = bsPos->getPosition();
        delete msg;
    }
    else if(msg->isName("BS_MS_POSITIONS"))  {
        //save the postitions of the mobile stations
        BsMsPositions *msPos = (BsMsPositions *) msg;
        int fromBsId = neighbourIdMatching->getDataStrId(msPos->getBsId());
        //cout << "Ms positions from bs " << fromBsId << " arrived!" << endl;
        for(unsigned int i = 0; i < msPos->getPositionsArraySize(); i++)  {
            msPositions[fromBsId][i] = msPos->getPositions(i);
        }
        if(bsId == 3 && fromBsId == 3){
			string out2 = "Positions.txt";
			ofstream output2(out2);
			for(int i = 0; i < numberOfMobileStations; i++){
				//output2 << msPositions[fromBsId][i].x << " " << msPositions[fromBsId][i].y
				 //<< " " <<  1/pow(10, (22*log10( sqrt( pow(msPositions[fromBsId][i].x - 500.0,2) + pow(msPositions[fromBsId][i].y - 500.0,2) ) ) + 28 + 20*log10(2.5)) /10) << " \t " << 1/pow(10, (22*log10( sqrt( pow(msPositions[fromBsId][i].x - 500.0,2) + pow(msPositions[fromBsId][i].y - 1000.0,2) ) ) + 28 + 20*log10(2.5)) /10) << " \t " << 1/pow(10, (22*log10( sqrt( pow(msPositions[fromBsId][i].x - 500.0,2) + pow(msPositions[fromBsId][i].y - 0.0,2) ) ) + 28 + 20*log10(2.5)) /10) << " \t " << 1/pow(10, (22*log10( sqrt( pow(msPositions[fromBsId][i].x - 933.0,2) + pow(msPositions[fromBsId][i].y - 750.0,2) ) ) + 28 + 20*log10(2.5)) /10) << " \t " << 1/pow(10, (22*log10( sqrt( pow(msPositions[fromBsId][i].x - 67.0,2) + pow(msPositions[fromBsId][i].y - 750.0,2) ) ) + 28 + 20*log10(2.5)) /10) << " \t " << 1/pow(10, (22*log10( sqrt( pow(msPositions[fromBsId][i].x - 933.0,2) + pow(msPositions[fromBsId][i].y - 250.0,2) ) ) + 28 + 20*log10(2.5)) /10) << " \t " << 1/pow(10, (22*log10( sqrt( pow(msPositions[fromBsId][i].x - 67.0,2) + pow(msPositions[fromBsId][i].y - 250.0,2) ) ) + 28 + 20*log10(2.5)) /10) << endl;
			}
			output2.close();
		}
        if(simTime() >= 1 && this->getIndex()==0){
		// The channel instance is shared among all BsChannel instances,
		// thus only one BsChannel actually needs to call the update 
		// method.
		//channel->updateChannel(msPositions); 
	}
        delete msPos;
    }
    else if(msg->arrivedOn("fromMs"))  {
		//std::cout << "received fromMs Message." << std::endl;
        assert(msg->isName("DATA_BUNDLE") == true);

        //the channel receives the packet in a bundle
        DataPacketBundle *bundle = (DataPacketBundle *) msg;
	//sendDelayed(bundle, tti - epsilon, "toPhy");

        vector<double> instSINR;
        
	int currentRessourceBlock = bundle->getRBs(0);
	vector<double> power;
	vector<Position> pos;
	vector<int> bsId;

	bsId.push_back(bundle->getBsId());
	bsId.push_back(bundle->getBsId());

	NeighbourMap *map = neighbourIdMatching->getNeighbourMap();
	for(NeighbourMap::iterator it = map->begin(); it != map->end(); it++)  {
		if(it->first == bundle->getBsId())
			continue; //skip the own bs; cant interfere
		int interfererId = schedules[(it->second).first][currentRessourceBlock];

		if(interfererId != -1)  { //skip unused slots
			bsId.push_back(it->first);

			//cout << "BS " << it->first << ", MS " << interfererId << " is an interferer for the packet!"  << endl;
		}
	}
	instSINR.push_back(channel->calcSINR(currentRessourceBlock,power,pos,bsId,true,bundle->getMsId()));

	double effSINR = getEffectiveSINR(instSINR,eesm_beta_values);
	//cout << "Effektive SINR (Up): " << effSINR << endl;
	double bler = getBler(bundle->getCqi(), effSINR, this);
	//cout << "Block Error Rate(Up): " << bler << endl;
	vec bler_(1);
	bler_.set(0,bler);
	double per = getPer(bler_);

	if(uniform(0,1) > per){
		sendDelayed(bundle, tti - epsilon, "toPhy");
	}else{
		delete bundle;
	}

    }
}

BsChannel::~BsChannel()  {
    /*ev << "Bs Channel---------------------------------------------" << endl;
    NeighbourMap *map1 = neighbourIdMatching->getNeighbourMap();
    for(NeighbourMap::iterator i = map1->begin(); i != map1->end(); i++)  {
        for(int j = 0; j < numberOfMobileStations; ++j)  {
            ev << "BS " << i->first << ", MS " << j << ": " << msPositions[(i->second).first][j].x << ", " << msPositions[(i->second).first][j].y << endl;
        }
    }
    ev << "-------------------------------------------------------" << endl;*/

    /*ev << "Bs Channel Schedules---------------------------------------" << endl;
    NeighbourMap *map2 = neighbourIdMatching->getNeighbourMap();
     for(NeighbourMap::iterator i = map2->begin(); i != map2->end(); i++)  {
        for(int j = 0; j < upResBlocks; ++j)  {
            ev << "BS " << i->first << ", Slot " << j << ": " << schedules[(i->second).first][j] << endl;
        }
    }
    ev << "-------------------------------------------------------" << endl;*/

    for(int i = 0; i < neighbourIdMatching->numberOfNeighbours(); ++i)  {
        delete[] msPositions[i];
	delete[] schedulePower[i];
        delete[] schedules[i];
    }
    delete[] msPositions;
    delete[] schedulePower;
    delete[] schedules;
    delete neighbourIdMatching;
    delete[] scheduleDirection;
    delete[] maxPower;
    if(this->getIndex() == 0){
	    delete channel;
    }
}
