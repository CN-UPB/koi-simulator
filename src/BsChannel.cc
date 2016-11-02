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

#include "BsChannel.h"
#include "ExpChannel.h"
#include "FactoryChannel.h"
#include "KoiData_m.h"
#include "METISChannel.h"
#include "SINR_m.h"
#include "PositionExchange_m.h"
#include "PointerExchange_m.h"
#include "Schedule_m.h"
#include "util.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>

Define_Module(BsChannel);

void BsChannel::initialize()  {
	maxNumberOfNeighbours = par("maxNumberOfNeighbours");
	upResBlocks = par("upResourceBlocks");
	downResBlocks = par("downResourceBlocks");
	numberOfMobileStations = par("numberOfMobileStations");
	tti = par("tti");
	epsilon = par("epsilon");
	initOffset = par("initOffset");
	bsId = par("bsId");

	//find the neighbours and store the pair (bsId, position in data structures) in a map
	cModule *cell = getParentModule()->getParentModule();
	neighbourIdMatching = new NeighbourIdMatching(bsId, maxNumberOfNeighbours, cell);

	// EESM Beta values for effective SINR
	string eesm_beta = par("eesm_beta");
	eesm_beta_values = vec(eesm_beta);

	// This counter counts how many neighbour have already transmitted all 
	// necessary position information.
	init_counter = 0;

	// Counts the received schedules from the own mac.
	scheduleCatch = false;

	int ch = par("channelModel");

	switch(ch){
		case 0:
			channel = new METISChannel();
			break;
		case 1:
			channel = new ExpChannel();
			break;
		case 2:
			channel = new FactoryChannel();
		default:
			throw std::invalid_argument("Invalid channelModel value");
	}

	// Save the CURRENT Direction of the schedule for next TTI for all Neighbours
	scheduleDirection = new int[neighbourIdMatching->numberOfNeighbours()];
	maxPower = new double[neighbourIdMatching->numberOfNeighbours()];
	schedulePower = new double*[neighbourIdMatching->numberOfNeighbours()];
	for(int i = 0; i < neighbourIdMatching->numberOfNeighbours(); i++){
		schedulePower[i] = new double[numberOfMobileStations];
	}

	// Prepare MS positions vector with correct sizes for each cell
	msPositions.resize(neighbourIdMatching->numberOfNeighbours(), 
			vector<Position>());
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
		PointerExchange *PtrMessage = new PointerExchange("POINTER_EXCHANGE2");
		PtrMessage->setPtr(channel);
		send(PtrMessage, "toPhy"); //set to 999*tti originally
	}
	scheduleAt(simTime()+initOffset-epsilon, new cMessage("SINR_ESTIMATION"));
}

void BsChannel::handleMessage(cMessage *msg)  {
	if(msg->isName("POINTER_EXCHANGE2"))  {
		//cout << "Received Channel Pointer!" << endl;
		PointerExchange *PtrMessage = (PointerExchange*) msg;
		Channel *old = channel;
		channel = PtrMessage->getPtr();
		if(old!=channel){
			// when this BS channel makes use of the agreed-upon
			// channel instance, we need to delete the instance 
			// created in it's own initialize method.
			delete old;
		}
		delete msg;
	}
	else if(msg->isName("BS_POSITION_MSG")) {
		PositionExchange *bsPos = (PositionExchange *) msg;
		neighbourPositions[bsPos->getId()] = bsPos->getPosition();
		init_counter++;
		if(this->getIndex()==0 && init_counter==2*maxNumberOfNeighbours){
			channel->init(this, msPositions, neighbourPositions);
		}
		delete msg;
	}
	else if(msg->isName("BS_MS_POSITIONS"))  {
		//save the postitions of the mobile stations
		BsMsPositions *msPos = (BsMsPositions *) msg;
		msPositions[msPos->getBsId()].resize(msPos->getPositionsArraySize());
		for(unsigned int i = 0; i < msPos->getPositionsArraySize(); i++)  {
			msPositions[msPos->getBsId()][i] = msPos->getPositions(i);
		}
		init_counter++;
		if(this->getIndex()==0 && init_counter==2*maxNumberOfNeighbours){
			channel->init(this, msPositions, neighbourPositions);
		}
		delete msPos;
	}
	else if(msg->isName("DEBUG")){
		if(this->getIndex()==0){
			// Only do debug output for the first BSChannel. The values 
			// are all the same for all BsChannel instances of any given 
			// Base Station.
			// Forward DEBUG message to the channel implementation
			channel->handleMessage(msg);
		}
	}
	else if(msg->isName("SINR_ESTIMATION")){
		SINR *bsSINREst = new SINR();
		bsSINREst->setBsId(bsId);
		// Special value to note that this is the SINR for 
		// a base station
		bsSINREst->setMsId(-1);
		// We only need the down SINR estimate, because the 
		// base station only ever uses DOWN resource blocks.
		bsSINREst->setDownArraySize(downResBlocks);
		for(int i = 0; i < downResBlocks; i++){
			bsSINREst->setDown(i,
					channel->calcAvgDownSINR(i,1.0));
		}
		// Route message to BS via MsPhy and MsMac
		send(bsSINREst,"toPhy");
		scheduleAt(simTime() + tti, msg);
	}
	else if(msg->getKind()==MessageType::transInfo){
		// We should only add the new transInfo to the channel for the BsChannel
		// with the index 0. Otherwise, as there is only one channel per cell but 
		// multiple BsChannels, we would add the same trans info message to the 
		// same channel instance multiple times. Clearly, a bad idea.
		if(this->getIndex()==0){
			TransInfo *inf = dynamic_cast<TransInfo*>(msg);
			channel->addTransInfo(inf);
		}
		else{
			// This is not the BsChannel which adds the trans info message, so 
			// delete it.
			delete msg;
		}
	}
	else if(msg->getKind()==MessageType::clearTransInfo){
		// Clear all TransInfo messages from the channel
		channel->clearTransInfo();
	}
	else if(msg->arrivedOn("fromMs"))  {
		assert(msg->getKind() == MessageType::koidata);
		KoiData *packet = dynamic_cast<KoiData*>(msg);
		// Set Scheduled flag to false again, the packet has been transmitted 
		// and now needs to be scheduled anew for the next transmission leg.
		packet->setScheduled(false);

		vector<double> instSINR;
		int currentRessourceBlock = packet->getResourceBlock();
		instSINR.push_back(channel->calcUpSINR(currentRessourceBlock,
					packet->getSrc(),
					packet->getTransPower()));

		/**
			double effSINR = getEffectiveSINR(instSINR,eesm_beta_values);
		//cout << "Effektive SINR (Up): " << effSINR << endl;
		double bler = getBler(packet->getCqi(), effSINR, this);
		//cout << "Block Error Rate(Up): " << bler << endl;
		vec bler_(1);
		bler_.set(0,bler);
		double per = getPer(bler_);
		 **/
		/**
			if(uniform(0,1) > per){
			sendDelayed(bundle, tti - epsilon, "toPhy");
			}else{
			delete bundle;
			}
		 **/
		//For now, all packets are send successfully
		sendDelayed(packet, tti - epsilon, "toPhy");

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
			delete[] schedulePower[i];
			delete[] schedules[i];
    }
    delete[] schedulePower;
    delete[] schedules;
    delete neighbourIdMatching;
    delete[] scheduleDirection;
    delete[] maxPower;
    if(this->getIndex() == 0){
	    delete channel;
    }
}
