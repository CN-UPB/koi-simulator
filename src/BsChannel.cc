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
#include "KoiData_m.h"
#include "SINR_m.h"
#include "PositionExchange_m.h"
#include "PointerExchange_m.h"
#include "VisibilityRegionMessage_m.h"
#include "ClusterMessage_m.h"
#include "Schedule_m.h"
#include "util.h"
#include "METISChannel.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>

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

    // Prepare transmission info lists for every UP ressource block
    transInfos.resize(upResBlocks);

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
		PtrExchange Pointer;
		Pointer.ptr = (uintptr_t) channel;
		PointerExchange *PtrMessage = new PointerExchange("POINTER_EXCHANGE2");
		PtrMessage->setPtr(Pointer);
		send(PtrMessage, "toPhy"); //set to 999*tti originally
		
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
	else if(msg->isName("BS_POSITION_MSG")) {
		PositionExchange *bsPos = (PositionExchange *) msg;
		neighbourPositions[bsPos->getId()] = bsPos->getPosition();
		init_counter++;
		if(this->getIndex()==0 && init_counter==2*maxNumberOfNeighbours){
			channel->init(this, msPositions, neighbourPositions);
		}
		delete msg;
	}
    else if(msg->getKind()==MessageType::transInfo){
    	TransInfo *info = dynamic_cast<TransInfo*>(msg);
	transInfos[info->getRb()].push_front(info);
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
		// Write out SINR values for all possible links in this cell
		// First, the downlinks
		ofstream downSinr;
		std::string fname("sinr_table_down_"+std::to_string(bsId)+".dat");
		downSinr.open(fname,ofstream::trunc);
		outputDownSINR(downSinr);
		downSinr.close();
		// Second, the uplinks
		ofstream upSinr;
		fname = "sinr_table_up_"+std::to_string(bsId)+".dat";
		upSinr.open(fname,ofstream::trunc);
		outputUpSINR(upSinr);
		upSinr.close();
	}
    }
    else if(msg->arrivedOn("fromMs"))  {
        assert(msg->getKind() == MessageType::koidata);
        KoiData *packet = dynamic_cast<KoiData*>(msg);

        vector<double> instSINR;
	int currentRessourceBlock = packet->getResourceBlock();
	instSINR.push_back(channel->calcUpSINR(currentRessourceBlock,
              transInfos[currentRessourceBlock],
              packet->getSrc(),
              packet->getTransPower()));

	double effSINR = getEffectiveSINR(instSINR,eesm_beta_values);
	//cout << "Effektive SINR (Up): " << effSINR << endl;
	double bler = getBler(packet->getCqi(), effSINR, this);
	//cout << "Block Error Rate(Up): " << bler << endl;
	vec bler_(1);
	bler_.set(0,bler);
	double per = getPer(bler_);

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

std::ostream& BsChannel::outputDownSINR(std::ostream& out){
	double transPower = 1.0;
	int downResBlocks = par("downResourceBlocks");
	// Add header for downlink SINR table
	out << "MS\t" << "RB\t" << "SINR\t" << std::endl;
	// Generate transmission info messages for all neighbouring base 
	// stations. We assume that all neighbouring BS are sending and 
	// we generate for all downlink ressource blocks. We also assume 
	// that the transmission power used is always 1.
	forward_list<TransInfo> inf;
	for(int i=0; i<maxNumberOfNeighbours; i++){
		if(i!=bsId){
			TransInfo info;
			info.setRb(0);
			info.setPower(transPower);
			info.setBsId(i);
			info.setMessageDirection(MessageDirection::down);
			inf.push_front(info);
		}
	}
	for(int i=0; i<numberOfMobileStations; i++){
		forward_list<TransInfo*> currInf;
		for(auto c:inf){
			currInf.push_front(c.dup());
		}
		for(int r=0; r<downResBlocks; r++){
			out << i << "\t" 
				<< r << "\t" 
				<< channel->calcDownSINR(r,currInf,i,transPower)
				<< std::endl;
		}
	}
	return out;
}

std::ostream& BsChannel::outputUpSINR(std::ostream& out){
	double transPower = 1.0;
	// Add header for uplink SINR table
	out << "MS\t" << "RB\t" << "SINR\t" << std::endl;
	// Generate transmission info messages for all neighbouring mobile 
	// stations. We assume that all neighbouring MS are sending and 
	// we generate for all uplink ressource blocks. We also assume 
	// that the transmission power used is always 1.
	forward_list<TransInfo> inf;
	for(int j=0; j<maxNumberOfNeighbours; j++){
		if(j!=bsId){
			int numMs = neighbourIdMatching->getNumberOfMS(j);
			for(int i=0; i<numMs; i++){
				TransInfo info;
				info.setRb(0);
				info.setPower(transPower);
				info.setBsId(j);
				info.setMsId(i);
				info.setMessageDirection(MessageDirection::up);
				inf.push_front(info);
			}
		}
	}
	for(int i=0; i<numberOfMobileStations; i++){
		forward_list<TransInfo*> currInf;
		for(auto c:inf){
			currInf.push_front(c.dup());
		}
		for(int r=0; r<upResBlocks; r++){
			out << i << "\t" 
				<< r << "\t" 
				<< channel->calcUpSINR(r,currInf,i,transPower)
				<< std::endl;
		}
	}
	return out;
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
