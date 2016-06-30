/*
 * MsChannel.cc
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerlin
 */

#include "MsChannel.h"
#include "DataPacket_m.h"
#include "DataPacketBundle_m.h"
#include "SINR_m.h"
#include "PositionExchange_m.h"
#include "ChannelExchange_m.h"
#include "SimpleChannelCalc.h"
#include "PointerExchange_m.h"
#include "util.h"
#include <math.h>
#include <cmath>
#include <itpp/itbase.h>

using std::vector;
using std::forward_list;

Define_Module(MsChannel);

void MsChannel::initialize()  {
    maxNumberOfNeighbours = par("maxNumberOfNeighbours");
    packetLoss = par("packetLoss");
    bsId = par("bsId");
    epsilon = par("epsilon");
    tti = par("tti");
    msId = par("msId");
    downResourceBlocks = par("downResourceBlocks");
    upResourceBlocks = par("upResourceBlocks");
    msPosition.x = 6;
    msPosition.y = 21;
    //find the neighbours and store the pair (bsId, position in data structures) in a map
    cModule *cell = getParentModule()->getParentModule();
    neighbourIdMatching = new NeighbourIdMatching(bsId, maxNumberOfNeighbours, cell);

    bsPositions = new Position[neighbourIdMatching->numberOfNeighbours()];
    
    // EESM Beta values for effective SINR
    string eesm_beta = par("eesm_beta");
    eesm_beta_values = vec(eesm_beta);

    // Instantiate a transmission info list for each down/up ressource block
    transInfosDown.resize(downResourceBlocks);
    transInfosUp.resize(upResourceBlocks);
    
    scheduleAt(simTime() + 1000 * tti + epsilon, new cMessage("SINR_ESTIMATION")); //originally set to 1000*tti + epsilon
}

simtime_t MsChannel::getProcessingDelay(cMessage *msg)  {
    if(msg->isName("DATA_BUNDLE"))
        return tti - 2 * epsilon;
    else
        return 0;
}

void MsChannel::handleMessage(cMessage *msg)  {
	if(msg->isName("SINR_ESTIMATION"))  {
		vector<double> power;
		vector<Position> pos;
		vector<int> bsId_;

		bsId_.push_back(bsId);
		bsId_.push_back(bsId);

		//interferer are all the neighbouring bs
		NeighbourMap *map = neighbourIdMatching->getNeighbourMap();
		for(NeighbourMap::iterator it = map->begin(); it != map->end(); ++it)  {
			if(it->first != bsId){
				bsId_.push_back(it->first);
			}


		}
		if(!channel){
			std::cout << "ERROR NULL POINTER!" << std::endl;
		}
		vec sinr = channel->calcSINR(power, pos, bsId_, false, msId);
		//cout << "SINR: " << sinr << endl;
		//cout << "Nr: " << power.size() << " " << pos.size() << " " << bsId_.size() << endl;

		delete msg;
		scheduleAt(simTime() + tti, new cMessage("SINR_ESTIMATION"));

		SINR *sinrMessage = new SINR("SINR_ESTIMATION");
		sinrMessage->setSINRArraySize(downResourceBlocks);
		//std::cout << "SINR at RB 0 = " << sinr(0) << std::endl;
		for(int i = 0; i < downResourceBlocks; i++){
			sinrMessage->setSINR(i,sinr(i));
		}
		//std::cout << "Send SINR Estimation from MS " << msId << std::endl;
		//std::cout << sinr << std::endl;
		send(sinrMessage,"toPhy");
	}
	else if(msg->isName("CHANNEL_INFO"))  {
		// Whenever you receive a message called CHANNEL_INFO forward it to channel.
		// This way channel can realise arbitrary communication if necessary.
		channel->handleMessage(msg);
	}
	else if(msg->isName("POINTER_EXCHANGE2"))  {
		//cout << "channel arrived at MS" << endl;
		PointerExchange *PtrMessage = (PointerExchange*) msg;
		channel = (Channel*) (PtrMessage->getPtr()).ptr;
		delete msg;
	}
	else if(msg->isName("MS_POS_UPDATE"))  {
		//save the position update of the mobile stations
		PositionExchange *posEx = (PositionExchange *) msg;
		msPosition = posEx->getPosition();
		delete msg;
	}
	else if(msg->isName("BS_POSITION_MSG"))  {
		PositionExchange *bsPos = (PositionExchange *) msg;
		int dataStrPos = neighbourIdMatching->getDataStrId(bsPos->getId());
		bsPositions[dataStrPos] = bsPos->getPosition();
		delete msg;
	}
	else if(msg->getKind()==MessageType::transInfo){
		TransInfo *info = dynamic_cast<TransInfo*>(msg);
		// The transInfo lists are sorted by RB by transmission direction
		switch(info->getMessageDirection()){
			case MessageDirection::d2dUp:
			case MessageDirection::up:
				transInfosUp[info->getRb()].push_front(info);
				break;
			case MessageDirection::d2dDown:
			case MessageDirection::down:
				transInfosDown[info->getRb()].push_front(info);
				break;

		}
	}
	else if(msg->isName("DATA_BUNDLE"))  {
		//the channel receives the packet in a bundle
		DataPacketBundle *bundle = (DataPacketBundle *) msg;
		// Just forward the packet for now, without error checking etc

		vector<double> instSINR;
		int currentRessourceBlock = bundle->getRBs(0);

		switch(bundle->getMessageDirection()){
			case MessageDirection::down:
				instSINR.push_back(channel->calcDownSINR(currentRessourceBlock,transInfosDown[currentRessourceBlock],msId,bundle->getTransPower()));
				break;
			case MessageDirection::d2dDown:
				instSINR.push_back(channel->calcD2DSINR(
							currentRessourceBlock,
							transInfosDown[currentRessourceBlock],
							bundle->getMsId(),
							msId,MessageDirection::d2dDown,
							bundle->getTransPower()));
				break;
			case MessageDirection::d2dUp:
				instSINR.push_back(channel->calcD2DSINR(
							currentRessourceBlock,
							transInfosUp[currentRessourceBlock],
							bundle->getMsId(),
							msId,MessageDirection::d2dUp,
							bundle->getTransPower()));
				break;
		}
		double effSINR = getEffectiveSINR(instSINR,eesm_beta_values);
		double bler = getBler(bundle->getCqi(), effSINR, this);
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
		// For now, all packets are received successfully
		sendDelayed(bundle, tti - epsilon, "toPhy");
	}
	
}

MsChannel::~MsChannel()  {
    delete[] bsPositions;
    delete neighbourIdMatching;
}
