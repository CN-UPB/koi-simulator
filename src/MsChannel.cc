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
    bsId = par("bsId");
    epsilon = par("epsilon");
    initOffset = par("initOffset");
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
    
    scheduleAt(simTime()+initOffset-epsilon, new cMessage("SINR_ESTIMATION")); //originally set to 1000*tti + epsilon
}

simtime_t MsChannel::getProcessingDelay(cMessage *msg)  {
    if(msg->isName("DATA_BUNDLE"))
        return tti - 2 * epsilon;
    else
        return 0;
}

void MsChannel::handleMessage(cMessage *msg)  {
	if(msg->isName("SINR_ESTIMATION")){
		SINR *sinrMessage = new SINR();
                sinrMessage->setBsId(bsId);
                sinrMessage->setMsId(msId);

                // Set SINR estimation to the average SINR value over all 
                // possible transmission recipients in the previous tti.
		sinrMessage->setDownArraySize(downResourceBlocks);
		for(int i = 0; i < downResourceBlocks; i++){
			sinrMessage->setDown(i,
                            channel->calcAvgD2DDownSINR(i,transInfosDown[i],msId,1.0));
		}
		sinrMessage->setUpArraySize(upResourceBlocks);
		for(int i = 0; i < upResourceBlocks; i++){
			sinrMessage->setUp(i,
                            channel->calcAvgUpSINR(i,transInfosUp[i],msId,1.0));
		}
                if(msId==0){
                  // We need to compute the SINR estimate for the local base 
                  // station too, because the BS does not have all the 
                  // necessary transmission information available. 
                  // To prevent doing unncessary work, only each cell's 
                  // mobile station with the ID 0 computes the BS's 
                  // SINR estimate.
                  //
                  // TODO Think of a better way to do this
                  
                  SINR *bsSINREst = new SINR();
                  bsSINREst->setBsId(bsId);
                  // Special value to note that this is the SINR for 
                  // a base station
                  bsSINREst->setMsId(-1);
                  // We only need the down SINR estimate, because the 
                  // base station only ever uses DOWN resource blocks.
                  bsSINREst->setDownArraySize(downResourceBlocks);
                  for(int i = 0; i < downResourceBlocks; i++){
                    bsSINREst->setDown(i,
                        channel->calcAvgDownSINR(i,transInfosDown[i],1.0));
                  }
                  // Route message to BS via MsPhy and MsMac
                  send(bsSINREst,"toPhy");
                }
                // Route extimate to MsMac via MsPhy
		send(sinrMessage,"toPhy");
		scheduleAt(simTime() + tti, msg);
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
