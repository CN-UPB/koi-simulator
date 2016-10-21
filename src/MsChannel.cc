/*
 * MsChannel.cc
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerlin
 */

#include "KoiData_m.h"
#include "MsChannel.h"
#include "SINR_m.h"
#include "PositionExchange_m.h"
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

	// File for SINR value storage
	std::string fname("sinr-ms-"+std::to_string(bsId)+"-"+std::to_string(msId));
	sinrFile = std::move(getResultFile(fname));
	sinrFile << "TTI\t" 
		<< "Cell\t" << "MS\t" << "RB\t" << "SINR" << std::endl;
	
	scheduleAt(simTime()+initOffset-epsilon, new cMessage("SINR_ESTIMATION")); //originally set to 1000*tti + epsilon
}

void MsChannel::finish(){
	sinrFile.close();
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
                            channel->calcAvgD2DDownSINR(i,msId,1.0));
		}
		sinrMessage->setUpArraySize(upResourceBlocks);
		for(int i = 0; i < upResourceBlocks; i++){
			sinrMessage->setUp(i,
                            channel->calcUpSINR(i,msId,1.0));
		}
		// Route estimate to MsMac via MsPhy
		send(sinrMessage,"toPhy");
		scheduleAt(simTime() + tti, msg);
		if(simTime()>initOffset){
			auto val = (simTime()-initOffset)/tti;
			int tti = std::floor(val);
			for(int i = 0; i < upResourceBlocks; i++){
				sinrFile << tti << "\t" << bsId << "\t" << msId << "\t" << i << "\t"
					<< channel->calcUpSINR(i,msId,1.0) << std::endl;
			}
		}
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
	else if(msg->getKind()==MessageType::koidata)  {
		KoiData *packet = (KoiData *) msg;
                // Set Scheduled to false, as the packet now need to be
                // scheduled anew.
                packet->setScheduled(false);
		// Just forward the packet for now, without error checking etc
		vector<double> instSINR;
		int currentRessourceBlock = packet->getResourceBlock();

		switch(packet->getMessageDirection()){
			case MessageDirection::down:
				instSINR.push_back(channel->calcDownSINR(currentRessourceBlock,msId,packet->getTransPower()));
				break;
			case MessageDirection::d2dDown:
				instSINR.push_back(channel->calcD2DSINR(
							currentRessourceBlock,
							packet->getSrc(),
							msId,MessageDirection::d2dDown,
							packet->getTransPower()));
				break;
			case MessageDirection::d2dUp:
				instSINR.push_back(channel->calcD2DSINR(
							currentRessourceBlock,
							packet->getSrc(),
							msId,MessageDirection::d2dUp,
							packet->getTransPower()));
				break;
		}
                /**
		double effSINR = getEffectiveSINR(instSINR,eesm_beta_values);
		double bler = getBler(packet->getCqi(), effSINR, this);
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
		// For now, all packets are received successfully
		sendDelayed(packet, tti - epsilon, "toPhy");
	}
	
}

MsChannel::~MsChannel()  {
    delete[] bsPositions;
    delete neighbourIdMatching;
}
