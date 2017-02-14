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

using namespace omnetpp;
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
	debug = par("debug");
	d2dActive = par("d2dActive");
	downResourceBlocks = par("downResourceBlocks");
	upResourceBlocks = par("upResourceBlocks");
	numMSAntenna = par("NumMsAntenna");
	numBSAntenna = par("NumBsAntenna");
	coding.init(par("MCSTable"),tti.dbl(),par("bandwidthPerRB"));

	// EESM Beta values for effective SINR
	string eesm_beta = par("eesm_beta");
	eesm_beta_values = vec(eesm_beta);
	if(debug){
		// File for SINR value storage
		std::string fname("sinr-ms-"+std::to_string(bsId)+"-"+std::to_string(msId));
		sinrFile = std::move(getResultFile(fname));
		sinrFile << "TTI\t" 
			<< "Cell\t" << "MS\t" << "RB\t" << "SINR" << std::endl;
	}
	
	scheduleAt(simTime()+initOffset-epsilon, new cMessage("SINR_ESTIMATION")); //originally set to 1000*tti + epsilon
}

void MsChannel::finish(){
	if(debug){
		sinrFile.close();
	}
}

void MsChannel::handleMessage(cMessage *msg)  {
	if(msg->isName("SINR_ESTIMATION")){
		SINR *sinrMessage = new SINR();
		sinrMessage->setBsId(bsId);
		sinrMessage->setMsId(msId);

		// Set SINR and rate estimation to the average SINR value over all 
		// possible transmission recipients in the previous tti.
		double sinr = 0.0;
		if(d2dActive){
			// Values for DOWN resource blocks are only needed with D2D, otherwise
			// the MS don't use the DOWN RBs.
			sinrMessage->setDownArraySize(downResourceBlocks);
			sinrMessage->setRDownArraySize(downResourceBlocks);
			for(int i = 0; i < downResourceBlocks; i++){
				sinr = channel->calcAvgD2DDownSINR(i,msId,1.0);
				sinrMessage->setDown(i,sinr);
				sinrMessage->setRDown(i,coding.getRBCapacity(sinr,numMSAntenna,
							numBSAntenna));
			}
		}
		sinrMessage->setUpArraySize(upResourceBlocks);
		sinrMessage->setRUpArraySize(upResourceBlocks);
		for(int i = 0; i < upResourceBlocks; i++){
			sinr = channel->calcUpSINR(i,msId,1.0);
			sinrMessage->setUp(i,sinr);
			sinrMessage->setRUp(i,coding.getRBCapacity(sinr,numMSAntenna,
						numBSAntenna));
		}
		// Route estimate to MsMac via MsPhy
		send(sinrMessage,"toPhy");
		scheduleAt(simTime() + tti, msg);
		if(simTime()>initOffset && debug){
			auto val = (simTime()-initOffset)/tti;
			int tti = std::floor(val);
			for(int i = 0; i < upResourceBlocks; i++){
				sinrFile << tti << "\t" << bsId << "\t" << msId << "\t" << i << "\t"
					<< sinrMessage->getUp(i) << std::endl;
			}
		}
	}
	else if(msg->isName("POINTER_EXCHANGE2"))  {
		//cout << "channel arrived at MS" << endl;
		PointerExchange *PtrMessage = (PointerExchange*) msg;
		channel = PtrMessage->getPtr();
		delete msg;
	}
	else if(msg->getKind()==MessageType::koidata)  {
		KoiData *packet = (KoiData *) msg;
                // Set Scheduled to false, as the packet now need to be
                // scheduled anew.
                packet->setScheduled(false);
		// Just forward the packet for now, without error checking etc
		/**
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
		double effSINR = getEffectiveSINR(instSINR,eesm_beta_values);
		double bler = getBler(packet->getCqi(), effSINR, this);
		vec bler_(1);
		bler_.set(0,bler);
		double per = getPer(bler_);
		if(uniform(0,1) > per){
			sendDelayed(bundle, tti - epsilon, "toPhy");
		}else{
			delete bundle;
		}
		**/
		// For now, all packets are received successfully
		sendDelayed(packet, epsilon, "toPhy");
	}
	else if(msg->getKind()==MessageType::longTermSinrEst){
		SINR* longtermEst = dynamic_cast<SINR*>(msg);
		longtermEst->setBsId(bsId);
		longtermEst->setMsId(msId);

		// Set SINR and rate estimation to the average SINR value over all 
		// possible transmission recipients in the previous tti.
		double sinr = 0.0;
		if(d2dActive){
			// Values for DOWN resource blocks are only needed with D2D, otherwise
			// the MS don't use the DOWN RBs.
			longtermEst->setDownArraySize(downResourceBlocks);
			longtermEst->setRDownArraySize(downResourceBlocks);
			for(int i = 0; i < downResourceBlocks; i++){
				sinr = channel->calcLongtermDownSINR(i,msId,1.0);
				longtermEst->setDown(i,sinr);
				longtermEst->setRDown(i,coding.getRBCapacity(sinr,numMSAntenna,
							numBSAntenna));
			}
		}
		longtermEst->setUpArraySize(upResourceBlocks);
		longtermEst->setRUpArraySize(upResourceBlocks);
		for(int i = 0; i < upResourceBlocks; i++){
			sinr = channel->calcLongtermUpSINR(i,msId,1.0);
			longtermEst->setUp(i,sinr);
			longtermEst->setRUp(i,coding.getRBCapacity(sinr,numMSAntenna,
						numBSAntenna));
		}
		// Route estimate to MsMac via MsPhy
		send(longtermEst,"toPhy");
	}
}

MsChannel::~MsChannel()  {
}
