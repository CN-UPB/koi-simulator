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
#include "ResultFileExchange_m.h"
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
		sinrFile = getResultFile(fname);
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
			sinrMessage->setMcsDownArraySize(downResourceBlocks);
			for(int i = 0; i < downResourceBlocks; i++){
				sinr = channel->calcAvgD2DDownSINR(i,msId,1.0);
				int mcs;
				unsigned cap;
				std::tie(mcs,cap) = coding.getRBCapacity(sinr,numMSAntenna,
						numMSAntenna);
				sinrMessage->setRDown(i,cap);
				sinrMessage->setMcsDown(i,mcs);
				sinrMessage->setDown(i,sinr);
			}
		}
		sinrMessage->setUpArraySize(upResourceBlocks);
		sinrMessage->setRUpArraySize(upResourceBlocks);
		sinrMessage->setMcsUpArraySize(upResourceBlocks);
		for(int i = 0; i < upResourceBlocks; i++){
			sinr = channel->calcUpSINR(i,msId,1.0);
			int mcs;
			unsigned cap;
			std::tie(mcs,cap) = coding.getRBCapacity(sinr,numMSAntenna,
					numBSAntenna);
			sinrMessage->setUp(i,sinr);
			sinrMessage->setRUp(i,cap);
			sinrMessage->setMcsUp(i,mcs);
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
	else if(msg->isName("MCS_FILE")){
		ResultFileExchange* mcs = dynamic_cast<ResultFileExchange*>(msg);
		mcs_file = mcs->getPtr();
		delete mcs;
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
		// Log the MCS this Packet ought to have used and the best one it could
		// use based on the actual SINR value.
		double sinr = channel->calcDownSINR(packet->getResourceBlock(),
				msId,
				packet->getTransPower());
		unsigned cap;
		int mcs;
		std::tie(mcs,cap) = coding.getRBCapacity(sinr,numBSAntenna,numMSAntenna);
		*mcs_file << -1 << "\t" << packet->getMcs() << "\t"
			<< mcs << std::endl;
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
			longtermEst->setMcsDownArraySize(downResourceBlocks);
			for(int i = 0; i < downResourceBlocks; i++){
				sinr = channel->calcLongtermDownSINR(i,msId,1.0);
				int mcs;
				unsigned cap;
				std::tie(mcs,cap) = coding.getRBCapacity(sinr,numMSAntenna,
						numMSAntenna);
				longtermEst->setRDown(i,cap);
				longtermEst->setMcsDown(i,mcs);
				longtermEst->setDown(i,sinr);
			}
		}
		longtermEst->setUpArraySize(upResourceBlocks);
		longtermEst->setRUpArraySize(upResourceBlocks);
		longtermEst->setMcsUpArraySize(upResourceBlocks);
		for(int i = 0; i < upResourceBlocks; i++){
			sinr = channel->calcLongtermUpSINR(i,msId,1.0);
			int mcs;
			unsigned cap;
			std::tie(mcs,cap) = coding.getRBCapacity(sinr,numMSAntenna,
					numBSAntenna);
			longtermEst->setUp(i,sinr);
			longtermEst->setRUp(i,cap);
			longtermEst->setMcsUp(i,mcs);
		}
		// Route estimate to MsMac via MsPhy
		send(longtermEst,"toPhy");
	}
}

MsChannel::~MsChannel()  {
}
