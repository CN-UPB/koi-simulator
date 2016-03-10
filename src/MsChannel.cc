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

Define_Module(MsChannel);

void MsChannel::initialize()  {
    maxNumberOfNeighbours = par("maxNumberOfNeighbours");
    useSimpleChannelCalc = par("useSimpleChannelCalc");
    simpleChannelCalcNops = par("simpleChannelCalcNops");
    packetLoss = par("packetLoss");
    currentChannel = par("currentChannel");
    bsId = par("bsId");
    epsilon = par("epsilon");
    tti = par("tti");
    msId = par("msId");
    downResourceBlocks = par("downResourceBlocks");

    	msPosition.x = 6;
    	msPosition.y = 21;
    

    //find the neighbours and store the pair (bsId, position in data structures) in a map
    cModule *cell = getParentModule()->getParentModule();
    neighbourIdMatching = new NeighbourIdMatching(bsId, maxNumberOfNeighbours, cell);

    bsPositions = new Position[neighbourIdMatching->numberOfNeighbours()];
    bsChannels = new int[neighbourIdMatching->numberOfNeighbours()];
    
    // EESM Beta values for effective SINR
	string eesm_beta = par("eesm_beta");
	eesm_beta_values = vec(eesm_beta);
    
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
                
		//complex packet loss with fading and pathloss
		int ownDataStrPos = neighbourIdMatching->getDataStrId(bsId);
		
		pos.push_back(msPosition);
		pos.push_back(bsPositions[ownDataStrPos]);
		power.push_back(1.0);
		power.push_back(1.0);
		bsId_.push_back(bsId);
		bsId_.push_back(bsId);
		
		//interferer are all the neighbouring bs
		NeighbourMap *map = neighbourIdMatching->getNeighbourMap();
		for(NeighbourMap::iterator it = map->begin(); it != map->end(); ++it)  {
			if(it->first == bsId)
				continue; //skip the own bs; cant interfere; skip unused cell connections
				
			//only the bss that send on the same channel are interferer
			if(currentChannel == bsChannels[(it->second).first]){
				pos.push_back(bsPositions[(it->second).first]);
				power.push_back(1.0);
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
    else if(msg->isName("BS_CHANNEL_UPDATE"))  {
        ChannelExchange *chnEx = (ChannelExchange *) msg;
        int dataStrPos = neighbourIdMatching->getDataStrId(chnEx->getId());
        bsChannels[dataStrPos] = chnEx->getChannel();
        delete msg;
    }
	else if(msg->isName("BS_POSITION_MSG"))  {
        PositionExchange *bsPos = (PositionExchange *) msg;
        int dataStrPos = neighbourIdMatching->getDataStrId(bsPos->getId());
        bsPositions[dataStrPos] = bsPos->getPosition();
        delete msg;
    }
    else if(msg->arrivedOn("fromBs"))  {
        assert(msg->isName("DATA_BUNDLE") == true);

        //the channel receives the packet in a bundle
        DataPacketBundle *bundle = (DataPacketBundle *) msg;
	// Just forward the packet for now, without error checking etc
	sendDelayed(bundle, tti - epsilon, "toPhy");
		
		/*
        bool forwardPacket[numberOfPackets];

        for(unsigned int i = 0; i < bundle->getPacketsArraySize(); ++i)  {
            forwardPacket[i] = true;

            if(useSimpleChannelCalc)
                forwardPacket[i] = SimpleChannelCalc::calc(simpleChannelCalcNops, uniform(0,1), packetLoss);
            else  {
                //complex packet loss with fading and pathloss
                DataPacket packet = bundle->getPackets(i);
                int ownDataStrPos = neighbourIdMatching->getDataStrId(packet.getBsId());
                channelCalc.setSenderPosition(bsPositions[ownDataStrPos],1.0);
                channelCalc.setTargetPosition(msPosition);
                //interferer are all the neighbouring bs
                NeighbourMap *map = neighbourIdMatching->getNeighbourMap();
                for(NeighbourMap::iterator it = map->begin(); it != map->end(); ++it)  {
                    if(it->first == packet.getBsId())
                        continue; //skip the own bs; cant interfere; skip unused cell connections
                    //only the bss that send on the same channel are interferer
                    if(currentChannel == bsChannels[(it->second).first])
                        channelCalc.addInterfererPosition(bsPositions[(it->second).first],1.0);
                }
                forwardPacket[i] = channelCalc.isPacketOK(packet.getArrivalTime(), packet.getBitLength());
                channelCalc.clearInterfererPostitions();
            }

            if(forwardPacket[i])
                packetsToForward++; //for forward packet bundle
        }
        */

        vector<double> instSINR;
        int ownDataStrId = neighbourIdMatching->getDataStrId(bundle->getBsId());
        for(uint i = 0; i < bundle->getRBsArraySize(); i++){
			int currentRessourceBlock = bundle->getRBs(i);
		
			vector<double> power;
			vector<Position> pos;
			vector<int> bsId;
			
			bsId.push_back(bundle->getBsId());
			bsId.push_back(bundle->getBsId());
			pos.push_back(msPosition);
			pos.push_back(bsPositions[ownDataStrId]);
			power.push_back(1.0);
			power.push_back(1.0);

			NeighbourMap *map = neighbourIdMatching->getNeighbourMap();
			for(NeighbourMap::iterator it = map->begin(); it != map->end(); it++)  {
				if(it->first == bundle->getBsId())
					continue; //skip the own bs; cant interfere
				
                if(currentChannel == bsChannels[(it->second).first]){
					pos.push_back(bsPositions[(it->second).first]);
					power.push_back(1.0);
					bsId.push_back(it->first);
				}
			}

			instSINR.push_back(channel->calcSINR(currentRessourceBlock,power,pos,bsId,false, msId));
		}
		//cout << "-----------------------" << endl;
		/**
		double effSINR = getEffectiveSINR(instSINR,eesm_beta_values);
		double bler = getBler(bundle->getCqi(), effSINR, this);
		vec bler_(1);
		bler_.set(0,bler);
		double per = getPer(bler_);
		
		if(uniform(0,1) > per){
			sendDelayed(bundle, tti - epsilon, "toPhy");
		}else{
			delete bundle;
		}
		**/
	}
	
	/*
	if(packetsToForward > 0)  {
		//make a new bundle with the packets that are not lost
		DataPacketBundle *forwardBundle = new DataPacketBundle("DATA_BUNDLE");
		forwardBundle->setMsId(bundle->getMsId());
		forwardBundle->setBsId(bundle->getBsId());
		forwardBundle->setPacketsArraySize(packetsToForward);

		int curForwardBunldePos = 0;
		for(unsigned int i = 0; i < bundle->getPacketsArraySize(); ++i)  {
			if(forwardPacket[i])  { //forward the message to the ms
				//ev << "MsChannel sends packet " << i << " to the mac layer" << endl;
				forwardBundle->setPackets(curForwardBunldePos, bundle->getPackets(i));
			}
			else  {
				//ev << "Packet " << i << " is lost!" << endl;
			}
		}
		sendDelayed(forwardBundle, tti - epsilon, "toPhy");
		*/
	//}
}

MsChannel::~MsChannel()  {
    /*if(bsId == 0)  {
    ev << "Basestations----------------" << endl;
    NeighbourMap *map1 = neighbourIdMatching->getNeighbourMap();
    for(NeighbourMap::iterator i = map1->begin(); i != map1->end(); i++)  {
            ev << i->first << ": " << bsPositions[(i->second).first].x << ", " << bsPositions[(i->second).first].y << endl;
    }
    ev << "-------------------------------------------------------" << endl;
    }*/
    /*std::cout << "Basestations----------------" << std::endl;
    NeighbourMap *map1 = neighbourIdMatching->getNeighbourMap();
    for(NeighbourMap::iterator i = map1->begin(); i != map1->end(); i++)  {
        std::cout << i->first << ": " << bsChannels[(i->second).first] << std::endl;
    }
    std::cout << "-------------------------------------------------------" << std::endl;*/
    delete[] bsPositions;
    delete[] bsChannels;
    delete neighbourIdMatching;
}
