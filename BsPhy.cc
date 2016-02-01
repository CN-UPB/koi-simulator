/*
 * BsPhy.cc
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#include "BsPhy.h"
#include "DataPacket_m.h"
#include "DataPacketBundle_m.h"
#include "PointerExchange_m.h"

Define_Module(BsPhy);

void BsPhy::initialize()  {
    numberOfMobileStations = par("numberOfMobileStations");
}

void BsPhy::handleMessage(cMessage *msg)  {
	if(msg->isName("POINTER_EXCHANGE2")){
		send(msg->dup(), "toMac");
        for(int i = 0; i < numberOfMobileStations; i++)  {
            send(msg->dup(), "toChannel", i);
        }
		delete msg;
	}
	else if(msg->isName("CHANNEL_INFO")){
		send(msg, "toMac");
	}
	else if(msg->isName("SINR")){
		std::cout << "rev sinr" << std::endl;
		delete msg;
	}
	else if(msg->isName("POINTER_EXCHANGE")){
        send(msg->dup(), "toMac");
        for(int i = 0; i < numberOfMobileStations; i++)  {
            send(msg->dup(), "toChannel", i);
        }
		delete msg;
	}
	else if(msg->isName("BS_POSITION_MSG"))  {
        for(int i = 0; i < numberOfMobileStations; i++)  {
            send(msg->dup(), "toChannel", i);
        }
        delete msg;
    }
    else if(msg->isName("TEST_VR"))  {
        send(msg, "toMac");
    }
    else if(msg->isName("VR_RETURN"))  {
        send(msg, "toMac");
    }
    else if(msg->isName("CLUSTER_INFO"))  {
        send(msg, "toMac");
    }
    else if(msg->isName("BS_CHANNEL_UPDATE"))  {
        for(int i = 0; i < numberOfMobileStations; i++)  {
            send(msg->dup(), "toChannel", i);
        }
        delete msg;
    }
    else if(msg->arrivedOn("fromMac"))  {
        assert(msg->isName("DATA_BUNDLE"));
        //forward the packet to the channel of the ms
        DataPacketBundle *bundle = (DataPacketBundle *) msg;
        //estimate the right data channel
        int msId = bundle->getMsId();
        //ev << "Forwarding packet to ms " << msId << "data channel" << endl;
        send(bundle, "toChannel", msId);
    }
    else if(msg->arrivedOn("fromChannel"))  {
        //forward the packet from the bs channel to the bs mac
        assert(msg->isName("DATA_BUNDLE"));
        //ev << "Forwarding packet bundle to BsMac" << endl;
        send(msg, "toMac");
    }
}
