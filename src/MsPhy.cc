/*
 * MsPhy.cc
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#include "MsPhy.h"
#include "SINR_m.h"
#include "MessageTypes.h"

Define_Module(MsPhy);

void MsPhy::initialize()  {
}

void MsPhy::handleMessage(cMessage *msg)  {
    if(msg->isName("SINR_ESTIMATION"))  {
        send(msg, "toMac");
    }
    else if(msg->isName("MS_POS_UPDATE"))  {
        send(msg, "toMsChannel");
    }
    //currently it only forward the packets
    else if(msg->arrivedOn("fromMac"))  {
	if(msg->getKind()==MessageType::transInfo){
		send(msg,"toMsChannel");
	}
	else{
		send(msg, "toChannel");
	}
    }
    else if(msg->arrivedOn("fromChannel"))  {
        //ev << "Forwarding packet/packets to the mac layer" << endl;
        send(msg, "toMac");
    }
}
