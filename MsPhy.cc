/*
 * MsPhy.cc
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#include "MsPhy.h"
#include "SINR_m.h"

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
        //ev << "Forwarding packet/packets to the data channel" << endl;
        send(msg, "toChannel");
    }
    else if(msg->arrivedOn("fromChannel"))  {
        //ev << "Forwarding packet/packets to the mac layer" << endl;
        send(msg, "toMac");
    }
}
