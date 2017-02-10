/*
 * BsPhy.cc
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#include "BsPhy.h"
#include "KoiData_m.h"

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
	else if(msg->arrivedOn("fromMac"))  {
		//forward the packet to the channel of the ms
		KoiData *packet = dynamic_cast<KoiData*>(msg);
		send(packet, "toChannel", packet->getDest());
	}
	else if(msg->arrivedOn("fromChannel"))  {
		//forward the packet from the bs channel to the bs mac
		//ev << "Forwarding packet bundle to BsMac" << endl;
		send(msg, "toMac");
	}
}
