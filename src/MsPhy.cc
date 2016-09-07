/*
 * MsPhy.cc
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#include "MsPhy.h"
#include "KoiData_m.h"
#include "MessageTypes.h"

Define_Module(MsPhy);

void MsPhy::initialize()  {
}

void MsPhy::handleMessage(cMessage *msg)  {
    if(msg->isName("MS_POS_UPDATE"))  {
        send(msg, "toMsChannel");
    }
    //currently it only forward the packets
    else if(msg->arrivedOn("fromMac"))  {
	switch(msg->getKind()){
		case MessageType::transInfo:
			send(msg,"toMsChannel");
			break;
		case MessageType::koidata:
			KoiData *packet = dynamic_cast<KoiData*>(msg);
			switch(packet->getMessageDirection()){
				case MessageDirection::up:
					send(msg, "toChannel");
					break;
				case MessageDirection::d2dDown:
				case MessageDirection::d2dUp:
					send(msg,"toMs",packet->getDest());
					break;
			}
	    }
    }
    else if(msg->arrivedOn("fromChannel"))  {
        //ev << "Forwarding packet/packets to the mac layer" << endl;
        send(msg, "toMac");
    }
}
