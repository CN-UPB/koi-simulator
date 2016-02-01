/*
 * App.cc
 *
 * Simple App that's generate traffic for the mobile stations
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#include "App.h"
#include "DataPacket_m.h"
#include "DataPacketBundle_m.h"
#include <iostream>

Define_Module(App);

inline simtime_t App::packetDistribution()  {
    return simTime() + par("packetDistribution");
}

void App::initialize()  {
    packetLength = par("packetLength");
    bsId = par("bsId");
    msId = par("msId");
    initOffset = par("initOffset");
    numberOfPackets = 0;
    epsilon = par("epsilon");
    runtime = par("runtime");
    countPackets = false;
    tti = par("tti");

    //self message that is rescheduled every time to produce traffic
    scheduleAt(packetDistribution() + initOffset, new cMessage("GEN_TRAFFIC"));
    //self message which starts the packet counting for eval
    scheduleAt(initOffset + (3 * tti) - (1.5 * epsilon), new cMessage("START_PACKET_COUNTING"));
}

void App::handleMessage(cMessage *msg)  {
    if(msg->isSelfMessage())  {
        //GEN_TRAFFIC
        if(msg->isName("GEN_TRAFFIC"))  {
            // Full Buffer Model: Send One Packet, Mac will ensure that buffer is full with similar packets
			DataPacket *packet = new DataPacket("DATA");
			packet->setByteLength(packetLength);
			packet->setBsId(bsId);
			packet->setMsId(msId);
			send(packet, "toMac");
            scheduleAt(simTime() + tti, msg);
        }
        else if(msg->isName("START_PACKET_COUNTING"))  {
            countPackets = true;
            msg->setName("STOP_PACKET_COUNTING");
            scheduleAt(runtime - 1.5 * epsilon, msg);
        }
        else if(msg->isName("STOP_PACKET_COUNTING"))  {
            countPackets = false;
            delete msg;
        }
    }
    else if(msg->arrivedOn("fromMac"))  { //traffic from MS or BS
        DataPacketBundle *bundle = (DataPacketBundle *) msg;
        //ev << "DataPacketBundle arrived in App MS: " << bundle->getMsId() << " BS: " << bundle->getBsId() << endl;
        if(countPackets)
            numberOfPackets += bundle->getPacketsArraySize();
        delete bundle;
    }
}

App::~App()  {
    //std::cout << "BS " << bsId << "; MS " << msId << ": " << numberOfPackets << " packets!" << std::endl;
}
