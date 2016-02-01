/*
 * Partition.cc
 *
 *  Created on: Jul 31, 2013
 *      Author: Sascha Schmerling
 */

#include "Partition.h"
#include <iostream>

//used by interval based sync
Register_Class(Partition); //only one partition per model

Partition::Partition()  {
    currentTTI = 0;
    tti = 1e-03; //FIXME get out of config
    epsilon = 1e-09; //FIXME get out of config
    initOffset = 1; //FIXME get out of config
    initPassed = false;
    nextBarrier = 0;
}

Partition::~Partition()  {
}

simtime_t Partition::getNextBarrier(simtime_t nextMsgTime, int partitionId)  {
    if(nextMsgTime<=initOffset /*!initPassed*/)  {
        initPassed = true;
        return nextMsgTime; //BS position exchanges at init time 0
    }
    else  {
        if(currentTTI == 0)  {
            currentTTI++; //first real simulation time; first position exchanges from the MSs
            return initOffset;
        }
        else  {
            if(nextMsgTime > nextBarrier)  {
                currentTTI++;
                nextBarrier = initOffset + (currentTTI * tti) - (2 * epsilon);
            }
            assert(nextBarrier >= simTime());
            return nextBarrier;
        }
    }
}

void Partition::processIncomingMsg(cMessage *msg)  {
    return;
}

void Partition::processOutgoingMsg(cMessage *msg)  {
    return;
}

void Partition::insertMsgMinBarrier(cMessage *msg, simtime_t minBarrier)  {
    return;
}

void Partition::updateMsgMinBarrier(cMessage *msg, simtime_t minBarrier)  {
    return;
}

void Partition::deleteMsgMinBarrier(cMessage *msg)  {
    return;
}
