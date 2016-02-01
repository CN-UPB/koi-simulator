/*
 * Partition.h
 *
 *  Created on: Jul 31, 2013
 *      Author: Sascha Schmerling
 */

#pragma once

#include <cpartition.h>

class Partition : public cPartition  {
    private:
        int currentTTI;
        simtime_t tti;
        simtime_t epsilon;
        simtime_t initOffset;
        simtime_t nextBarrier;
        bool initPassed;

    public:
        Partition();
        virtual ~Partition();
        simtime_t getNextBarrier(simtime_t nextMsgTime, int partitionId);
        void processIncomingMsg(cMessage *msg);
        void processOutgoingMsg(cMessage *msg);
        void insertMsgMinBarrier(cMessage *msg, simtime_t minBarrier);
        void updateMsgMinBarrier(cMessage *msg, simtime_t minBarrier);
        void deleteMsgMinBarrier(cMessage *msg);
};
