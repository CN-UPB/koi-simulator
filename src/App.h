/*
 * App.cc
 *
 * Simple App that's generate traffic for the mobile stations
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#pragma once

#include "includes.h"

class App : public cSimpleModule  {
    private:
        int packetLength;
        int bsId;
        int msId;
        double initOffset;
        int numberOfPackets;
        bool countPackets;
        simtime_t runtime;
        simtime_t epsilon;
        simtime_t tti;

        inline simtime_t packetDistribution();

    protected:
        virtual void initialize();
        virtual void handleMessage(cMessage *msg);

    public:
        ~App();
};
