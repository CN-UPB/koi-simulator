/*
 * StopWatch.h
 *
 *  Created on: Jul 31, 2013
 *      Author: Sascha Schmerling
 */

#pragma once

#include "includes.h"
#include <cstopwatch.h>

class StopWatch : public cSimpleModule {
    private:
        cStopWatch watch;
        char result[100];

        double tti;
        double initOffset;
        double epsilon;
        double runtime;

    public:
        StopWatch();
        virtual ~StopWatch();
        virtual void initialize();
        virtual void handleMessage(cMessage *msg);
};
