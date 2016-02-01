/*
 * StopWatch.cc
 *
 *  Created on: Jul 31, 2013
 *      Author: Sascha Schmerling
 */

#include "StopWatch.h"
#include <iostream>
#include <fstream>

Define_Module(StopWatch);

StopWatch::StopWatch() {
}

void StopWatch::initialize()  {
    epsilon = par("epsilon");
    runtime = par("runtime");
    tti = par("tti");
    initOffset = par("initOffset");

    cMessage *selfMsg = new cMessage("START_WATCH");
    scheduleAt(initOffset + (5 * tti) - (1.5 * epsilon), selfMsg);
}

void StopWatch::handleMessage(cMessage *msg)  {
    if(msg->isSelfMessage())  {
        if(msg->isName("START_WATCH"))  {
            msg->setName("STOP_WATCH");
            scheduleAt(runtime - 1.5 * epsilon, msg);
            watch.start();
        }
        else if(msg->isName("STOP_WATCH"))  {
            watch.stop();
            watch.getTime(result, 100);
            delete msg;
        }
    }
}

StopWatch::~StopWatch()  {
    std::cout << "StopWatch result: " << result << std::endl;
    // To get times when run several times via script
	std::ofstream myfile;
	myfile.open ("stopwatch.txt", std::ofstream::app);
	myfile << result << "\n";
	myfile.close();
}
