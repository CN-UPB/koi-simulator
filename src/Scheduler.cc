/*
 * Scheduler.cc
 *
 *  Created on: Jul 10, 2014
 *      Author: Thomas Prinz
 */
 
#include "Scheduler.h"

using namespace omnetpp;

void Scheduler::init(cModule *module){}
void Scheduler::setSINR(mat sinr){}
Schedule* Scheduler::calcUpSchedule(cQueue *packetQueue, std::vector<int> *transmitRequests){
    return new Schedule;
}
std::vector<int> Scheduler::calcDownSchedule(cQueue *packetQueue){
    std::vector<int> result;
    return result;
}
Scheduler::~Scheduler(){}
