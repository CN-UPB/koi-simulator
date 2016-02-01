/*
 * bestCQIScheduler.cc
 *
 *  Created on: Jul 10, 2014
 *      Author: Thomas Prinz
 */
 
#include "bestCQIScheduler.h"

using namespace std;
using namespace itpp;

void bestCQIScheduler::init(cModule *module){
	std::cout << "Init Random Scheduler test" << std::endl;
	numberOfRessourceBlocks = module->par("upResourceBlocks");
	numberOfMobileStations = module->par("numberOfMobileStations");
}
		
void bestCQIScheduler::setSINR(mat sinr){
	sinr_ = sinr;
}

Schedule* bestCQIScheduler::calcUpSchedule(cQueue *packetQueue, vector<int> *transmitRequests){
	Schedule *schedule = new Schedule("SCHEDULE");
	schedule->setScheduleDirection(0);
	vector<int> packetCount(numberOfMobileStations);
    
	//from the ms to bs
	int numberOfSendingMs = 0;
	int numberOfPacketsToReceive = 0;
	//vector that stores the ids of the sending ms; dynamic size numberOfSendingMs
	vector<int> sendingMs(numberOfMobileStations);
	for(int i = 0; i < numberOfMobileStations; ++i)  {
		numberOfPacketsToReceive += transmitRequests->at(i);
		if(transmitRequests->at(i) != 0)  {
			sendingMs.at(numberOfSendingMs) = i;
			numberOfSendingMs++;
		}
	}
	schedule->setUpScheduleArraySize(numberOfRessourceBlocks);
	schedule->setPowerAdaptationArraySize(numberOfRessourceBlocks);
	for(int i = 0; i < numberOfRessourceBlocks; ++i)  {
		schedule->setUpSchedule(i, -1); //init the slot as free
		schedule->setPowerAdaptation(i, 1.0 / numberOfRessourceBlocks); //power equal for all RBs

		if(numberOfPacketsToReceive == 0)
			continue; //no packets to receive; nothing to do for this slot

		//BEST CQI/SINR
		std::cout << "Up Schedule Computation (Best CQI/SINR): " << std::endl;
		double max, maxIdx;
		max = 0;
		maxIdx = 0;
		for(uint j = 0; j < sendingMs.size(); j++){
			if(sinr_(i,j) > max){
				max = sinr_(i,j);
				maxIdx = j;
			}
		}
		//find a ms for the slot from the packet queues
		int fromMsId = sendingMs.at(maxIdx);
		std::cout << "RB " << i << " is scheduled to MS " << fromMsId << std::endl;
		//assign the ms the slot in the schedule
		schedule->setUpSchedule(i, fromMsId);
		numberOfPacketsToReceive--; //one packets less to receive
		transmitRequests->at(fromMsId)--; //one packet less in the queue
		//if the queue is now empty remove the ms from sending ms
		if(transmitRequests->at(fromMsId) == 0)  {
			sendingMs.erase(sendingMs.begin() + maxIdx);
			numberOfSendingMs--;
		}
	}
	return schedule;
}

std::vector<int> bestCQIScheduler::calcDownSchedule(cQueue *packetQueue){
	vector<int> packetCount(numberOfMobileStations);
	vector<int> downSchedule(numberOfRessourceBlocks);
	vector<int> sendingApps(numberOfMobileStations);
	
	int numberOfPacketsToSend = 0;
    int numberOfSendingApps = 0;
    //save the length of the queues; and the total number of packets
    for(int i = 0; i < numberOfMobileStations; ++i)  {
        int length = packetQueue[i].length();
        numberOfPacketsToSend += length;
        packetCount.at(i) = length;

        if(length != 0)  {
            sendingApps.at(numberOfSendingApps) = i;
            numberOfSendingApps++;
        }
    }
	
	for(int i = 0; i < numberOfRessourceBlocks; ++i)  {
		downSchedule.at(i) = -1; //mark the slot as free

		if(numberOfPacketsToSend == 0)
			continue; //no packets to send; nothing to do

		// BEST CQI/SINR
		std::cout << "Down Schedule Computation (Best CQI/SINR): " << std::endl;
		double max, maxIdx;
		max = 0;
		maxIdx = 0;
		for(uint j = 0; j < sendingApps.size(); j++){
			if(sinr_(i,j) > max){
				max = sinr_(i,j);
				maxIdx = j;
			}
		}
		int toMsId = sendingApps.at(maxIdx);
		downSchedule.at(i) = toMsId;
		std::cout << "RB " << i << " is scheduled to MS " << toMsId << std::endl;
		numberOfPacketsToSend--; //one packets less to send
		packetCount.at(toMsId)--; //one packet less in the queue
		if(packetCount.at(toMsId) == 0)  {
			sendingApps.erase(sendingApps.begin()+ maxIdx);
			numberOfSendingApps--;
		}
	}
	return downSchedule;
}

bestCQIScheduler::~bestCQIScheduler(){}
