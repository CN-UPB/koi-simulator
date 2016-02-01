/*
 * proportionalFairScheduler.cc
 *
 *  Created on: Jul 10, 2014
 *      Author: Thomas Prinz
 */
 
#include "proportionalFairScheduler.h"
#include <time.h>
#include <stdlib.h>

using namespace std;
using namespace itpp;

void proportionalFairScheduler::init(cModule *module){
	//std::cout << "Init Prop Fair Scheduler test" << std::endl;
	numberOfRessourceBlocks = module->par("upResourceBlocks");
	numberOfMobileStations = module->par("numberOfMobileStations");
	windowSize = module->par("SlidingWindowSize");
	smoothingFactor = module->par("SmoothingFactor");
	useSinr = module->par("useSinr");
	bsId = module->par("bsId");
	
	string t = module->par("AverageType");
	if(t.compare("Exp Smoothing") == 0){
		type = EXP_AVERAGE;
	}else if(t.compare("Total Average") == 0){
		type = TOTAL_AVERAGE;
	}else if(t.compare("Sliding Window") == 0){
		type = SLIDING_WINDOW;
	}else{
		//default
		type = SLIDING_WINDOW;
	}
	
	uniform = module->par("uniform");
	offset = module->par("offset");
	N_s = module->par("N_s");
	N_c = module->par("N_c");
	runtime = module->par("runtime");
	
	// Used to distinguish output from different Runs
	t = static_cast<long long> (time(NULL));
	
    sinr_ = zeros(numberOfMobileStations,numberOfRessourceBlocks);
    
    averageSinr = ones(numberOfMobileStations,numberOfRessourceBlocks);
    averageRate = ones(numberOfMobileStations,numberOfRessourceBlocks);
    averageTotal = zeros(numberOfMobileStations,numberOfRessourceBlocks);
    SinrSum = zeros(numberOfMobileStations,numberOfRessourceBlocks);
    totalRate = new double[numberOfMobileStations];
    
    for(int i = 0; i < numberOfMobileStations; i++){
		totalRate[i] = 0;
	}
	
	for(int i = 0; i < windowSize; i++){
		sinrHistory.push_back(zeros(numberOfMobileStations,numberOfRessourceBlocks));
	}
}
		
void proportionalFairScheduler::setSINR(mat sinr){
	sinr_ = sinr;
	//std::cout << sinr << std::endl;
	//std::cout << "SINR set: " << sinr << std::endl;
	//std::cout << "SINR set: " << sinr_ << std::endl;
	rate = SINR_to_datarate(sinr);
}

mat proportionalFairScheduler::SINR_to_datarate(mat sinr){
	// Simplification for now. (use util.h for real spectral efficiency)
	return itpp::log2(1 + sinr);
}

// Returns the average according to the chosen type
double proportionalFairScheduler::getAverage(int msID, int rbID){
	switch(type){
		case EXP_AVERAGE:
			return averageSinr(msID,rbID);
		case SLIDING_WINDOW:
			return (SinrSum(msID,rbID) / windowSize);
		case TOTAL_AVERAGE:
			return (averageTotal(msID,rbID) / windowSize);
		default:
			return 1;
	}
}

void proportionalFairScheduler::updateSINR(){	
	
	// Update smoothed Average(SINR):
	averageSinr = (1 - smoothingFactor) * averageSinr + smoothingFactor * sinr_;
	
	// Update smoothed Average(RATE):
	//averageRate = (1 - smoothingFactor) * averageRate + smoothingFactor * SINR_to_datarate(sinr_);
	
	// Update Sliding Window:
	mat sinr_new = inv_dB(sinr_);
	//std::cout << sinr_new << std::endl;
	SinrSum += sinr_new;
	SinrSum -= sinrHistory.front();
	sinrHistory.pop_front();
	sinrHistory.push_back(sinr_new);
	//std::cout << "SINR Updated! SINR sum: " << SinrSum << std::endl;
	
	//cout << sinr_ << endl;
	
	// Update Total Average (TODO: actual transferred data):
	averageTotal = averageTotal + sinr_;
}

Schedule* proportionalFairScheduler::calcUpSchedule(cQueue *packetQueue, vector<int> *transmitRequests){
	//std::cout << "Up Schedule Computation (Proportional Fair (Smoothed Average/SINR)): " << std::endl;
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

		// TODO: Real Datarate insteat of SINR / Add Sliding Window variant
		// Proportional Fair (Smoothed Average/SINR)
		
		double max, maxIdx;
		max = 0;
		maxIdx = 0;
		for(uint j = 0; j < sendingMs.size(); j++){
			double prio = sinr_(i,j) / getAverage(j,i);
			if(prio > max){
				max = prio;
				maxIdx = j;
			}
		}
		int toMsId = sendingMs.at(maxIdx);
		//assign the ms the slot in the schedule
		schedule->setUpSchedule(i, toMsId);
		numberOfPacketsToReceive--; //one packets less to receive
		transmitRequests->at(toMsId)--; //one packet less in the queue
		//if the queue is now empty remove the ms from sending ms
		if(transmitRequests->at(toMsId) == 0)  {
			sendingMs.erase(sendingMs.begin() + maxIdx);
			numberOfSendingMs--;
		}
	}
	return schedule;
}

std::vector<int> proportionalFairScheduler::calcDownSchedule(cQueue *packetQueue){
	//std::cout << "Down Schedule Computation (Proportional Fair (Smoothed Average/SINR)): " << std::endl;
	vector<int> packetCount(numberOfMobileStations);
	vector<int> downSchedule(numberOfRessourceBlocks,-1);
	vector<int> sendingApps(numberOfMobileStations);
	
	ofstream output_datarate;
	static string out = "/home/sahari/updated_METIS/abstractLTEChannelModel/results/BS_" + std::to_string( (long long) bsId ) + "_Schedule_" + std::to_string( t ) + ".txt";
	output_datarate.open (out, fstream::app);
	
    //save the length of the queues; and the total number of packets
    /*
    for(int i = 0; i < numberOfMobileStations; ++i)  {
        int length = packetQueue[i].length();
        numberOfPacketsToSend += length;
        packetCount.at(i) = length;

        if(length != 0)  {
            sendingApps.at(numberOfSendingApps) = i;
            numberOfSendingApps++;
        }
    }*/
    
	for(int i = 0; i < numberOfRessourceBlocks; ++i)  {
		double max, maxIdx;
		max = 0;
		maxIdx = 0;
		double prio = 0;
		for(int j = 0; j < numberOfMobileStations; j++){
			// DB Scale, therefore subtract instead of divide?!
			prio = inv_dB(sinr_(j,i)) / getAverage(j,i);
			//std::cout << "Prio: " << prio << std::endl;
			//cout << "SINR: " << SinrSum(i,j) << " for MS: " << j << endl;
			//cout << "Scaled SINR Value: " << prio << endl;
			if(prio > max){
				max = prio;
				maxIdx = j;
			}
		}
		int toMsId = maxIdx;
		downSchedule.at(i) = toMsId;
		
		//std::cout << "RB: " << i << " Scheduled MS: " << toMsId << std::endl;
		
		//output << "RB: " << i << " SINR: " << sinr_(maxIdx,i) << " Average SINR: " << getAverage(maxIdx,i) << " Scaled SINR: " << max << " MS: " << toMsId << " Datarate: " << datarate << " Simtime: " << simTime() << std::endl;
	}
		// Sum up total rate (uniform case)
		if(uniform){
			for(int i = 0; i < numberOfMobileStations; i++){
				int rbSum = 0;
				int lowestCQI = 15;
				
				// Find the number of ressource blocks scheduled and
				// the CQI for the RB with the lowest CQI
				for(int j = 0; j < numberOfRessourceBlocks; j++){
					if(downSchedule.at(j) == i){
						rbSum++;
						if(lowestCQI > SINR_to_CQI(sinr_(i,j))){
							lowestCQI = SINR_to_CQI(sinr_(i,j));
						}
					}
				}
				
				// start recording after offset
				if(simTime() >= offset){
					datarate = rbSum * getSpectralEfficiency(lowestCQI) * N_s * N_c;

					//std::cout << "MS: " << i << " Efficiency: " << getSpectralEfficiency(lowestCQI) << " Lowest CQI: " << lowestCQI << " N_s: " << N_s << " N_c: " << N_c << " RB Sum: " << rbSum << " Datarate: " << datarate << std::endl;
					output_datarate << "MS: " << i << " Datarate: " << datarate << " Simtime: " << simTime() << " Counter: " << std::round( simTime().dbl() *1000.0*4.0) << std::endl;
					totalRate[i] += datarate;
					
				}
			}
		// Sum up total rate (non-uniform case; each RB scheduled separately)
		}else{
			for(int i = 0; i < numberOfRessourceBlocks; i++){
				if(simTime() >= offset){
					datarate = getSpectralEfficiency(SINR_to_CQI(sinr_(downSchedule.at(i),i))) * N_s * N_c;
					output_datarate << "RB: " << i << " Datarate: " << datarate << " Simtime: " << simTime() << std::endl;
					totalRate[downSchedule.at(i)] += datarate;
				}
			}
		}

		/*
		numberOfPacketsToSend--; //one packets less to send
		packetCount.at(toMsId)--; //one packet less in the queue
		if(packetCount.at(toMsId) == 0)  {
			sendingApps.erase(sendingApps.begin()+ maxIdx);
			numberOfSendingApps--;
			cout << "deleted Sending App " << maxIdx << endl;
		}	
		* */
	output_datarate.close();
	return downSchedule;
}

proportionalFairScheduler::~proportionalFairScheduler(){
	string filename = "/home/sahari/updated_METIS/abstractLTEChannelModel/results/BS_" + std::to_string( (long long) bsId ) + "_Totalrate.txt";
	
	ofstream output(filename);
	if( !output.is_open() ){
		std::cout << "WARNING: Unable to write totalrate to file (cannot open file)" << std::endl;	
	}else{
		double realRuntime = runtime - offset;
		for(int i = 0; i < numberOfMobileStations; i++){
			output << totalRate[i] / realRuntime << endl;
		}
	}
	output.close();
}
