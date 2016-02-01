/*
 * proportionalFairScheduler.h
 *
 *  Created on: Jul 10, 2014
 *      Author: Thomas Prinz
 */
 
#pragma once
 
#include <itpp/itbase.h>
#include <omnetpp.h>
#include <vector>
#include "Schedule_m.h"
#include "Scheduler.h"
#include "util.h"

using namespace itpp;
using namespace std;

enum SCHEDULETYPE{EXP_AVERAGE,SLIDING_WINDOW,TOTAL_AVERAGE};

 class proportionalFairScheduler : public Scheduler{
	 protected:
		int numberOfRessourceBlocks;
		int numberOfMobileStations;
		mat sinr_;
		mat rate;
		mat averageSinr;
		mat SinrSum;
		mat averageRate;
		mat averageTotal;
		double *totalRate;
		list<mat> sinrHistory;
		double smoothingFactor;
		int windowSize;
		bool useSinr; // If this is false, datarate is used
		SCHEDULETYPE type;
		int bsId;
		bool uniform;
		double offset;
		int N_s;
		int N_c;
		double runtime;
		double datarate;
		long long t;
		
		mat SINR_to_datarate(mat sinr);
		double getAverage(int msID, int rbID);
		
	 public:
		void init(cModule *module);
		void setSINR(mat sinr);
		void updateSINR();
		Schedule* calcUpSchedule(cQueue *packetQueue, std::vector<int> *transmitRequests);
		std::vector<int> calcDownSchedule(cQueue *packetQueue);
		~proportionalFairScheduler();
 };
