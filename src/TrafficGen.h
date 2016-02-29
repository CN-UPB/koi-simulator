/**
 * \file TrafficGen.h
 *
 * The KOI traffic generator.
 *
 * Creates packets of type KoiData, with the current mobile station as the 
 * source and one of a number of communication partners as the destination.
 * Packets can be generated with a fixed frequency, or interarrival times 
 * drawn from an exponential distribution.
 */

#pragma once

#include <vector>
#include <omnetpp.h>

class TrafficGen: public cSimpleModule{
	private:
		int bsId;
		std::vector<int> commPartners;
		double deadline;
		double initOffset;
		int msId;
		int packetLength;
		double period;
		bool periodicTraffic;
	
	protected:
		virtual void initialize();
		virtual void handleMessage(cMessage *msg);
	
	public:
		~TrafficGen() = default;
};
