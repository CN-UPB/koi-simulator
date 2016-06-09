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

#include "includes.h"
#include <unordered_map>

class TrafficGen: public cSimpleModule{
	private:
		struct StreamDef{
			unsigned long streamId;
			int destBsId;
			int destMsId;
			double period;
			bool d2d;
		};
		int bsId;
		std::unordered_map<unsigned long,StreamDef> streams;
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
