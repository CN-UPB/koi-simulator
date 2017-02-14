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
#include <fstream>
#include <mutex>
#include <unordered_map>

class TrafficGen: public omnetpp::cSimpleModule{
	private:
		struct StreamDef{
			StreamDef(unsigned long streamId,
					int destBsId,
					int destMsId,
					double period,
					double deadline,
					bool d2d)
				:streamId(streamId),
				destBsId(destBsId),
				destMsId(destMsId),
				period(period),
				deadline(deadline),
				d2d(d2d)
			{}
			StreamDef() = default;
			unsigned long streamId;
			int destBsId;
			int destMsId;
			double period;
			double deadline;
			bool d2d;
		};
		int bsId;
		std::unordered_map<unsigned long,StreamDef> streams;
		double initOffset;
		double tti;
		int msId;
		int packetLength;
		bool periodicTraffic;
		std::ofstream* delays;
		/**
		 * Flag signaling that comm table loading should only be conducted once
		 */
		static std::once_flag tFlag;
		static omnetpp::cXMLElement* commTable;
		bool d2dActive;
		static std::vector<StreamDef> parseCommTable(int bsId, int msId);
		static void loadComTable(const std::string& fpath);
	
	protected:
		virtual void initialize();
		virtual void handleMessage(omnetpp::cMessage *msg);
	
	public:
		~TrafficGen();
};
