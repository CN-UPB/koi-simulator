/**
 * @file ExpChannel.h
 * @class ExpChannel ExpChannel.h
 * @brief This channel draws coefficients from an exponential distribution.
 */

#include "Channel.h"
#include "includes.h"
#include "Position.h"
#include <fstream>
#include <vector>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/exponential_distribution.hpp>

class ExpChannel: public Channel{
	private:
		double expMean;
		double plExp;
		omnetpp::simtime_t initOffset;
		std::ofstream downValues;
		std::ofstream upValues;
		boost::random::mt19937 randEng;
		boost::random::exponential_distribution<double> distExp;

		double pathgain(Position sender, Position receiver);
		void recomputeCoefficients(
				const std::vector<std::vector<Position>>& msPositions);
	
	public:
		void handleMessage(omnetpp::cMessage* msg);
		bool init(omnetpp::cSimpleModule* module,
				const std::vector<std::vector<Position>>& msPositions, 
				const std::map<int,Position>& neighbourPositions);

		void recomputePerTTIValues();
		void updateChannel(const std::vector<std::vector<Position>>& msPos);

		virtual ~ExpChannel();
};
