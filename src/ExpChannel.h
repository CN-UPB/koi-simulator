/**
 * @file ExpChannel.h
 * @class ExpChannel ExpChannel.h
 * @brief This channel draws coefficients from an exponential distribution.
 */

#include "Channel.h"
#include "includes.h"
#include "Position.h"
#include <vector>

class ExpChannel: public Channel{
	private:
		double expMean;
		double plExp;

		double pathloss(Position sender, Position receiver);
		void recomputeCoefficients(
				const std::vector<std::vector<Position>>& msPositions);
	
	public:
		void handleMessage(cMessage* msg);
		bool init(cSimpleModule* module,
				const std::vector<std::vector<Position>>& msPositions, 
				std::map<int,Position>& neighbourPositions);

		void updateChannel(const std::vector<std::vector<Position>>& msPos);

		virtual ~ExpChannel(){}
};
