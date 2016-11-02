/**
 * @file FactoryChannel.h
 * @class FactoryChannel FactoryChannel.h
 * @brief Channel model based on Tanghe et al.
 */

#include "Channel.h"
#include "includes.h"
#include "Position.h"
#include <fstream>
#include <vector>

class FactoryChannel: public Channel{
	private:
		double expMean;
		double plExp;
		double d0;
		double pl0;
		simtime_t initOffset;
		std::ofstream downValues;
		std::ofstream upValues;
		std::vector<std::vector<Position>> msPos;

		double pathgain(Position sender, Position receiver);
		double shadowfading();
		void recomputeCoefficients(
				const std::vector<std::vector<Position>>& msPositions);
	
	public:
		void clearTransInfo();
		void handleMessage(cMessage* msg);
		bool init(cSimpleModule* module,
				const std::vector<std::vector<Position>>& msPositions, 
				std::map<int,Position>& neighbourPositions);

		void updateChannel(const std::vector<std::vector<Position>>& msPos);

		virtual ~FactoryChannel();
};
