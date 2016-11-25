/**
 * @file FactoryChannel.h
 * @class FactoryChannel FactoryChannel.h
 * @brief Channel model based on Tanghe et al.
 */

#include "Channel.h"
#include "includes.h"
#include "Position.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/normal_distribution.hpp>

#include <fstream>
#include <vector>

class FactoryChannel: public Channel{
	private:
		double plExp;
		double expMean;
		double d0;
		double pl0;
		double kMean;
		double kSigma;
		double shSigma;
		double bsGain;
		double msGain;
		double transPower;
		simtime_t initOffset;
		std::ofstream downValues;
		std::ofstream upValues;
		std::vector<std::vector<Position>> msPos;
		boost::random::mt19937 randEng;
		boost::random::exponential_distribution<double> distExp;
		boost::random::normal_distribution<double> distNorm;
		/**
		 * Shadowing table UP, [BS]x[MS]
		 */
		std::vector<std::vector<double>> shUp;
		/**
		 * Shadowing table DOWN, [MS]x[BS]
		 */
		std::vector<std::vector<double>> shDown;
		/**
		 * Shadowing table D2D DOWN, [CELL]x[MS]x[MS]
		 */
		std::vector<std::vector<std::vector<double>>> shD2DDown;
		/**
		 * Shadowing table D2D UP, [CELL]x[MS]x[MS]
		 */
		std::vector<std::vector<std::vector<double>>> shD2DUp;

		double pathgain(Position sender, Position receiver);
		double fadingRicean(double pl, double gainTx, double gainRx);
		double fadingExponential();
		double shadowing();
		void recomputeCoefficients(
				const std::vector<std::vector<Position>>& msPositions);
		void generateShadowing(
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
