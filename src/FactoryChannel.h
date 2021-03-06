/**
 * @file FactoryChannel.h
 * @class FactoryChannel FactoryChannel.h
 * @brief Channel model based on Tanghe et al.
 */

#include "Channel.h"
#include "includes.h"
#include "Position.h"
#include "VecNd.h"

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
		omnetpp::simtime_t initOffset;
		std::ofstream downValues;
		std::ofstream upValues;
		boost::random::mt19937 randEng;
		boost::random::exponential_distribution<double> distExp;
		boost::random::normal_distribution<double> distNorm;
		/**
		 * Shadowing table UP, [BS]x[MS]
		 */
		VectorNd<double,2> shUp;
		/**
		 * Shadowing table DOWN, [MS]x[BS]
		 */
		VectorNd<double,2> shDown;
		/**
		 * Shadowing table D2D DOWN, [CELL]x[MS]x[MS]
		 */
		VectorNd<double,3> shD2DDown;
		/**
		 * Shadowing table D2D UP, [CELL]x[MS]x[MS]
		 */
		VectorNd<double,3> shD2DUp;

		double pathgain(Position sender, Position receiver);
		double fadingExponential();
		double shadowing();
		void recomputeCoefficients(
				const std::vector<std::vector<Position>>& msPositions);
		void generateShadowing(
				const std::vector<std::vector<Position>>& msPositions);
	
	public:
		void handleMessage(omnetpp::cMessage* msg) override;
		bool init(omnetpp::cSimpleModule* module,
				const std::vector<std::vector<Position>>& msPositions, 
				const std::map<int,Position>& neighbourPositions) override;

		void recomputePerTTIValues() override;
		void updateChannel(const std::vector<std::vector<Position>>& msPos) override;
		double calcLongtermUpSINR(int rb, int msId, double transPower) override;
		double calcLongtermDownSINR(int rb, int msId, double transPower) override;
		~FactoryChannel() override;
};
