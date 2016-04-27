 /**
 * @file   METISChannel.h
 * @Author Thomas Prinz (thomas.prinz@rwth-aachen.de)
 * @date   Januar, 2015
 * @brief  METIS Channel Subclass header file
 *
 * Header file that contains the METIS Channel subclass
 */
 
#pragma once

#include "includes.h"
#include <unordered_map>
#include <cmath>
#include "Position.h"
#include "NeighbourIdMatching.h"
#include "Channel.h"
#include <algorithm>
#include <vector>

using std::vector;

class METISChannel : public Channel{
	private:
		// METIS Channel Parameters
		double freq_c;						/*!< center/carrier frequence */
		double heightUE;					/*!< Height of the user equipments */
		double heightBS;					/*!< Height of the base stations */
		static constexpr double speedOfLight = 299792458.0;					/*!< speed of light value */
		double xPos;						/*!< BS Position x value */
		double yPos;						/*!< BS Position y value */
		int N_cluster_LOS;
		int N_cluster_NLOS;
		int numOfRays_LOS;
		int numOfRays_NLOS;
		int timeSamples;					/*!< Number of TTIs until Position is updated (Number of Time Samples for Channel Model) */
		double **timeVector;					/*!< Time vector für scm computation */
		double ***channelGain;					/*!< Final channel gain within time axis */
		int numberOfMobileStations;				/*!< Number of MSs within this BS */
		int bsId;						/*!< Unique ID of according BS */
		double tti;						/*!< Transmission Time Interval */
		int maxNumberOfNeighbours;				/*!< Max Number of Neighbour BS that interfere with this one */
		map <int,Position> neighbourPositions;			/*!< Positions of Neighbour BS */
		NeighbourIdMatching *neighbourIdMatching;
		cSimpleModule *initModule;				/*!< Pointer to OMNeT module for intermodule communication */
		double ***SINRtable;					/*!< Table to save precomputed channel values (linear) */
		double ****SINRneighbour;				/*!< Table to save precomputed neighbourvalues values (linear) */
		int upRBs;						/*!< Number of up resource blocks*/
		int downRBs;						/*!< Number of down resource blocks */
		int SINRcounter;					/*!< If position resend intervall > 1, it counts the current TTI */
		int NumBsAntenna;					/*!< Number of Base Station Antenna */
		int NumMsAntenna;					/*!< Number of Mobile Station Antenna */
		double ***OwnBsAntennaPosition;				/*!< Position vector of Base Station antenna */
		int numOfInterferers;					/*!< Number of actual interferers, based on network layout and neighbour distance */
		double vel;
		double XPR_Mean_LOS;
		double XPR_Std_LOS;
		double XPR_Mean_NLOS;
		double XPR_Std_NLOS;
		bool initialized;					/*!< True iff METISChannel::init has been called */
		
		//! Calculates the pathloss for a given distance.
		double CalcPathloss(double dist2D, double dist3D, bool LOS);
		
		//! Maps the cluster number to scaling factors for azimuth angles
		double C_AS(int numCluster, bool LOS, int i,
				vector<vector<double>> sigma_kf_LOS);
		
		//! Maps the cluster number to scaling factors for zenith angles
		double C_ZS(int numCluster, bool LOS,
				vector<vector<double>> sigma_kf_LOS);
		
		//! Calculates the probability of a link being a LOS link.
		bool LineOfSight(double dist2D);
		
		//! Calculate the mean of Zenith spread of departure.
		double mean_ZSD(double dist2D, double heightUE, bool LOS);
		
		//! Calculate the sigma of Zenith spread of departure.
		double sigma_ZSD(double meanZSD, bool LOS);

		//! Recalculate the METIS large scale params for the given Rx/Tx
		void recomputeLargeScaleParameters(const vector<Position>& senders,
				const vector<Position>& receivers, 
				vector<vector<double>>& sigma_ds_LOS,
				vector<vector<double>>& sigma_asD_LOS,
				vector<vector<double>>& sigma_asA_LOS,
				vector<vector<double>>& sigma_zsD_LOS,
				vector<vector<double>>& sigma_zsA_LOS,
				vector<vector<double>>& sigma_sf_LOS,
				vector<vector<double>>& sigma_kf_LOS,
				vector<vector<double>>& sigma_ds_NLOS,
				vector<vector<double>>& sigma_asD_NLOS,
				vector<vector<double>>& sigma_asA_NLOS,
				vector<vector<double>>& sigma_zsD_NLOS,
				vector<vector<double>>& sigma_zsA_NLOS,
				vector<vector<double>>& sigma_sf_NLOS,
				vector<vector<double>>& sigma_kf_NLOS);

		//! Recalculate all position dependent values, e.g. SINR
		void recomputeMETISParams(Position **msPositions);

		//! Generate the spatial correlation between the MS for LOS links.
		void generateAutoCorrelation_LOS(const vector<Position>& senders,
				const vector<Position>& receivers,
				vector<vector<vector<double>>>& correlation);

		//! Generate the spatial correlation between the MS for NLOS links.
		void generateAutoCorrelation_NLOS(const vector<Position>& senders,
				const vector<Position>& receivers,
				vector<vector<vector<double>>>& correlation);
        
	public:
		//! Constructor of METIS Channel subclass.
		METISChannel(){
			bsId = -1;
			initialized = false;
		}
		
		//! Initialize the METIS channel through ini access via OMNeT++ module pointer.
		bool init(cSimpleModule* module, Position** msPositions, std::map <int,Position> neighbourPositions);

		//! Generate the channel coefficients every positionResendInterval
		//void METISChannel::calcChannel_METIS();
		
		//! Allows the OMNeT++ module to pass messages to this METIS channel class.
		void handleMessage(cMessage* msg);
		
		//! Computes the pathloss for a given distance using an arbitrary model.
		double calcPathloss(double dist);
		
		//! Computes the Termal Noise ("johnson nyquist noise")
		double getTermalNoise(double temp, double bandwidth);
		
		//! Calculates the current SINR for given interferers and given RB.
		double calcSINR(int RB, vector<double> &power, vector<Position> &pos, vector<int> &bsId_, bool up, int msId);
		
		//! Calculates the current SINR for given interferers and all RB.
		vec calcSINR(vector<double> &power, vector<Position> &pos, vector<int> &bsId_, bool up, int msId);
		
		//! Updates the MS position if velocity > 0. The interval in which the postion is updated can be set within omnet.ini
		void updateChannel(Position** msPos);
		
		//! Destructor of METIS Channel subclass.
		virtual ~METISChannel();
};
