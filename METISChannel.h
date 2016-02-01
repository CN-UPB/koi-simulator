 /**
 * @file   METISChannel.h
 * @Author Thomas Prinz (thomas.prinz@rwth-aachen.de)
 * @date   Januar, 2015
 * @brief  METIS Channel Subclass header file
 *
 * Header file that contains the METIS Channel subclass
 */
 
#pragma once

#include <omnetpp.h>
#include <unordered_map>
#include <cmath>
#include "Position.h"
#include "NeighbourIdMatching.h"
#include "Channel.h"
#include <algorithm>


class METISChannel : public Channel{
	private:
		// METIS Channel Parameters
		double freq_c;						/*!< center/carrier frequence */
		double heightUE;					/*!< Height of the user equipments */
		double heightBS;					/*!< Height of the base stations */
		double speedOfLight;					/*!< speed of light value */
		double xPos;						/*!< BS Position x value */
		double yPos;						/*!< BS Position y value */
		double **sigma_ds_LOS;					/*!< Delay Spread LOS sigma */
		double **sigma_asD_LOS;					/*!< an integer value */
		double **sigma_asA_LOS;					/*!< an integer value */
		double **sigma_zsD_LOS;					/*!< an integer value */
		double **sigma_zsA_LOS;					/*!< an integer value */
		double **sigma_sf_LOS;					/*!< Shadow Fading LOS sigma */
		double **sigma_kf_LOS;					/*!< an integer value */
		double **sigma_ds_NLOS;					/*!< Delay Spread NLOS sigma */
		double **sigma_asD_NLOS;				/*!< an integer value */
		double **sigma_asA_NLOS;				/*!< an integer value */
		double **sigma_zsD_NLOS;					/*!< an integer value */
		double **sigma_zsA_NLOS;					/*!< an integer value */
		double **sigma_sf_NLOS;					/*!< Shadow Fading NLOS sigma */
		//double **autoCorrelation;				/*!< The spatial correlation values of the MS links. */
		double **autoCorrelation_LOS;				/*!< The spatial correlation values of the MS for LOS links. */
		double **autoCorrelation_NLOS;				/*!< The spatial correlation values of the MS for NLOS links. */
		double *bs_antenna_bearing;				/*!< bearing angles of the 3 BS sectors */
		double *bs_antenna_downtilt;				/*!< downtilt angles of the 3 BS sectors */
		double *bs_antenna_slant;				/*!< Slant angles of the 3 BS sectors */
		double *ms_antenna_bearing;				/*!< bearing angle of the MS */
		double *ms_antenna_downtilt;				/*!< downtilt angle of the MS */
		double *ms_antenna_slant;				/*!< Slant angle of the MS */
		double ***clusterDelays;				/*!< The delays in seconds for each cluster [#MS]x[#BS]x[#Cluster] */
		double ***clusterDelays_LOS;				/*!< The delays in seconds for each cluster in LOS case [#MS]x[#BS]x[#Cluster] */
		double ***clusterPowers;				/*!< The cluster power distribution for each link [#MS]x[#BS]x[#Cluster] */
		double ***rayPowers;					/*!< The ray power for each cluster */
		double ****azimuth_ASA;					/*!< The angles of arrival in azimuth plane */
		double ****azimuth_ASD;					/*!< The angles of departure in azimuth plane */
		double ****elevation_ASA;				/*!< The angles of arrival in azimuth plane */
		double ****elevation_ASD;				/*!< The angles of departure in azimuth plane */
		double *****randomPhase;				/*!< The random subpath phases uniformly from [0,2pi) */
		double **randomPhase_LOS;					/*!< The random phase for the LOS component, taken uniformly from [0,2pi) */
		double ****Xn_m;					/*!< Cross polarization values per ray */
		double **AoA_LOS_dir;					/*!< azimuth angle of arrival of true geometric LOS direction */
		double **AoD_LOS_dir;					/*!< azimuth angle of departure of true geometric LOS direction */
		double **ZoA_LOS_dir;					/*!< zenith angle of arrival of true geometric LOS direction */
		double **ZoD_LOS_dir;					/*!< zenith angle of departure of true geometric LOS direction */
		int timeSamples;					/*!< Number of TTIs until Position is updated (Number of Time Samples for Channel Model) */
		double **timeVector;					/*!< Time vector für scm computation */
		double ***channelGain;					/*!< Final channel gain within time axis */
		bool **LOSCondition;					/*!< Stores wether each of the links is LOS or NLOS */
		int numberOfMobileStations;				/*!< Number of MSs within this BS */
		int bsId;						/*!< Unique ID of according BS */
		double tti;						/*!< Transmission Time Interval */
		Position *MSPos;					/*!< Position of the MS */
		double *MSVelMag;					/*!< Magnitude of MS Velocity */
		double **MSVelDir;					/*!< Direction of MS Velocity */
		int maxNumberOfNeighbours;				/*!< Max Number of Neighbour BS that interfere with this one */
		map <int,Position> neighbourPos;			/*!< Positions of Neighbour BS */
		cSimpleModule *initModule;				/*!< Pointer to OMNeT module for intermodule communication */
		double ***SINRtable;					/*!< Table to save precomputed channel values (linear) */
		double ****SINRneighbour;				/*!< Table to save precomputed neighbourvalues values (linear) */
		int upRBs;						/*!< Number of up resource blocks*/
		int downRBs;						/*!< Number of down resource blocks */
		int SINRcounter;					/*!< If position resend intervall > 1, it counts the current TTI */
		int NumTxAntenna;					/*!< Number of Transmitter Antenna */
		int NumRxAntenna;					/*!< Number of Receiver Antenna */
		double ***RxAntennaPosition;				/*!< Position vector of Receiver antenna */
		double ***TxAntennaPosition;				/*!< Position vector of Transmitter antenna */
		int numOfInterferers;					/*!< Number of actual interferers, based on network layout and neighbour distance */
		
		//! Calculates the pathloss for a given distance.
		double CalcPathloss(double dist2D, double dist3D, bool LOS);
		
		//! Maps the cluster number to scaling factors for azimuth angles
		double C_AS(int numCluster, bool LOS, int i);
		
		//! Maps the cluster number to scaling factors for zenith angles
		double C_ZS(int numCluster, bool LOS);
		
		//! Calculates the probability of a link being a LOS link.
		bool LineOfSight(double dist2D);
		
		//! Calculate the mean of Zenith spread of departure.
		double mean_ZSD(double dist2D, double heightUE, bool LOS);
		
		//! Calculate the sigma of Zenith spread of departure.
		double sigma_ZSD(double meanZSD, bool LOS);

		//! Generate the spatial correlation between the MS for LOS links.
		void generateAutoCorrelation_LOS();

		//! Generate the spatial correlation between the MS for NLOS links.
		void generateAutoCorrelation_NLOS();
        
	public:
		//! Constructor of METIS Channel subclass.
		METISChannel(){bsId = -1;}
		
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
		~METISChannel();
};
