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
#include "Channel.h"
#include "MetisRays.h"
#include "Position.h"
#include "VecNd.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <ostream>
#include <tuple>
#include <unordered_map>
#include <vector>

using std::vector;
using std::array;
using std::tuple;

class METISChannel : public Channel{
	private:
		// METIS Channel Parameters
		vector<Position> bsPositions;
		VectorNd<std::complex<double>,4> delayDownTable;
		VectorNd<std::complex<double>,5> delayUpTable;
		VectorNd<std::complex<double>,5> delayD2DUpTable;
		VectorNd<std::complex<double>,5> delayD2DDownTable;
		double epsilon;
		double freq_c;						/*!< center/carrier frequence */
		double heightUE;					/*!< Height of the user equipments */
		double heightBS;					/*!< Height of the base stations */
		/**
		 * Wavenumber
		 *
		 * 2*PI/(C/Carrier Freq)
		 */
		double k_0;
		VectorNd<bool,2> losDownTable;
		VectorNd<bool,3> losUpTable;
		VectorNd<bool,3> losD2DTable;
		static double ray_offset[20];				/* Ray offset. Table 7.6 METIS D1.2 */
		int N_cluster_LOS;
		int N_cluster_NLOS;
		int numOfRays_LOS;
		int numOfRays_NLOS;
		int NumBsAntenna;					/*!< Number of Base Station Antenna */
		int NumMsAntenna;					/*!< Number of Mobile Station Antenna */
		vector<vector<array<double,3>>> bsAntennaPositions;				/*!< Position vector of Base Station antenna */
		int numOfInterferers;					/*!< Number of actual interferers, based on network layout and neighbour distance */
		VectorNd<Position,2> msPos;
		VectorNd<RayCluster,5> precompDownTable;
		VectorNd<RayCluster,6> precompUpTable;
		VectorNd<RayCluster,6> precompD2DTable;
		/**
		 * MS velocity in m/s
		 */
		double velocity;
		double wavelength;
		double XPR_Mean_LOS;
		double XPR_Std_LOS;
		double XPR_Mean_NLOS;
		double XPR_Std_NLOS;
		
		/**
		 * @brief Calculate Antenna positions for the given transmitters
		 */
		vector<vector<array<double,3>>> computeAntennaPos(
				const vector<Position>& transmitterPos,
				int numAntennas,
				double heightAntennas);
		
		//! Calculates the pathloss for a given distance.
		double CalcPathloss(double dist2D, double dist3D, bool LOS);
		
		//! Maps the cluster number to scaling factors for azimuth angles
		double C_AS(int numCluster, bool LOS, int i,
				VectorNd<double,2> sigma_kf_LOS);
		
		//! Maps the cluster number to scaling factors for zenith angles
		double C_ZS(int numCluster, bool LOS,
				VectorNd<double,2> sigma_kf_LOS);
		
		//! Calculates the probability of a link being a LOS link.
		bool LineOfSight(double dist2D);
		
		//! Calculate the mean of Zenith spread of departure.
		double mean_ZSD(double dist2D, double heightUE, bool LOS);
		
		//! Calculate the sigma of Zenith spread of departure.
		double sigma_ZSD(double meanZSD, bool LOS);

		//! Recalculate the METIS large scale params for the given Rx/Tx
		void recomputeLargeScaleParameters(const vector<Position>& senders,
				const vector<Position>& receivers, 
				VectorNd<double,2>& sigma_ds_LOS,
				VectorNd<double,2>& sigma_asD_LOS,
				VectorNd<double,2>& sigma_asA_LOS,
				VectorNd<double,2>& sigma_zsD_LOS,
				VectorNd<double,2>& sigma_zsA_LOS,
				VectorNd<double,2>& sigma_sf_LOS,
				VectorNd<double,2>& sigma_kf_LOS,
				VectorNd<double,2>& sigma_ds_NLOS,
				VectorNd<double,2>& sigma_asD_NLOS,
				VectorNd<double,2>& sigma_asA_NLOS,
				VectorNd<double,2>& sigma_zsD_NLOS,
				VectorNd<double,2>& sigma_zsA_NLOS,
				VectorNd<double,2>& sigma_sf_NLOS);

		/**
		 * @brief Precompute all METIS values for clusters/delays
		 *
		 * Cluster ray sum components as well as delays can be precomputed 
		 * at the start of each simulation, instead of being recomputed for each 
		 * TTI.
		 */
		void precomputeMETISValues(const vector<vector<Position>>& msPositions);

		//! Generate the spatial correlation between the MS for LOS links.
		void generateAutoCorrelation_LOS(const vector<Position>& senders,
				const vector<Position>& receivers,
				VectorNd<double,3>& correlation);

		//! Generate the spatial correlation between the MS for NLOS links.
		void generateAutoCorrelation_NLOS(const vector<Position>& senders,
				const vector<Position>& receivers,
				VectorNd<double,3>& correlation);

		/**
		 * @brief Compute line of sight for each sender/receiver pair
		 *
		 * For each pair of sender/receiver, the METIS computations 
		 * need the information on whether there is a line of sight 
		 * between them or not. The probabilities are taken from 
		 * section 7.3.10 of METIS D1.2.
		 *
		 * @param sendPos The 2D positions of the senders
		 * @param receivePos The 2D positions of the receivers
		 * @return A two dimensional [#receivers]x[#senders] vector
		 * 		where [i][j] is `TRUE` iff there is a line of 
		 * 		sight between receiver `i` and sender `j`
		 */
		VectorNd<bool,2> genLosCond(const vector<Position>& sendPos,
				const vector<Position>& receivePos);
        
		/**
		 * @brief Recompute per-cluster powers
		 *
		 */
		VectorNd<double,3> genClusterPowers(const VectorNd<bool,2>& LOSCondition,
				const VectorNd<double,3>& clusterDelays,
				const VectorNd<double,2>& sigmaDS_LOS,
				const VectorNd<double,2>& sigmaDS_NLOS,
				const VectorNd<double,2>& sigmaKF_LOS
				);

		/**
		 * @brief Recompute per-ray powers
		 */
		VectorNd<double,3> recomputeRayPowers(const VectorNd<bool,2>& LOSCondition,
				VectorNd<double,3>& clusterPowers
				);

		/**
		 * @brief Recompute angle directions
		 */
		tuple<VectorNd<double,2>,VectorNd<double,2>,VectorNd<double,2>,VectorNd<double,2>> 
		recomputeAngleDirection(
				const vector<Position>& receivers,
				const vector<Position>& senders,
				double heightSenders,
				double heightReceivers
				);
				
		/**
		 * @brief Recompute ray azimuth angles
		 *
		 * This method works for angles of arrival as well as angles 
		 * of departure, depending on which angle spread and angle 
		 * direction values are provided.
		 */
		VectorNd<double,4>
		recomputeAzimuthAngles(const VectorNd<bool,2>& LOSCondition,
				const VectorNd<double,2>& sigma_as_LOS,
				const VectorNd<double,2>& sigma_as_NLOS,
				const VectorNd<double,2>& sigma_kf,
				const VectorNd<double,3>& clusterPowers,
				const VectorNd<double,2>& angleDir,
				const bool arrival
				);

		/**
		 * @brief Recompute ray zenith angles
		 *
		 * This method works for angles of arrival as well as angles 
		 * of departure, depending on which angle spread and angle 
		 * direction values are provided.
		 */
		VectorNd<double,4> recomputeZenithAngles(
				const VectorNd<bool,2>& LOSCondition,
				const VectorNd<double,2>& sigma_zs_LOS,
				const VectorNd<double,2>& sigma_zs_NLOS,
				const VectorNd<double,2>& sigma_kf,
				const VectorNd<double,3>& clusterPowers,
				const VectorNd<double,2>& angleDir,
				const bool arrival
				);

		/**
		 * @brief Generate random phases
		 */
		tuple<VectorNd<double,4>,VectorNd<double,2>>
		genRandomPhases(const VectorNd<bool,2>& LOSCondition);

		/**
		 * @brief Generate cross polarization values
		 */
		VectorNd<double,4> genCrossPolarization(VectorNd<bool,2>& LOSCondition);

		/**
		 * Precompute values for ray clusters
		 */
		VectorNd<RayCluster,5> precomputeRayValues(
				VectorNd<bool,2>& LOSCondition,
				const VectorNd<double,2>& sigma_kf,
				int numReceiverAntenna,
				int numSenderAntenna,
				const VectorNd<double,3>& clusterPowers,
				const VectorNd<double,4>& azimuth_ASA,
				const VectorNd<double,4>& azimuth_ASD,
				const VectorNd<double,4>& elevation_ASA,
				const VectorNd<double,4>& elevation_ASD,
				const VectorNd<array<double,3>,2>& receiverAntennaPos,
				const VectorNd<array<double,3>,2>& senderAntennaPos,
				const VectorNd<double,4>& randomPhase,
				const VectorNd<double,2>& randomPhase_LOS,
				const VectorNd<double,2>& AoA_LOS_dir,
				const VectorNd<double,2>& ZoA_LOS_dir,
				const VectorNd<double,2>& AoD_LOS_dir,
				const VectorNd<double,2>& ZoD_LOS_dir
				);

		/**
		 * @brief Precompute per-cluster delays
		 *
		 * For each receiver/sender pair, the METIS model generates 
		 * several rays with normal distributed path delay. See
		 * section 7.3.12 of the METIS D1.2 document.
		 *
		 * @param LOSCondition A [#receivers]x[#senders] vector 
		 * 			containing LOS conditions 
		 * 			(see genLosCond())
		 * @param sigmaDS_LOS The delay spread variance for line of sight 
		 * 			(see recomputeLargeScaleParameters() 
		 * 			and 3GPP Doc 25996)
		 * @param sigmaDS_NLOS The delay spread variance for non line 
		 * 			of sight (see recomputeLargeScaleParameters()
		 * 			and 3GPP Doc 25996)
		 * @param sigmaKF_LOS The variance of the Ricean K-factor for 
		 * 			line of sight (see recomputeLargeScaleParameters() 
		 * 			and 3GPP Doc 25996)
		 * @return A 3D vector [#receivers]x[#senders]x[#ray clusters] 
		 */
		VectorNd<double,3> precomputeClusterDelays(
				const VectorNd<bool,2>& LOSCondition,
				const VectorNd<double,2>& sigmaDS_LOS,
				const VectorNd<double,2>& sigmaDS_NLOS,
				const VectorNd<double,2>& sigmaKF_LOS
				);

		/**
		 * @brief Add subcluster delays for first two clusters
		 *
		 * See METIS D1.2 equation 7-60.
		 */
		VectorNd<std::complex<double>,4> addClusterDelayOffsets(
				VectorNd<double,3>& delays,
				bool up,
				size_t numRb);

		/**
		 * @brief Compute coefficients for given receivers/senders
		 */
		VectorNd<double,3> computeCoeffs(
				const VectorNd<bool,2>& LOSCondition,
				const vector<Position>& receiverPos,
				const vector<Position>& senderPos,
				double heightReceivers,
				double heightSenders,
				bool up,
				int numRBs,
				int numReceiverAntenna,
				int numSenderAntenna,
				const VectorNd<RayCluster,5>& rayClusters,
				const VectorNd<std::complex<double>,4>& delays
				);

		/**
		 * @brief Precompute the downlink (BS->MS) cluster/delay values
		 *
		 * After executing this method, the precomp and delay tables will be
		 * populated with all values which can be computed at the beginning of the
		 * simulation.
		 */
		void precomputeDownValues(const vector<Position>& msPositions,
				const vector<Position>& bsPositions);

		/**
		 * @brief Precompute the uplink (MS->BS) cluster/delay values
		 *
		 * After executing this method, the precomp and delay tables will be
		 * populated with all values which can be computed at the beginning of the
		 * simulation.
		 */
		void precomputeUpValues(const vector<vector<Position>>& msPositions,
				const vector<Position>& bsPositions);
		
		/**
		 * @brief Precompute the D2D (MS->MS) cluster/delay
		 *
		 * After executing this method, the precomp and delay tables will be
		 * populated with all values which can be computed at the beginning of the
		 * simulation.
		 */
		void precomputeD2DValues(const vector<vector<Position>>& msPositions);

	public:
		//! Initialize the METIS channel through ini access via OMNeT++ module pointer.
		bool init(cSimpleModule* module,
				const vector<vector<Position>>& msPositions, 
				std::map<int,Position>& neighbourPositions);

		//! Allows the OMNeT++ module to pass messages to this METIS channel class.
		void handleMessage(cMessage* msg);
		
		//! Updates the MS position if velocity > 0. The interval in which the postion is updated can be set within omnet.ini
		void updateChannel(const vector<vector<Position>>& msPos);

		//! Destructor of METIS Channel subclass.
		virtual ~METISChannel(){}
		/**
		 * @brief Compute per-TTI coefficients for the current time
		 *
		 * This method fills the coeff* table members with the coefficient values
		 * for the current TTI. This method makes use of the precomputed values.
		 *
		 * This means: The precomputeMETISValues method needs to be called before
		 * this method is called.
		 */
		void recomputePerTTIValues();
};
