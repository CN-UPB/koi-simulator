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
#include <cmath>
#include "Position.h"
#include "Channel.h"
#include "VecNd.h"

#include <algorithm>
#include <array>
#include <complex>
#include <ostream>
#include <tuple>
#include <unordered_map>
#include <vector>

using std::vector;
using std::array;
using std::tuple;

class Ray{
	private:
		double dirAoA;
		std::complex<double> expArrival;
		std::complex<double> expDeparture;
		std::complex<double> pol;
	
	public:
		Ray(const double dirAoA, const std::complex<double> expArrival,
				const std::complex<double> expDeparture, const std::complex<double> pol)
				:dirAoA(dirAoA),expArrival(expArrival),
				expDeparture(expDeparture),pol(pol){}

		static Ray initialize(
				double azimuthASA,
				double azimuthASD,
				double zenithASA,
				double zenithASD,
				double k_0,
				const array<double,3>& senderAntennaPos,
				const array<double,3>& receiverAntennaPos,
				const vector<double>& randomPhase
				);
		virtual std::complex<double> value(double t, double moveAngle,
				double velocity, double k_0);

};

class LOSRay: public Ray{
	public:
		LOSRay(const double dirAoA, const std::complex<double> expArrival,
				const std::complex<double> expDeparture, const std::complex<double> pol)
				:Ray(dirAoA,expArrival,expDeparture,pol){}
		static LOSRay initialize(
				double dirAoA,
				double dirAoD,
				double dirZoA,
				double dirZoD,
				double k_0,
				const array<double,3>& senderAntennaPos,
				const array<double,3>& receiverAntennaPos,
				double randomPhase
				);
};

class METISChannel : public Channel{
	private:
		// METIS Channel Parameters
		double freq_c;						/*!< center/carrier frequence */
		double heightUE;					/*!< Height of the user equipments */
		double heightBS;					/*!< Height of the base stations */
		static double ray_offset[20];				/* Ray offset. Table 7.6 METIS D1.2 */
		int N_cluster_LOS;
		int N_cluster_NLOS;
		int numOfRays_LOS;
		int numOfRays_NLOS;
		int NumBsAntenna;					/*!< Number of Base Station Antenna */
		int NumMsAntenna;					/*!< Number of Mobile Station Antenna */
		vector<vector<array<double,3>>> bsAntennaPositions;				/*!< Position vector of Base Station antenna */
		int numOfInterferers;					/*!< Number of actual interferers, based on network layout and neighbour distance */
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

		//! Recalculate all position dependent values, e.g. SINR
		void recomputeMETISParams(const vector<vector<Position>>& msPositions);

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
		 * @brief Recompute per-cluster delays
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
		 * @return A pair of [#receivers]x[#senders]x[#ray clusters] vectors 
		 * 		where the first component contains the spreads 
		 * 		for the LOS case and second one containing those 
		 * 		for the NLOS case.
		 */
		std::tuple<VectorNd<double,3>,VectorNd<double,3>>
		recomputeClusterDelays(const VectorNd<bool,2>& LOSCondition,
				const VectorNd<double,2>& sigmaDS_LOS,
				const VectorNd<double,2>& sigmaDS_NLOS,
				const VectorNd<double,2>& sigmaKF_LOS);
        
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
		tuple<VectorNd<double,5>,VectorNd<double,2>>
		genRandomPhases(const VectorNd<bool,2>& LOSCondition);

		/**
		 * @brief Generate cross polarization values
		 */
		VectorNd<double,4> genCrossPolarization(VectorNd<bool,2>& LOSCondition);

		/**
		 * @brief Compute ray sum for a full cluster
		 */
		void computeRaySumCluster(
				size_t numRays,
				double prefactor,
				double k_0,
				const vector<double>& zenithASA,
				const vector<double>& zenithASD,
				const vector<double>& azimuthASA,
				const vector<double>& azimuthASD,
				const array<double,3>& senderAntennaPos,
				const array<double,3>& receiverAntennaPos,
				size_t receiverAntennaIndex,
				size_t senderAntennaIndex,
				const vector<vector<double>>& randomPhase,
				vector<int> *subcluster,
				VectorNd<std::complex<double>,2>& raySum
				);

		/**
		 * @brief Compute ray sums for given receivers/senders
		 */
		tuple<VectorNd<std::complex<double>,5>, VectorNd<std::complex<double>,5>>
				computeRaySums(VectorNd<bool,2>& LOSCondition,
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
						const VectorNd<double,5>& randomPhase,
						const VectorNd<double,2>& randomPhase_LOS,
						const VectorNd<double,2>& AoA_LOS_dir,
						const VectorNd<double,2>& ZoA_LOS_dir,
						const VectorNd<double,2>& AoD_LOS_dir,
						const VectorNd<double,2>& ZoD_LOS_dir
						);

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
				const VectorNd<std::complex<double>,5>& raySum,
				const VectorNd<std::complex<double>,5>& raySum_LOS,
				const VectorNd<double,3>& clusterDelays,
				const VectorNd<double,3>& clusterDelays_LOS
				);

		/**
		 * @brief Compute the downlink (BS->MS) coefficients
		 *
		 * After executing this method, the coeffDown table will hold 
		 * the coefficients for the links from all BS to the 
		 * local MS.
		 */
		void recomputeDownCoefficients(const vector<Position>& msPositions,
				const vector<Position>& bsPositions);

		/**
		 * @brief Compute the uplink (MS->BS) coefficients
		 *
		 * After executing this method, the coeffUp table will hold 
		 * the coefficients for the links from all MS to the local BS
		 */
		void recomputeUpCoefficients(const vector<vector<Position>>& msPositions,
				const vector<Position>& bsPositions);
		
		/**
		 * @brief Compute the D2D (MS->MS) coefficients
		 *
		 * After executing this method, the coeffUpD2D and coeffDownD2D tables will hold 
		 * the coefficients for the links from all MS to all local MS
		 */
		void recomputeD2DCoefficients(const vector<vector<Position>>& msPositions);

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
};
