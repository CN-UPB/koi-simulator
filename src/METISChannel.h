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
#include "Channel.h"
#include <algorithm>
#include <vector>
#include <tuple>
#include <array>
#include <ostream>

using std::vector;
using std::array;
using std::tuple;

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
		int timeSamples;					/*!< Number of TTIs until Position is updated (Number of Time Samples for Channel Model) */
		vector<vector<vector<vector<double>>>> coeffDownTable;				/*!< Table to save downlink coefficients */
		vector<vector<vector<vector<vector<double>>>>> coeffUpTable;				/*!< Table to save uplink coefficients */
		vector<vector<vector<vector<vector<double>>>>> coeffDownD2DTable;				/*!< Table to save D2D DOWN Rb coefficients */
		vector<vector<vector<vector<vector<double>>>>> coeffUpD2DTable;				/*!< Table to save D2D UP Rb coefficients */
		int SINRcounter;					/*!< If position resend intervall > 1, it counts the current TTI */
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
				vector<vector<double>>& sigma_sf_NLOS);

		//! Recalculate all position dependent values, e.g. SINR
		void recomputeMETISParams(const vector<vector<Position>>& msPositions);

		//! Generate the spatial correlation between the MS for LOS links.
		void generateAutoCorrelation_LOS(const vector<Position>& senders,
				const vector<Position>& receivers,
				vector<vector<vector<double>>>& correlation);

		//! Generate the spatial correlation between the MS for NLOS links.
		void generateAutoCorrelation_NLOS(const vector<Position>& senders,
				const vector<Position>& receivers,
				vector<vector<vector<double>>>& correlation);

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
		vector<vector<bool>> genLosCond(const vector<Position>& sendPos,
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
		std::tuple<vector<vector<vector<double>>>,vector<vector<vector<double>>>>
		recomputeClusterDelays(const vector<vector<bool>>& LOSCondition,
				const vector<vector<double>>& sigmaDS_LOS,
				const vector<vector<double>>& sigmaDS_NLOS,
				const vector<vector<double>>& sigmaKF_LOS);
        
		/**
		 * @brief Recompute per-cluster powers
		 *
		 */
		vector<vector<vector<double>>> genClusterPowers(const vector<vector<bool>>& LOSCondition,
				const vector<vector<vector<double>>>& clusterDelays,
				const vector<vector<double>>& sigmaDS_LOS,
				const vector<vector<double>>& sigmaDS_NLOS,
				const vector<vector<double>>& sigmaKF_LOS
				);

		/**
		 * @brief Recompute per-ray powers
		 */
		vector<vector<vector<double>>> recomputeRayPowers(const vector<vector<bool>>& LOSCondition,
				vector<vector<vector<double>>>& clusterPowers
				);

		/**
		 * @brief Recompute angle directions
		 */
		tuple<vector<vector<double>>,vector<vector<double>>,vector<vector<double>>,vector<vector<double>>> 
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
		vector<vector<vector<vector<double>>>>
		recomputeAzimuthAngles(const vector<vector<bool>>& LOSCondition,
				const vector<vector<double>>& sigma_as_LOS,
				const vector<vector<double>>& sigma_as_NLOS,
				const vector<vector<double>>& sigma_kf,
				const vector<vector<vector<double>>>& clusterPowers,
				const vector<vector<double>>& angleDir,
				const bool arrival
				);

		/**
		 * @brief Recompute ray zenith angles
		 *
		 * This method works for angles of arrival as well as angles 
		 * of departure, depending on which angle spread and angle 
		 * direction values are provided.
		 */
		vector<vector<vector<vector<double>>>> recomputeZenithAngles(
				const vector<vector<bool>>& LOSCondition,
				const vector<vector<double>>& sigma_zs_LOS,
				const vector<vector<double>>& sigma_zs_NLOS,
				const vector<vector<double>>& sigma_kf,
				const vector<vector<vector<double>>>& clusterPowers,
				const vector<vector<double>>& angleDir,
				const bool arrival
				);

		/**
		 * @brief Generate random phases
		 */
		tuple<vector<vector<vector<vector<vector<double>>>>>,vector<vector<double>>>
		genRandomPhases( const vector<vector<bool>>& LOSCondition);

		/**
		 * @brief Generate cross polarization values
		 */
		vector<vector<vector<vector<double>>>> genCrossPolarization(
				vector<vector<bool>>& LOSCondition);

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
				vector<vector<vector<std::complex<double>>>>& raySum
				);

		/**
		 * @brief Compute ray sums for given receivers/senders
		 */
		tuple<vector<vector<vector<vector<vector<vector<std::complex<double>>>>>>>,
			vector<vector<vector<vector<vector<vector<std::complex<double>>>>>>>>
				computeRaySums(vector<vector<bool>>& LOSCondition,
						const vector<vector<double>>& sigma_kf,
						int numReceiverAntenna,
						int numSenderAntenna,
						const vector<vector<vector<double>>>& clusterPowers,
						const vector<vector<vector<vector<double>>>>& azimuth_ASA,
						const vector<vector<vector<vector<double>>>>& azimuth_ASD,
						const vector<vector<vector<vector<double>>>>& elevation_ASA,
						const vector<vector<vector<vector<double>>>>& elevation_ASD,
						const vector<vector<array<double,3>>>& receiverAntennaPos,
						const vector<vector<array<double,3>>>& senderAntennaPos,
						const vector<vector<vector<vector<vector<double>>>>>& randomPhase,
						const vector<vector<double>>& randomPhase_LOS,
						const vector<vector<double>>& AoA_LOS_dir,
						const vector<vector<double>>& ZoA_LOS_dir,
						const vector<vector<double>>& AoD_LOS_dir,
						const vector<vector<double>>& ZoD_LOS_dir
						);

		/**
		 * @brief Compute coefficients for given receivers/senders
		 */
		vector<vector<vector<vector<double>>>> computeCoeffs(
				const vector<vector<bool>>& LOSCondition,
				const vector<Position>& receiverPos,
				const vector<Position>& senderPos,
				double heightReceivers,
				double heightSenders,
				bool up,
				int numRBs,
				int numReceiverAntenna,
				int numSenderAntenna,
				const vector<vector<vector<vector<vector<vector<std::complex<double>>>>>>>& raySum,
				const vector<vector<vector<vector<vector<vector<std::complex<double>>>>>>>& raySum_LOS,
				const vector<vector<vector<double>>>& clusterDelays,
				const vector<vector<vector<double>>>& clusterDelays_LOS
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

		/**
		 * @brief Output the Up coefficient table to out stream
		 */
		std::ostream& printCoeffUpTables(std::ostream& out);

		/**
		 * @brief Output the Down coefficient table to out stream
		 */
		std::ostream& printCoeffDownTables(std::ostream& out);

		/**
		 * @brief Calculate interference for transmission
		 */
		double calcInterference(std::forward_list<TransInfo*>& interferers,
				int rb,
				int receiverId,
				int SINRCounter,
				MessageDirection dir);

	public:
		//! Initialize the METIS channel through ini access via OMNeT++ module pointer.
		bool init(cSimpleModule* module,
				const vector<vector<Position>>& msPositions, 
				std::map<int,Position>& neighbourPositions);

		//! Allows the OMNeT++ module to pass messages to this METIS channel class.
		void handleMessage(cMessage* msg);
		
		//! Computes the pathloss for a given distance using an arbitrary model.
		double calcPathloss(double dist);
		
		//! Computes the Termal Noise ("johnson nyquist noise")
		double getTermalNoise(double temp, double bandwidth);
		
		double calcUpSINR(int RB, 
				int msId,
				double transPower);

		double calcDownSINR(int RB, 
				int msId,
				double transPower);

		double calcD2DSINR(int RB, 
				int sendMsID,
				int receiveMsId,
				MessageDirection direction,
				double transPower);

		double calcAvgUpSINR(int RB, 
				int msId,
				double transPower);

		double calcAvgDownSINR(int RB, 
				double transPower);

		double calcAvgD2DDownSINR(int RB, 
				int msId,
				double transPower);

		
		//! Updates the MS position if velocity > 0. The interval in which the postion is updated can be set within omnet.ini
		void updateChannel(const vector<vector<Position>>& msPos);

		//! Destructor of METIS Channel subclass.
		virtual ~METISChannel();
};
