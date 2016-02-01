/*
 * Cost2100Channel.h
 *
 *  Created on: Jul 15, 2014
 *      Author: Thomas Prinz
 * 
 */
 
#pragma once

#include <omnetpp.h>
#include <unordered_map>
#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>
#include "Position.h"
#include "NeighbourIdMatching.h"
#include "VisibilityRegion.h"
#include "cluster.h"
#include "VisibilityRegionMessage_m.h"
#include "ClusterMessage_m.h"
#include "Channel.h"
#include "PointerExchange_m.h"

using namespace std;
using namespace itpp;

class Cost2100Channel : public Channel{
	private:
	    // Cost 2100 parameters
        double power_factor;
        double LOS_cutoff_distance;
        double VR_radius_los;
        double TR_radius_los;
        double lambda;
        double vrRadius;
        double trRadius;
        int tx_num;
        int rx_num;
        double startfrequency;
		double endfrequency;
		double carrierFreq;
		double subcarrier_bandwidth;
		double excess_delay;
		double cluster_power;
		int bsId;
		int maxNumberOfNeighbours;
		int numberOfNeighbours;
		int init_counter;
		int numberOfVR;
		simtime_t tti;
		map <int,Position> neighbourPos;
        int numberOfMobileStations;
		
		// Save the init module for later access
		cSimpleModule *initModule;
		
		// Cost 2100 Environment
		// The Visibility Region belonging to this BS
        //vector<VisibilityRegion> VR;
        std::map<int,vector<VisibilityRegion>> VR;
        // The Line of sight Visibility Region
        VisibilityRegion LOS_VR;
        // The local BS cluster
        Cluster localcluster_bs;
        // The local clusters for each MS
        vector<Cluster> localclusters_ms;
        // All Single and Twin Remote Clusters
        unordered_map <int, Cluster> remoteCluster;
        // We need the VR from other BSs to create the remote cluster
        vector<VisibilityRegion> ForeignVR;
        
        double attenuation(double decay_factor, double tau_0, double tau_b, double tau);
		void calc_MIMO_Channel_Matrix( 	int num_rx_antennae, int num_tx_antennae, mat const &rx_positions, mat const &tx_positions,
										double azimuthOfArrival, double azimuthOfDeparture, double elevationOfArrival, double elevationOfDeparture,	 double lambda, cmat &result);
		mat compute_channel_matrix(double pathloss, Cluster &c, int num_mpc, int num_dmc, vec const &bs_pos, vec const &ms_pos, double center_freq, double VR_gain);
		cmat polarization_matrix();
		double calc_VRGain(double lambda, double TR_Radius, double VR_Radius, vec VR_Center, vec MS_Pos);
		vec H(Cluster &local_bs, vector<Cluster> &local_ms, unordered_map <int, Cluster> &cluster,std::map<int,vector<VisibilityRegion>> &vr, int num_mpc, int num_dmc, int bsId, Position &ms_pos, Position &bs_Pos, VisibilityRegion &LOS_VR, int RB = -1);

	public:
		// Constructor
		Cost2100Channel(){bsId = -1;}
		// Initialize your Channel through ini access via module pointer.
		bool init(cSimpleModule* module, Position** msPositions, std::map <int,Position> neighbourPositions);
		// It may be necessary for the Channel to receive Message from other LPs.
		void handleMessage(cMessage* msg);
		// Computes the pathloss for a given distance using an arbitrary model.
		double calcPathloss(double dist);
		// Computes the Termal Noise.
		double getTermalNoise(double temp, double bandwidth);
		// Calculates the current SINR for given interferers and given RB.
		double calcSINR(int RB, vector<double> &power, vector<Position> &pos, vector<int> &bsId_, bool up, int msId);
		// Calculates the current SINR for given interferers and given RB.
		vec calcSINR(vector<double> &power, vector<Position> &pos, vector<int> &bsId_, bool up, int msId);
		// Updates the Channel if necessary for moving MS
		void updateChannel(Position** msPos);
		//Destructor
		~Cost2100Channel();
};
