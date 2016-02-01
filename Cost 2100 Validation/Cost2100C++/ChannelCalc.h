/*
 * ChannelCalc.h
 *
 *  Created on: Jul 9, 2013
 *      Author: Sascha Schmerling
 * Last edited: May 30, 2014
 *      Author: Thomas Prinz
 */

#pragma once

#include <assert.h>
#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include "Position.h"
#include "cluster.h"
#include "VisibilityRegion.h"

using namespace std;
using namespace itpp;

double dist(vec pos1, vec pos2);
vec Cart_to_Sph(vec input);
vec Sph_to_Cart(vec input);
mat rotate_matrix(const vec &phi_theta);

class ChannelCalc  {
    private:
        const static double speedOfLightVac;
        const static double carrierFreq;
        const static double bsAntennaHeight;
        const static double msAntennaHeight;
        
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
		double subcarrier_bandwidth;
		double excess_delay;
		double cluster_power;
        
        Position senderPosition;
        Position targetPosition;
        vector<Position> interfererPositions;

		double attenuation(double decay_factor, double tau_0, double tau_b, double tau);
		cmat calc_MIMO_Channel_Matrix( 	int num_rx_antennae, int num_tx_antennae, mat rx_positions, mat tx_positions,
										double azimuthOfArrival, double azimuthOfDeparture, double elevationOfArrival, double elevationOfDeparture,	 double lambda);
		mat compute_channel_matrix(double pathloss, Cluster c, int num_mpc, int num_dmc, vec bs_pos, vec ms_pos, double center_freq, double VR_gain);
		cmat polarization_matrix();
		double calcPathLossFS(Position a, Position b);
		double calc_path_loss_cost2100(const double dist);
		double calcPathLossCost231(Position a, Position b);
		double calc_VRGain(double lambda, double TR_Radius, double VR_Radius, vec VR_Center, vec MS_Pos);
        double distanceBetween(Position a, Position b);
        double calcSINR();

    public:
		void cost_init( double p_factor, double cutoff, double VR_los, double TR_los, double lambda, double vr_r, double tr_r,
						double tx, double rx, double startfreq, double endfreq, double subcarrier_band, double excess, double cluster_p);
		vector<cmat> H(Cluster local_bs, Cluster local_ms, vector<Cluster> &cluster,vector<VisibilityRegion> &vr, int num_mpc, int num_dmc, vec ms_pos, int bsId, std::map <int,Position> &neighbourPositions, VisibilityRegion &LOS_VR);
        void setSenderPosition(Position p) { senderPosition = p; }
        Position getSenderPosition() { return senderPosition; }
        void setTargetPosition(Position p) { targetPosition = p; }
        Position getTargetPosition() { return targetPosition; }
        void clearInterfererPostitions() { interfererPositions.clear(); }
        void addInterfererPosition(Position p) { interfererPositions.push_back(p); }
};
