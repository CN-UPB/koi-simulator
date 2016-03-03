/*
 * Cluster.h
 * 
 * Datastructure for clusters
 *
 *  Created on: May 5, 2014
 *      Author: Thomas Prinz
 * 
 * Contains parts of the COST 2100 Implementation by Meifang Zhu (meifang.zhu@eit.lth.se)
 */

#ifndef CLUSTER_H
#define CLUSTER_H

#include "includes.h"
#include <itpp/itbase.h>
#include <vector>
#include "mpc.h"

using namespace itpp;

mat rotate_matrix(const vec &phi_theta); 

// Cluster type: Single Cluster, Twin Cluster, Local BS Cluster or Local MS Cluster
enum CLUSTER_TYPE {
	SINGLE, TWIN, LOCAL_AT_BS, LOCAL_AT_MS
};

// Spread type: Delay Spread, AoD Angular Spread or EoD Angular Spread to BS/MS side
enum SPREAD_TYPE {
	DELAY_AT_BS, DELAY_AT_MS, AOD_AT_BS, AOD_AT_MS, EOD_AT_BS, EOD_AT_MS
};
// Angle type: Azimuth or Elevation Angle at BS/MS side
enum ANGLE_TYPE {
	AZIMUTH_AT_BS, AZIMUTH_AT_MS, ELEVATION_AT_BS, ELEVATION_AT_MS
};

class Cluster{
public:
	//! Empty constructor
	Cluster() {}
	//! Default constructor
	Cluster(const CLUSTER_TYPE c_type, const vec &c_bs_pos, const vec &c_ms_pos, const vec &spread, double c_sf, double c_tau_link, const vec &angle, int id, MPC *oldmpc = NULL);
	//! Destructor
	virtual ~Cluster() {}

	//! Get cluster type
	CLUSTER_TYPE get_cluster_type() { return C_type; };
	//! Get cluster positions to BS/MS side
	vec get_cluster_pos(const SIDE_TYPE sd_type) { switch (sd_type) { case BS: return C_BS_pos;	case MS: return C_MS_pos; default: return zeros(3); } }
	//! Get cluster shadow fading
	double get_cluster_shadowing_fading() { return C_shadow_f; }
	//! Get cluster link delay
	double get_cluster_link_delay() { return C_tau_link; }
	//! Get cluster spread
	double get_cluster_spread(const SPREAD_TYPE s_type);
	//! Get cluster spread vector to BS/MS side
	vec get_cluster_spread_vec(const SIDE_TYPE sd_type);
	//! Get cluster angle
	double get_cluster_angle(const ANGLE_TYPE a_type);
	//! Get Cluster ID
	int get_cluster_ID(){return clusterId;}
	//! Get cluster angle vector to BS/MS side
	vec get_cluster_angle_vec(const SIDE_TYPE sd_type);
	//! Get cluster MPCs
	MPC get_cluster_MPC() { return mpc; }
	//! Get cluster DMCs
	MPC get_cluster_DMC() { return dmc; }

protected:
	CLUSTER_TYPE C_type;					//!< cluster type
	vec C_BS_pos;							//!< cluster position to BS side
	vec C_MS_pos;							//!< cluster position to MS side
	double C_delay_BS;						//!< cluster spatial delay spread to BS side
	double C_delay_MS;						//!< cluster spatial delay spread to MS side
	double C_AoD_BS;						//!< cluster AoD to BS side
	double C_AoD_MS;						//!< cluster AoD to MS side
	double C_EoD_BS;						//!< cluster EoD to BS side
	double C_EoD_MS;						//!< cluster EoD to MS side
	double C_shadow_f;						//!< cluster shadowing fading.
	double C_tau_link;						//!< cluster link delay.
	double C_Phi_BS;						//!< cluster azimuth angle to BS side
	double C_Theta_BS;						//!< cluster elevation angle to BS side
	double C_Phi_MS;						//!< cluster azimuth angle to MS side
	double C_Theta_MS;						//!< cluster elevation angle to MS side
	int clusterId;
    MPC mpc;
    MPC dmc;
};

void doPacking(cCommBuffer *buffer, Cluster &c);

void doUnpacking(cCommBuffer *buffer, Cluster &c);


#endif // CLUSTER_H
