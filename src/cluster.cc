/*
 * Cluster.cc
 * 
 * Datastructure for clusters
 *
 *  Created on: May 6, 2014
 *      Author: Thomas Prinz
 * 
 * Contains parts of the COST 2100 Implementation by Meifang Zhu (meifang.zhu@eit.lth.se)
 */
 
#include "cluster.h"

using namespace itpp;

/*
 * Apparently this function creates a rotation matrix, that does not rotate around
 * the X axes, rotates counterclockwise around the Y axes and clockwise (!)
 * around the z axes. The reason behind this is yet to be determined...
 * TODO: verify this actually does what we want...
 */
mat rotate_matrix(const vec &phi_theta) {
    double delta = 0;
    double phi = phi_theta.get(0);
    double theta = phi_theta.get(1);
    mat m = repmat(zeros(3), 1, 3);
    double sin_delta = sin(delta);
    double cos_delta = cos(delta);
    mat Tx = m;
    Tx.set(0, 0, 1);
    Tx.set(1, 1, cos_delta);
    Tx.set(1, 2, -1 * sin_delta);
    Tx.set(2, 1, sin_delta);
    Tx.set(2, 2, cos_delta);
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    mat Ty = m;
    Ty.set(0, 0, cos_theta);
    Ty.set(0, 2, sin_theta);
    Ty.set(1, 1, 1);
    Ty.set(2, 0, -1 * sin_theta);
    Ty.set(2, 2, cos_theta);
    double sin_phi = sin(phi);
    double cos_phi = cos(phi);
    mat Tz = m;
    Tz.set(0, 0, cos_phi);
    Tz.set(0, 1, sin_phi);
    Tz.set(1, 0, -1 * sin_phi);
    Tz.set(1, 1, cos_phi);
    Tz.set(2, 2, 1);
    return Tx * Ty * Tz;
}

// Serializing Datastructures
void doPacking(cCommBuffer *buffer, Cluster &c) {
    vec BS_Pos = c.get_cluster_pos(BS);
    vec MS_Pos = c.get_cluster_pos(MS);
    
    // Spread first
    buffer->pack(c.get_cluster_spread_vec(BS_MS)(0));
    buffer->pack(c.get_cluster_spread_vec(BS_MS)(1));
    buffer->pack(c.get_cluster_spread_vec(BS_MS)(2));
    buffer->pack(c.get_cluster_spread_vec(BS_MS)(3));
    buffer->pack(c.get_cluster_spread_vec(BS_MS)(4));
    buffer->pack(c.get_cluster_spread_vec(BS_MS)(5));
    
    // Positions second
    buffer->pack(BS_Pos(0));
    buffer->pack(BS_Pos(1));
    buffer->pack(BS_Pos(2));
    buffer->pack(MS_Pos(0));
    buffer->pack(MS_Pos(1));
    buffer->pack(MS_Pos(2));
    
    // Shadow Fading
    buffer->pack(c.get_cluster_shadowing_fading());
    
    // Link Delay
    buffer->pack(c.get_cluster_link_delay());
    
    // Angles
    buffer->pack(c.get_cluster_angle(AZIMUTH_AT_BS));
    buffer->pack(c.get_cluster_angle(ELEVATION_AT_BS));
    buffer->pack(c.get_cluster_angle(AZIMUTH_AT_MS));
    buffer->pack(c.get_cluster_angle(ELEVATION_AT_MS));
    
    // Type
    buffer->pack(c.get_cluster_type());
    buffer->pack(c.get_cluster_ID());
    
    // MPC
    mat bs_pos = c.get_cluster_MPC().get_MPC_pos(BS);
    mat ms_pos = c.get_cluster_MPC().get_MPC_pos(MS);
    
    // Mat Dimension is identical
    buffer->pack(bs_pos.rows());
    buffer->pack(bs_pos.cols());
    
    for(int i = 0; i < bs_pos.rows(); i++){
		for(int j = 0; j < bs_pos.cols(); j++){
			buffer->pack(bs_pos(i,j));
			buffer->pack(ms_pos(i,j));
		}
	}
}

// Deserializing Datastructures
void doUnpacking(cCommBuffer *buffer, Cluster &c) {
    vec spread, BS_Pos, MS_Pos, angles;
    double temp, sf, delay;
    CLUSTER_TYPE c_type;
    
    // Spread
    buffer->unpack(temp);
    spread.ins(0,temp);
    buffer->unpack(temp);
    spread.ins(1,temp);
    buffer->unpack(temp);
    spread.ins(2,temp);
    buffer->unpack(temp);
    spread.ins(3,temp);
    buffer->unpack(temp);
    spread.ins(4,temp);
    buffer->unpack(temp);
    spread.ins(5,temp);
    
    // Positions
    buffer->unpack(temp);
    BS_Pos.ins(0,temp);
    buffer->unpack(temp);
    BS_Pos.ins(1,temp);
    buffer->unpack(temp);
    BS_Pos.ins(2,temp);
    buffer->unpack(temp);
    MS_Pos.ins(0,temp);
    buffer->unpack(temp);
    MS_Pos.ins(1,temp);
    buffer->unpack(temp);
    MS_Pos.ins(2,temp);
    
    // Shadow Fading
    buffer->unpack(sf);
    
    // Link Delay
    buffer->unpack(delay);
    
    // Angles
    buffer->unpack(temp);
    angles.ins(0,temp);
    buffer->unpack(temp);
    angles.ins(1,temp);
    buffer->unpack(temp);
    angles.ins(2,temp);
    buffer->unpack(temp);
    angles.ins(3,temp);
    
    // Type
    int temp_type;
    buffer->unpack(temp_type);
    c_type = static_cast<CLUSTER_TYPE>(temp_type);
    
    int id;
    buffer->unpack(id);
    
    int a,b;
    
    buffer->unpack(a);
    buffer->unpack(b); 
    
    mat bs_pos(a,b);
    mat ms_pos(a,b);
    vec delay_vec = zeros(5);
    double entry;
    
    for(int i = 0; i < a; i++){
		for(int j = 0; j < b; j++){
			buffer->unpack(entry);
			bs_pos.set(i,j,entry);
			buffer->unpack(entry);
			ms_pos.set(i,j,entry);
		}
	}
	
	MPC mpc = MPC(bs_pos, ms_pos, 5, zeros(5));
    
	c = Cluster(c_type, BS_Pos, MS_Pos, spread, sf, delay, angles, id, &mpc);
}

Cluster::Cluster(const CLUSTER_TYPE c_type, const vec &c_bs_pos, const vec &c_ms_pos, const vec &spread, double c_sf, double c_tau_link, const vec &angle, int id, MPC *oldmpc) {
    C_type = c_type;
	C_BS_pos = c_bs_pos;
	C_MS_pos = c_ms_pos;
	C_delay_BS = spread.get(0);
	C_AoD_BS = spread.get(1);
	C_EoD_BS = spread.get(2);
	C_delay_MS = spread.get(3);
	C_AoD_MS = spread.get(4);
	C_EoD_MS = spread.get(5);
	C_shadow_f = c_sf;
	C_tau_link = c_tau_link;
	C_Phi_BS = angle.get(0);
	C_Theta_BS = angle.get(1);
	C_Phi_MS = angle.get(2);
	C_Theta_MS = angle.get(3);
	clusterId = id;
	
	/*
    std::cout << "Generated Cluster: " << c_type << std::endl;
    std::cout << "Spatial Spread of the ellipsoid at BS (x,y,z): " << spread.get(0,2) << std::endl;
    std::cout << "Spatial Spread of the ellipsoid at MS (x,y,z): " << spread.get(3,5) << std::endl;
    std::cout << "Spatial Spread Both (x,y,z): " << spread << std::endl;
    std::cout << "Shadow Fading value: " << c_sf << std::endl;
    std::cout << "Link Delay (0 for all cluster except TWIN): " << c_tau_link << "s" << std::endl;
    std::cout << "Angle at BS: " << angle.get(0,1) << std::endl;
    std::cout << "Angle as MS: " << angle.get(2,3) << std::endl;
    */
    
	int N_MPC = 5;
	vec bs_angle = angle.get(0,1);
	vec ms_angle = angle.get(2,3);
	
	mat MPC_BS_pos, MPC_MS_pos;
	mat tmp;
	mat c_bs_pos_mat = repmat(((mat) c_bs_pos).transpose(),N_MPC,1);
	mat c_ms_pos_mat = repmat(((mat) c_ms_pos).transpose(),N_MPC,1);
    
    // TODO: get from config.
    if(oldmpc == NULL){	
		// Multi Path Components
		// Uniform for local cluster, normal distributed for single/twin
		switch (c_type) {
			case LOCAL_AT_BS: // local cluster at BS
				tmp = randu(N_MPC, 3) * diag(spread.get(0,2) / 3); 
				MPC_BS_pos = tmp + c_bs_pos_mat;
				MPC_MS_pos = MPC_BS_pos;
				break;
			case SINGLE: // single cluster
				tmp = randn(N_MPC, 3) * diag(spread.get(0,2) / 3);
				MPC_BS_pos = tmp * rotate_matrix(bs_angle) + c_bs_pos_mat;
				MPC_MS_pos = MPC_BS_pos;
				break;
			case TWIN: // twin cluster
				tmp = randn(N_MPC, 3) * diag(spread.get(0,2) / 3);
				MPC_BS_pos = tmp * rotate_matrix(bs_angle) + c_bs_pos_mat;
				MPC_MS_pos = tmp * rotate_matrix(ms_angle) + c_ms_pos_mat;
				break;
			case LOCAL_AT_MS: // local cluster at MS
				tmp = randu(N_MPC, 3) * diag(spread.get(0,2) / 3);
				MPC_MS_pos = tmp + c_ms_pos_mat;
				MPC_BS_pos = MPC_MS_pos;
				break;
		}
		mpc = MPC(MPC_BS_pos, MPC_MS_pos, N_MPC, zeros(N_MPC));
	}else{
		mpc = *oldmpc;
	}
	
	// Dense/Diffuse Multipaths
	// Geometric Approach as proposed within the book "Adding Dense Multipath
	// Components to Geometry-Based MIMO Channel Models"
	// The radius parameter has to be adjusted based on measurements
	// The delay domain is Poisson distributed
	int num_dmc = 15;
	
	double radius_dmc = sqrt(pow(spread(0),2) + pow(spread(1),2) + pow(spread(2),2));
	double add_delay_mean = 0.0;
	vec dmc_spread;
	dmc_spread.ins(0,radius_dmc);
	dmc_spread.ins(0,radius_dmc);
	dmc_spread.ins(0,radius_dmc);
	
	c_bs_pos_mat = repmat(((mat) c_bs_pos).transpose(),num_dmc,1);
    c_ms_pos_mat = repmat(((mat) c_ms_pos).transpose(),num_dmc,1);
	
	switch(c_type){
    	case LOCAL_AT_BS: // local cluster at BS
			tmp = randu(num_dmc, 3) * diag(dmc_spread);
    		MPC_BS_pos = tmp + c_bs_pos_mat;
    		MPC_MS_pos = MPC_BS_pos;
    		break;
    	case SINGLE: // single cluster
			tmp = randn(num_dmc, 3) * diag(dmc_spread);
    		MPC_BS_pos = tmp * rotate_matrix(bs_angle) + c_bs_pos_mat;
    		MPC_MS_pos = MPC_BS_pos;
    		break;
    	case TWIN: // twin cluster
			tmp = randn(num_dmc, 3) * diag(dmc_spread);
	   		MPC_BS_pos = tmp * rotate_matrix(bs_angle) + c_bs_pos_mat;
	   		MPC_MS_pos = tmp * rotate_matrix(ms_angle) + c_ms_pos_mat;
	   		break;
    	case LOCAL_AT_MS: // local cluster at MS
			tmp = randu(num_dmc, 3) * diag(dmc_spread);
	   		MPC_MS_pos = tmp + c_ms_pos_mat;
	   		MPC_BS_pos = MPC_MS_pos;
	   		break;
	}
	
	// TODO: correct additional delay
	vec additional_delay = zeros(num_dmc);
	for(int i = 0; i < num_dmc; i++){
		additional_delay.ins(i,poisson(add_delay_mean));
	}
	
	dmc = MPC(MPC_BS_pos, MPC_MS_pos, num_dmc, additional_delay);
}

double Cluster::get_cluster_spread(const SPREAD_TYPE s_type) {
	double spread;
	switch (s_type) {
		case DELAY_AT_BS:
			spread = C_delay_BS;
			break;

		case DELAY_AT_MS:
			spread = C_delay_MS;
			break;

		case AOD_AT_BS:
			spread = C_AoD_BS;
			break;

		case AOD_AT_MS:
			spread = C_AoD_MS;

		case EOD_AT_BS:
			spread = C_EoD_BS;
			break;

		case EOD_AT_MS:
			spread = C_EoD_MS;
			break;

	}
	return spread;
}

vec Cluster::get_cluster_spread_vec(const SIDE_TYPE sd_type) {
	vec s_vec;
	switch (sd_type) {
		case BS:
			s_vec.ins(0, C_delay_BS);
			s_vec.ins(1, C_AoD_BS);
			s_vec.ins(2, C_EoD_BS);
			break;
		case MS:
			s_vec.ins(0, C_delay_MS);
			s_vec.ins(1, C_AoD_MS);
			s_vec.ins(2, C_EoD_MS);
			break;
		case BS_MS:
			s_vec.ins(0, C_delay_BS);
			s_vec.ins(1, C_AoD_BS);
			s_vec.ins(2, C_EoD_BS);
			s_vec.ins(3, C_delay_MS);
			s_vec.ins(4, C_AoD_MS);
			s_vec.ins(5, C_EoD_MS);
			break;
	}
	return s_vec;
};

double Cluster::get_cluster_angle(const ANGLE_TYPE a_type) {
	double angle;
	switch (a_type) {
		case AZIMUTH_AT_BS:
			angle = C_Phi_BS;
			break;

		case AZIMUTH_AT_MS:
			angle = C_Phi_MS;
			break;

		case ELEVATION_AT_BS:
			angle = C_Theta_BS;
			break;

		case ELEVATION_AT_MS:
			angle = C_Theta_MS;
			break;

	}
	return angle;
}

vec Cluster::get_cluster_angle_vec(const SIDE_TYPE sd_type) {
	vec a_vec;
	switch (sd_type) {
		case BS:
			a_vec.ins(0, C_Phi_BS);
			a_vec.ins(1, C_Theta_BS);
			break;

		case MS:
			a_vec.ins(0, C_Phi_MS);
			a_vec.ins(1, C_Theta_MS);
			break;

		case BS_MS:
			a_vec.ins(0, C_Phi_BS);
			a_vec.ins(1, C_Theta_BS);
			a_vec.ins(2, C_Phi_MS);
			a_vec.ins(3, C_Theta_MS);
			break;

	}
	return a_vec;
}
