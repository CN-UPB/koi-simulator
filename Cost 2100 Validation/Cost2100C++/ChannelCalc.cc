/*
 * ChannelCalc.h
 *
 *  Created on: Jul 9, 2013
 *      Author: Sascha Schmerling
 * Last edited: Jun 25, 2014
 *      Author: Thomas Prinz
 * 
 * All book references target the book "Pervasive Mobile and Ambient Wireless Communications"
 * 
 * Contains parts of the COST 2100 Implementation by Meifang Zhu (meifang.zhu@eit.lth.se)
 */

#include "ChannelCalc.h"

using namespace itpp;

const double ChannelCalc::speedOfLightVac = 3e8; // for validating reasons
const double ChannelCalc::carrierFreq = 1.8e9;
const double ChannelCalc::bsAntennaHeight = 50;
const double ChannelCalc::msAntennaHeight = 3;

/*
 * Cartesian: (x,y,z)
 * Spherical (Theta,Phi,r) [azimuth,elevation,r] (MatLab 'Convention')
 * The formula is identical to the Matlab intern one.
 */
vec Cart_to_Sph(vec input){
	vec output = zeros(3);
	output.set(0,atan2(input(1), input(0)));
	output.set(1,atan2(input(2), sqrt(pow(input(0),2) + pow(input(1),2))));
	output.set(2,sqrt(pow(input(0),2) + pow(input(1),2) + pow(input(2),2)));
	return output;
}

/*
 * Cartesian: (x,y,z)
 * Spherical (Theta,Phi,r) [azimuth,elevation,r] (MatLab 'Convention')
 * The formula is identical to the Matlab intern one.
 */
vec Sph_to_Cart(vec input){
	vec output = zeros(3);
	output.set(0,input(2) * cos(input(1)) * cos(input(0)));
	output.set(1,input(2) * cos(input(1)) * sin(input(0)));
	output.set(2,input(2) * sin(input(1)));
	return output;
}

/*
 * Initializes Variables for Cost2100 Computations
 */
void ChannelCalc::cost_init( double p_factor, double cutoff, double VR_los, double TR_los, double l, double vr_r, double tr_r,
							 double tx, double rx, double startfreq, double endfreq, double subcarrier_band, double excess, double cluster_p){
	power_factor = p_factor;
	LOS_cutoff_distance = cutoff;
	VR_radius_los = VR_los;
	TR_radius_los = TR_los;
	lambda = l;
	vrRadius = vr_r;
	trRadius = tr_r;
	tx_num = tx;
	rx_num = rx;
	startfrequency = startfreq;
	endfrequency = endfreq;
	subcarrier_bandwidth = subcarrier_band;			
	cluster_power = cluster_p;
	excess_delay = excess;		
}

/*
 * Distance between two points in 3D Space
 */
double dist(vec pos1, vec pos2){
	return sqrt(pow(pos1(0) - pos2(0),2) + pow(pos1(1) - pos2(1),2) + pow(pos1(2) - pos2(2),2));
}

/*
 * Computes attenuation according to Formula 3.40
 * tau_0 = Line of Sight Delay
 * tau_b = Cut off Delay
 * tau = Actual delay
 * *1e6 because parameter in 'us' and formula in 's'.
 */
double ChannelCalc::attenuation(double decay_factor, double tau_0, double tau_b, double tau){
	return std::max( exp(-1.0*decay_factor*(tau - tau_0)*1e6) , exp(-1.0*decay_factor*(tau_b - tau_0)*1e6) );
}

/*
 * Computes the Final Matrix over all active clusters.
 */
vector<cmat> ChannelCalc::H(Cluster local_bs, Cluster local_ms, vector<Cluster> &cluster,vector<VisibilityRegion> &vr, int num_mpc, int num_dmc, vec ms_pos, int bsId, std::map <int,Position> &neighbourPositions, VisibilityRegion &LOS_VR){
	std::cout << "Position of Mobilestation: " << ms_pos << " (lies within cell " << bsId << ")" << std::endl;
	
	vector<cmat> H;
	//Position bsPos = neighbourPositions[bsId];
	
	// BS Poshardcoded to (0/0)
	vec bsPosition;
	bsPosition.ins(0,0);
	bsPosition.ins(1,0);
	bsPosition.ins(2,0);
	
	std::cout << "Position of Basestation: " << bsPosition << std::endl;
	
	//! BEGIN Generate Signal Matrix
	
	vector<VisibilityRegion> signalVRs;
	vector<Cluster> signalCluster;
	vector<double> VRGain;
	double lambda = speedOfLightVac/285000000.0;
	
	for(uint i = 0; i < vr.size(); i++){
		vec VRPos;
		VRPos.ins(0,vr.at(i).x);
		VRPos.ins(1,vr.at(i).y);
		VRPos.ins(2,0);
		
		// TODO: get variables from ini files
		if(dist(VRPos,ms_pos) < 32.8 && (std::find(vr.at(i).bsIds.begin(), vr.at(i).bsIds.end(), bsId) != vr.at(i).bsIds.end())){
			double vrgain = calc_VRGain(lambda,16.8,32.8,VRPos,ms_pos);
			signalVRs.push_back(vr.at(i));
			VRGain.push_back(vrgain);
			std::cout << "VR number: " << i << std::endl;
			std::cout << "VR POS: " << VRPos << std::endl;
			std::cout << "VR Gain: " << vrgain << std::endl;
		}
	}
	
	for(uint i = 0; i < signalVRs.size(); i++){
		for(uint j = 0; j < cluster.size(); j++){
			if(signalVRs.at(i).clusterId == cluster.at(j).get_cluster_ID()){
				signalCluster.push_back(cluster.at(j));
			}
		}
	}
	
	//double pathloss = calc_path_loss_cost2100(dist(ms_pos,bsPosition));
	double pathloss = 3.746883049849024e-04;
	std::cout << "Pathloss (Cost 231 Walfish Ikegami Model): " << pathloss << std::endl;
	
	// Antennae spacing 1/2 wavelength
	mat antennae_pos = zeros(2,3);
	antennae_pos.set(1,2,lambda * 0.5);
	
	// Contains the double directional (Angles for BS and MS) impusle response (DDIR)
	mat channel_matrix;
	
	for(uint i = 0; i < signalCluster.size(); i++){
		mat channel = compute_channel_matrix(pathloss, signalCluster.at(i), num_mpc, num_dmc, bsPosition, ms_pos, 285000000.0, VRGain.at(i));
		//std::cout << "Current cluster: " << i << " " << channel << std::endl; 
		for(int j = 0; j < num_mpc + num_dmc; j++){
			channel_matrix.append_row(channel.get_row(j));
		}
		//std::cout << "cluster type: " << signalCluster.at(i).get_cluster_type() << std::endl;
		//std::cout << "Wavelength: " << lambda << std::endl;
		//std::cout << "Channel Matrix: " << std::endl << channel_matrix << std::endl << std::endl;
	}
	
	ofstream output;
	output.open("../Results_Validation/Results_Matlab_Channel_Info_CPP.txt");
	output.precision(15);
	for(int i = 0; i < channel_matrix.rows(); i++){
		for(int j = 0; j < channel_matrix.cols(); j++){
			if(channel_matrix(i,j) < 0){
				output << std::fixed << "  " << channel_matrix(i,j);
			} else{
				output << std::fixed << "   " << channel_matrix(i,j);
			}
		}
		output << endl << endl;
	}
	
	/*
	Local Clusters disabled for simplicity. Since they are computed exactly the same way as remote clusters, this isnt necessary.
	double vr_gain_local_ms = calc_VRGain(lambda,16.8,32.8,ms_pos,ms_pos);
	mat channel = compute_channel_matrix(pathloss, local_bs, num_mpc, num_dmc, bsPosition, ms_pos, 285000000.0, vr_gain_local_ms);
	for(int j = 0; j < num_mpc + num_dmc; j++){
		channel_matrix.append_row(channel.get_row(j));
	}
	
	double vr_gain_local_bs = calc_VRGain(lambda,16.8,32.8,bsPosition,ms_pos);
	channel = compute_channel_matrix(pathloss, local_ms, num_mpc, num_dmc, bsPosition, ms_pos, 285000000.0, vr_gain_local_bs);
	for(int j = 0; j < num_mpc + num_dmc; j++){
		channel_matrix.append_row(channel.get_row(j));
	}
	*/
	// Compute Line of Sight Component
	double power_los, power_factor_los;
	double d_ms_bs = dist(ms_pos,bsPosition);
	double LOS_cutoff_distance = 350.0;
	double VR_radius_los = 343.0;
	double TR_radius_los = 93.0;
	vec LOS_VR_pos;
	LOS_VR_pos.ins(0,LOS_VR.x);
	LOS_VR_pos.ins(1,LOS_VR.y);
	LOS_VR_pos.ins(2,0);
	
	double vr_los_gain = 0;
	// Power is Zero, if MS lies outside of the Line of Sight area
	if (d_ms_bs > LOS_cutoff_distance) {
		power_los = 0;
		power_factor_los = 0;
	} else {
		double d_ms_vr_los = dist(ms_pos, LOS_VR_pos);
		// Power is Zero, if MS lies outside of the Line of Sight VR
		if (d_ms_vr_los > VR_radius_los) {
			power_los = 0;
			power_factor_los = 0;
		} else {
			// The Line of Sight Power equals the sum off all cluster powers 
			// multiplied by a factor that is log normal distributed (Section 3.6.1.3)
			vr_los_gain = calc_VRGain(lambda, TR_radius_los, VR_radius_los, LOS_VR_pos, ms_pos); // Not used
			double power_other = 0;
			ifstream power ("power.txt");
			power >> power_factor_los;
			power_other = pow(sum(channel_matrix.get_col(5)), 2) + pow(sum(channel_matrix.get_col(6)), 2);
			power_los = abs(power_factor_los) * power_other;
		}
	}
	
	vec los_bs_sph = Cart_to_Sph(ms_pos - bsPosition);
	vec los_ms_sph = Cart_to_Sph(bsPosition - ms_pos);

	double d_bs_ms_xy = dist(bsPosition, ms_pos);		// distance BS to MS, in X-Y plane
	double tau_0_xy = d_bs_ms_xy / speedOfLightVac;							// delay of the LOS
	complex<double> phase_los(0,  -2 * M_PI * 285000000.0 * tau_0_xy);
	complex<double> channel_amp_los = sqrt(power_los) * exp(phase_los);
	double channel_amp_los_real = channel_amp_los.real();
	double channel_amp_los_imag = channel_amp_los.imag();
	vec channel_los = concat(los_bs_sph.get(0, 1), los_ms_sph.get(0, 1), to_vec(tau_0_xy), to_vec(channel_amp_los_real), to_vec(channel_amp_los_imag));
	channel_matrix.append_row(channel_los);
	
	// Generate own final Matrix for each ressourceblock frequency
	// In principle this Matrix is 7-dimensional: BS X MS X Cluster X MPC X RX_Antenna X TX_Antenna X Frequency
	// Since BS, MS, are fixed here and all MPC/Cluster are added up, we end up with a 3-dimensional Matrix
	// According to Formula 4.2
	
	double startfrequency = 275000000.0;
	double endfrequency = 295000000.0;
	double subcarrier_bandwidth = 78125.0;			
	int tx_num = 2;									// Number of transmitter antennae
	int rx_num = 2;									// Number of receiver antennae
	
	for(int f = startfrequency; f < endfrequency;f = f + subcarrier_bandwidth){
		cmat H_current(tx_num,rx_num);
		for(int i = 0; i < tx_num; i++){
			for(int j = 0; j < rx_num; j++){
				H_current.set(i,j,complex<double>(0,0));
			}
		}
		for (int m = 0; m < channel_matrix.rows(); m++ ) {
			complex<double> delay_phase(0, -2 * pi * f * channel_matrix.get(m, 4));
			complex<double> delay_response = exp(delay_phase);
			complex<double> channel_amplitude(channel_matrix.get(m, 5), channel_matrix(m, 6));
			
			cmat MIMO_mat = calc_MIMO_Channel_Matrix(2, 2, antennae_pos,antennae_pos, channel_matrix.get(m, 2), channel_matrix.get(m, 3), channel_matrix.get(m, 0), channel_matrix.get(m, 1), lambda);
			
			// Final formula 4.2
			H_current = H_current + (MIMO_mat * delay_response * channel_amplitude);
		}
		//std::cout << std::endl << "Composite Channel Matrix for frequency " << f << ": " << std::endl << H_current << std::endl;
		H.push_back(H_current);
	}

	std::cout << "[AngleBS, AngleMS, Delay, Amplitude_real, Amplitude_imag]" << std::endl << std::endl;
	//std::cout << channel_matrix << std::endl;
	std::cout << "Mobilestation Position: " << ms_pos << std::endl;
	std::cout << "Basestation Position: (0/0)" << std::endl;

	//! END Generate Signal Matrix
	return H;
}

/*
 * Computes a Matrix, that contains angles, delays and amplitude for all MPC's of a given cluster.
 * Formula 3.41 / 3.42 from the book.
 */
mat ChannelCalc::compute_channel_matrix(double pathloss, Cluster c, int num_mpc, int num_dmc, vec bs_pos, vec ms_pos, double center_freq, double VR_gain){
	mat mpc_bs_pos = c.get_cluster_MPC().get_MPC_pos(BS);
	mat mpc_ms_pos = c.get_cluster_MPC().get_MPC_pos(MS);
	mat channel_matrix;
	
	double cluster_delay = (dist(c.get_cluster_pos(BS), bs_pos) + dist(c.get_cluster_pos(MS), ms_pos)) / speedOfLightVac + c.get_cluster_link_delay();
	double LOS_delay = dist(ms_pos, bs_pos) / speedOfLightVac;
	
	double cluster_attenuation = attenuation(2.786127962522795, LOS_delay, 2.4e-6, cluster_delay);
	double cluster_shadow_fading = c.get_cluster_shadowing_fading();
	cout << "cluster_attenuation " << cluster_attenuation << endl;
	cout << "cluster_shadow_fading " << cluster_shadow_fading << endl;
	
	//cout << "Cluster Type: " << c.get_cluster_type() << endl;
	//cout << "C_LINK_DELAY: " << c.get_cluster_link_delay() << endl;
	
	// Specular MPC
	for(int i = 0; i < num_mpc; i++){
		double mpc_delay = (dist(mpc_bs_pos.get_row(i), bs_pos) + dist(mpc_ms_pos.get_row(i), ms_pos)) / speedOfLightVac + c.get_cluster_link_delay();
		complex<double> mpc_phase = exp(complex<double>(0,-2.0 * M_PI * mpc_delay * center_freq));
		
		// MPC complex amplitude according to formula 3.41
		cout.precision(16);
		cout << c.get_cluster_MPC().get_MPC_amplitude()(i) << endl;
		complex<double> mpc_amplitude = VR_gain * pathloss * c.get_cluster_MPC().get_MPC_amplitude()(i) * sqrt(cluster_attenuation*cluster_shadow_fading) * mpc_phase;
		//cout << endl << "Pathloss: " << pathloss << endl;
		//cout << "Rayleight Amp: " << c.get_cluster_MPC().get_MPC_amplitude()(i) << endl;
		//cout << "Attenuation * Fading: " << sqrt(cluster_attenuation*cluster_shadow_fading) << endl;
		//cout << "mpc_phase: " << mpc_phase << endl;
		//cout << "Result: " << mpc_amplitude << endl << endl;
		
		vec channel = concat(Cart_to_Sph(mpc_bs_pos.get_row(i) - bs_pos).get(0, 1), Cart_to_Sph(mpc_ms_pos.get_row(i) - ms_pos).get(0, 1), to_vec(mpc_delay), to_vec(mpc_amplitude.real()), to_vec(mpc_amplitude.imag()));
		
		// Contains MPC angle to BS, MPC angle to MS, MPC delay, MPC complex amplitude forall MPC of the given cluster.
		channel_matrix.append_row(channel);
	}
	
	//mpc_bs_pos = c.get_cluster_DMC().get_MPC_pos(BS);
	//mpc_ms_pos = c.get_cluster_DMC().get_MPC_pos(MS);
	// TODO: add delay
	/*
	// Dense/Diffuse Multipaths
	// Missing within the matlab implementation and therefore commented out for now.
	for(int i = 0; i < num_dmc; i++){
		// Exact Delay for this multipath component.
		//double dmc_delay = (distance(mpc_bs_pos.get_row(i), bs_pos) + distance(mpc_ms_pos.get_row(i), ms_pos)) / speedOfLightVac + c.get_cluster_link_delay() + c.get_cluster_DMC().get_DMC_delay_add()(i);
		double dmc_delay = (distance(mpc_bs_pos.get_row(i), bs_pos) + distance(mpc_ms_pos.get_row(i), ms_pos)) / speedOfLightVac + c.get_cluster_link_delay();
		
		complex<double> mpc_phase = exp(complex<double>(0,-2.0 * M_PI * dmc_delay * center_freq));
		
		// MPC complex amplitude according to formula 3.41
		complex<double> mpc_amplitude = pathloss * c.get_cluster_DMC().get_MPC_amplitude()(i) * sqrt(cluster_attenuation*cluster_shadow_fading) * mpc_phase;
		
		vec channel = concat(Cart_to_Sph(mpc_bs_pos.get_row(i) - bs_pos).get(0, 1), Cart_to_Sph(mpc_ms_pos.get_row(i) - ms_pos).get(0, 1), to_vec(dmc_delay), to_vec(mpc_amplitude.real()), to_vec(mpc_amplitude.imag()));
		
		// Contains MPC angle to BS, MPC angle to MS, MPC delay, MPC complex amplitude forall MPC of the given cluster.
		channel_matrix.append_row(channel);
	}
	*/
	return channel_matrix;
}

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

/*
 * Computes the VR Gain according to formula 3.35
 */
double ChannelCalc::calc_VRGain(double lambda, double TR_Radius, double VR_Radius, vec VR_Center, vec MS_Pos){
	double result = 0.5 - (1 / M_PI) * (atan((2.0 * sqrt(2.0) * (dist(VR_Center,MS_Pos) - VR_Radius + TR_Radius)) / sqrt(lambda * TR_Radius)));
	return result;
}

/*
 * Computes the general MIMO channel Matrix for arbitrary antennae
 * under the plane wave and balanced narrowband assumption excluding
 * the MPC amplitude according to formula 3.34.
 * Correctness verified for the following Scenario within matlab:
 * elevation angle always 90°, delta_y for antenna always 0. (2D Scenario)
 */
cmat ChannelCalc::calc_MIMO_Channel_Matrix(	int num_rx_antennae, int num_tx_antennae,
								mat rx_positions, mat tx_positions,
								double azimuthOfArrival, double elevationOfArrival,
								double azimuthOfDeparture, double elevationOfDeparture,	
								double lambda){
	
	//std::cout << std::endl << "Azimuth of Arrival: " << azimuthOfArrival << std::endl;
	//std::cout << "Elevation of Arrival: " << elevationOfArrival << std::endl;
	//std::cout << "Azimuth of Departure: " << azimuthOfDeparture << std::endl;
	//std::cout << "Elevation of Departure: " << elevationOfDeparture << std::endl << std::endl;
	
	// Correctness Check								
	assert( rx_positions.get_col(0).length() == num_rx_antennae );
	assert( tx_positions.get_col(0).length() == num_tx_antennae );
	assert( rx_positions.get_row(0).length() == 3 );
	assert( tx_positions.get_row(0).length() == 3 );
									
	// Compute unit direction of wave vector for arrival
	vec waveVectorDirection;
	waveVectorDirection.ins(0,sin(azimuthOfArrival) * cos(elevationOfArrival));
	waveVectorDirection.ins(1,sin(azimuthOfArrival) * sin(elevationOfArrival));
	waveVectorDirection.ins(2,cos(azimuthOfArrival));
	
	// Must be unit vector
	//std::cout << "RX: " << waveVectorDirection << std::endl;
	//assert( sqrt(pow(waveVectorDirection[0],2) + pow(waveVectorDirection[1],2) + pow(waveVectorDirection[2],2)) == 1 );
	
	// Compute the wave vector assuming a plane wave
	vec k_arrival =  ((2 * M_PI) / lambda) * waveVectorDirection;
	//std::cout << "k_arrival: " << k_arrival << std::endl;
	
	// Compute unit direction of wave vector for departure
	vec waveVectorDirection2;
	waveVectorDirection2.ins(0,sin(azimuthOfDeparture) * cos(elevationOfDeparture));
	waveVectorDirection2.ins(1,sin(azimuthOfDeparture) * sin(elevationOfDeparture));
	waveVectorDirection2.ins(2,cos(azimuthOfDeparture));
	
	// Must be unit vector
	//std::cout << "TX: " << waveVectorDirection2 << std::endl;
	//assert(sqrt(pow(waveVectorDirection2[0],2) + pow(waveVectorDirection2[1],2) + pow(waveVectorDirection2[2],2)) == 1);
	
	// Compute TX Steering vector
	vec k_departure =  ((2 * M_PI) / lambda) * waveVectorDirection2;
	//std::cout << "k_departure: " << k_departure << std::endl;
	
	// Compute Receiver (MS) Steering vector
	cvec rx_steering;
	for(int i = 0; i < rx_positions.get_col(0).length();i++){
		std::complex<double> entry(0,dot(k_arrival,rx_positions.get_row(i)));
		rx_steering.ins(i,exp(-1.0 * entry));
	}
	
	// Compute Transmitter (BS) Steering vector
	cvec tx_steering;
	for(int i = 0; i < tx_positions.get_col(0).length();i++){
		std::complex<double> entry(0,dot(k_departure,tx_positions.get_row(i)));
		tx_steering.ins(i,exp(-1.0 * entry));
	}
	//std::cout << "RX Steering: " << rx_steering << std::endl;
	//std::cout << "TX Steering: " << tx_steering << std::endl;
	cmat channel_matrix;
	//TODO remove double transpose (easy way to make the vector a matrix)
	channel_matrix = rx_steering.transpose().transpose() * tx_steering.transpose();
	return channel_matrix;
}

/* Calculates the distance between two points */
double ChannelCalc::distanceBetween(Position a, Position b)  {
    return sqrt(std::pow(b.x - a.x, 2) + std::pow(b.y - a.y, 2));
}

/*
 * We have a Macro Cell Scenario. According to page 129/131 in this case
 * the Cost231-Walfish-Ikegami Model is used. The formula used is defined 
 * in the book: "Digital Mobile Radio towards Future Generation Systems"
 * We assume the NLOS Case here. (Page 136ff from the book above.)
 * 
 * The following parameter restrictions exist:
 * Frequency: 800 ... 2000 MHz
 * Basestation height: 4 ... 50 m
 * Mobilestation height: 1 ... 3 m
 * Distance: 0.02 ... 5 km
 * 
 * Recommend default Values, if unknown:
 * Building Seperation: 20 ... 50 m
 * Road Width: 0.5 * Building Seperation
 * Roof height: 0 ... 3 m (0m if flat, 3m if pitched)
 * Rooftop Height: 3m * num_of_floors + roof height
 * 
 * This model is well suited if Basestation height >> rooftop height.
 * "Digital Mobile Radio towards Future Generation Systems" is avaible as free PDF!
 */
double ChannelCalc::calc_path_loss_cost2100(const double dist) {
	double pathloss;
	
	// For Testing the Cost 273 parameters defined on page 130 of the Cost2100 book are used.
	// The cell radius here is defined as 1000m which fits our scenario.
	// TODO: get from external
    double h_BS = 50;			// Basestation height
    double h_rooftop = 15;		// Rooftop height
    double freq_c = 2e9;		// center frequency
    double phi_road = 45;		// Road orientation
    double w_road = 25;			// Road width
    double w_street = 50;		// Building Seperation
    double h_MS = 1.5;			// Mobile station height
    
	double delta_base = h_BS - h_rooftop;
	
	// Free space Pathloss acording to formula 4.4.7 ("Digital Mobile Radio towards Future Generation Systems")
	double pl_0 = 32.45 + 20 * log10(dist / 1000) + 20 * log10(freq_c / 1e6); // Free space loss
	
	double pl_ori;
	
	// Road Orientation according to formula 4.4.9 ("Digital Mobile Radio towards Future Generation Systems")
	if (phi_road >= 0 && phi_road < 35) {
		pl_ori = -10+0.354*phi_road;
	} else if (phi_road >= 35 && phi_road < 55) {
		pl_ori = 2.5+0.075*(phi_road-35);
	} else if (phi_road >= 55 && phi_road < 90) {
		pl_ori = 4.0-0.114*(phi_road-55);
	}
	
	// "roof­top-­to-­street diffraction and scatter loss" according to formula 4.4.8 ("Digital Mobile Radio towards Future Generation Systems")
	double pl_lts = -16.9 - 10*log10(w_road) + 10 * log10(freq_c / 1e6) + 20 * log10(h_rooftop - h_MS) + pl_ori;
	
	// k_a = additional path loss, if BS attenae is belong rooftop of adjacent buildings
	// k_d = multi-screen diffraction loss versus distance
	// k_f = multi-screen diffraction loss versus frquency
	// Formulae 4.4.13, 4.4.14, 4.4.15 ("Digital Mobile Radio towards Future Generation Systems")
	double pl_bsh, k_a, k_d;
	if (h_BS > h_rooftop) {
		pl_bsh = -18 * log10(1 + delta_base);
		k_a = 54;
		k_d = 18;
	} else {
		pl_bsh = 0;
		k_d = 18 - 15 * delta_base / h_rooftop;
		if (dist >= 500) {
			k_a = 54 - 0.8 * delta_base;
		} else {
			k_a = 54 - 0.8 * delta_base * dist / 500;
		}
	}
	// Assume Metropolian city (for medium cities: k_f = -4 + 1.5 * (freq_c / 1e6 / 925-1) )
	double k_f = -4 + 0.7 * (freq_c / 1e6 / 925-1);
	
	// "multiple screen diffraction loss" according to formula 4.4.12 ("Digital Mobile Radio towards Future Generation Systems")
	double pl_msd = pl_bsh + k_a + k_d * log10(dist / 1000) + k_f * log10(freq_c / 1e6) - 9 * log10(w_street);
	
	// Add up different pieces of pathloss
	if (pl_msd + pl_lts >= 0 ) {
		pathloss = pl_0 + pl_msd + pl_lts;
	} else {
		pathloss = pl_0;
	}
	
	// We divide by 20 instead of ten, as this equals taking the squareroot afterwards.
	pathloss = pow(10.0, (-pathloss/20.0));

	// The squareroot of the amount of power remained.
	return pathloss;
}

/* Calculates the pathloss with the cost 231 model */
double ChannelCalc::calcPathLossCost231(Position a, Position b)  {
    double distance = distanceBetween(a, b);
    double log_f = std::log (carrierFreq / 1000000000) / 2.302;
    double C_H = 0.8 + ((1.11 * log_f) - 0.7) * msAntennaHeight - (1.56 * log_f);
    double log_BSH = std::log (bsAntennaHeight) / 2.303;

    // from the COST231 wiki entry
    // 2.303 is for the logarithm base change

    double loss_in_db = 46.3 + (33.9 * log_f) - (13.82 * log_BSH) - C_H +
                        ((44.9 - 6.55 * log_BSH) * std::log (distance) / 2.303);

    double loss = std::pow(10.0, loss_in_db / 10);

    //std::cout << "Distance: " << distance << "; PL: " << loss << std::endl;

    return loss;
}

/* Calculates the pathloss in free-space */
double ChannelCalc::calcPathLossFS(Position a, Position b)  {
    double distance = distanceBetween(a, b);
    //std::cout << "Distance: " << distance << std::endl;
    return std::pow((4 * M_PI * distance * carrierFreq) / speedOfLightVac, 2);
}

/* Uses the sender position and the interfernce positions */
double ChannelCalc::calcSINR()  {
    //the sending power of all devices is fix
    //calculate the pathloss between sender and target
    double senderInvPL = 1. / calcPathLossCost231(senderPosition, targetPosition);

    //calculate the sum of the inverse of all interferer pathloss at the target
    int numberOfInterferer = interfererPositions.size();
    double sumOfInvInterPL = 0;
    for(int i = 0; i < numberOfInterferer; i++)  {
        sumOfInvInterPL += (1. / calcPathLossCost231(targetPosition, interfererPositions.at(i)));
    }

    return senderInvPL / sumOfInvInterPL;
}

