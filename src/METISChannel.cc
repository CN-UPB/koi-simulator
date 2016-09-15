  /**
 * @file   METISChannel.cc
 * @Author Thomas Prinz (thomas.prinz@rwth-aachen.de)
 * @date   Januar, 2015
 * @brief  METIS Channel Subclass cc file
 *
 * Implementation of the METIS Channel subclass function.
 * Note on random numbers: This class is created once per BS.
 * Each BS is on a different node and per Definition by OMNeT that means
 * that it gets a different random seed. All MS of one BS are initialized
 * sequentially, which means their seed does not repeat as well.
 */

#include "METISChannel.h"
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <fstream>

using std::vector;
using std::tuple;

double getMSGain(double AoA, double EoA){ // Hertzian dipole; TODO: Change the return value to a vector for getting both theta and phi components
	//return (-1.0 * sin(EoA)); 
	return -1.0;				// For Z-oriented Hertzian dipole, F_phi is always zero and F_theta is -sin(theta), theta is in radians
}

double getBSGain(double AoD, double EoD){ // Hertzian dipole; TODO: Change the return value to a vector for getting both theta and phi components, EoD and AoD are in radians
	//return (-1.0 * sin(EoD));
	return -1.0;
}

/*
 * Cartesian: (x,y,z)
 * Spherical (Theta,Phi,r) [azimuth,elevation,r] (MatLab 'Convention')
 * The formula is identical to the Matlab intern one.
 */
inline vec Cart_to_Sph(vec const &input){
	vec output = zeros(3);
	output.set(0,atan((input(1) / input(0))));
	output.set(1,atan((input(2) / sqrt(pow(input(0),2) + pow(input(1),2)))));
	output.set(2,sqrt(pow(input(0),2) + pow(input(1),2) + pow(input(2),2)));
	return output;
}

/*
 * Cartesian: (x,y,z)
 * Spherical (Theta,Phi,r) [azimuth,elevation,r] (MatLab 'Convention')
 * The formula is identical to the Matlab intern one.
 */
inline vec Sph_to_Cart(vec const &input){
	vec output = zeros(3);
	output.set(0,input(2) * cos(input(1)) * cos(input(0)));
	output.set(1,input(2) * cos(input(1)) * sin(input(0)));
	output.set(2,input(2) * sin(input(1)));
	return output;
}

double METISChannel::ray_offset[20] = {	
		0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492, 
		0.3715, -0.3715, 0.5129, -0.5129, 0.6797, -0.6797,
		0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195,
		2.1551, -2.1551
	};

vector<vector<array<double,3>>> METISChannel::computeAntennaPos(
		const vector<Position>& transmitterPos,
		int numAntennas,
		double heightAntennas){
	vector<vector<array<double,3>>> antennaPos(transmitterPos.size(),
			vector<array<double,3>>(numAntennas));
	for(int i = 0; i < transmitterPos.size(); i++){
		for(int k = 0; k < (numAntennas/2); k++){
			antennaPos[i][k][0] = transmitterPos[i].x - (0.25 * wavelength * (numAntennas - 1 - (k*2.0)));
			antennaPos[i][k][1] = transmitterPos[i].y;
			antennaPos[i][k][2] = heightAntennas;
		}
		for(int k = (numAntennas/2); k < numAntennas ; k++){
			antennaPos[i][k][0] = antennaPos[i][k-1][0] + (0.5 * wavelength);
			antennaPos[i][k][1] = transmitterPos[i].y;
			antennaPos[i][k][2] = heightAntennas;
		}
	}
	return antennaPos;	
}

/**
* Function, that initializes all large scale and small scale parameters according to METIS specifications.
* @param module OMNeT++ module, which calls this function to allow later .ini access.
* @param msPositions Positions of the mobile stations.
* @param neighbourPositions Positions of the Neighbour BS.
* @return True if initialization was successful, false otherwise
*/
bool METISChannel::init(cSimpleModule* module, const vector<vector<Position>>& msPositions, std::map <int,Position> neighbourPositions){
	// Basic Initialization
	this->neighbourPositions = neighbourPositions;
	maxNumberOfNeighbours = module->par("maxNumberOfNeighbours");
	bsId = module->par("bsId");
	tti = module->par("tti");
    	numberOfMobileStations = module->par("numberOfMobileStations");
    	xPos = module->par("xPos");
	yPos = module->par("yPos");
	upRBs = module->par("upResourceBlocks");
	downRBs = module->par("downResourceBlocks");
	N_cluster_LOS = module->par("NumberOfClusters_LOS");
	N_cluster_NLOS = module->par("NumberOfClusters_NLOS");
	numOfRays_LOS = module->par("NumberOfRays_LOS");
	numOfRays_NLOS = module->par("NumberOfRays_NLOS");
   	initModule = module;
    	freq_c = module->par("CarrierFrequency");
    	SINRcounter = 0;
    	NumBsAntenna = module->par("NumBsAntenna");
    	NumMsAntenna = module->par("NumMsAntenna");
    	heightUE = module->par("OutdoorHeightUE");
   	heightBS = module->par("BsHeight");
	vel = module->par("Velocity");
	XPR_Mean_LOS = module->par("XPR_Mean_LOS");
	XPR_Std_LOS = module->par("XPR_Std_LOS");
	XPR_Mean_NLOS = module->par("XPR_Mean_NLOS");
	XPR_Std_NLOS = module->par("XPR_Std_NLOS");
    
    
    	// Find the neighbours and store the pair (bsId, position in data structures) in a map
    	cModule *cell = module->getParentModule()->getParentModule();
    	neighbourIdMatching = new NeighbourIdMatching(bsId, maxNumberOfNeighbours, cell);

	// Get Playground size from cell module:
	sizeX = cell->par("playgroundSizeX");
	sizeY = cell->par("playgroundSizeY");
    
    	// Actually, this counts the own BS as well, so substract 1 
    	numOfInterferers = neighbourIdMatching->numberOfNeighbours() - 1;
    
    
	// One Channel module per base station
    	// Half wavelength distance between antennas; give the position of Tx and Rx antennas in GCS
	// For even value of NumBsAntenna, the antenna elements will be equally spaced around the center of Tx
    	wavelength = speedOfLight / freq_c;
	vector<Position> tmpPos(neighbourPositions.size());
	for(size_t i = 0; i<neighbourPositions.size(); i++){
		tmpPos[i] = neighbourPositions[i];
	}
	bsAntennaPositions = computeAntennaPos(tmpPos,NumBsAntenna,
			heightBS);
		
	// Get position resend interval (Stationary MS assumed during this interval)
	timeSamples = module->par("positionResendInterval");
	// 4 Samples per TTI; for smooth Fourier transform
	timeSamples = timeSamples * 4; //scale according to positionResendInterval; timeSamples should be comparable to sim-time-limit
	timeVector = new double*[numberOfMobileStations];
	for(int i = 0; i < numberOfMobileStations; i++){
		timeVector[i] = new double[timeSamples];
		for(int t = 0; t < timeSamples; t++){
			timeVector[i][t] = 0.00025 * t;
		}
	}
	
	//compute initial SINR parameters
	recomputeMETISParams(msPositions);

	// There are a number of dynamically allocated member variables which 
	// are only allocated in this init method, which need only be freed 
	// in the destructor iff init has actually been called.
	initialized = true;
	
	return true;
}

void METISChannel::recomputeLargeScaleParameters(const vector<Position>& senders,
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
		vector<vector<double>>& sigma_sf_NLOS
		){
	double dist2D;
	std::cout << "Before Auto Correlation!" << std::endl;
	
	// Generate Autocorrelation
	vector<vector<vector<double>>> correlation;
	generateAutoCorrelation_LOS(senders,receivers,correlation);
	
	std::cout << "After Auto Correlation!" << std::endl;
       
	//TODO: Remove hardcoded, values taken after taking the square-root of the cross_matrix using sqrtm() in MATLAB
	double cross_matrix[7][7] = {
		{ 0.7249,   0.2510,   0.4678,  -0.0401,   0.1106,  -0.0903,  -0.4132},
    		{ 0.2510,   0.8598,   0.1428,   0.2893,   0.1575,  -0.2620,  -0.0143},
  		{ 0.4678,   0.1428,   0.8526,  -0.0093,  -0.0382,  -0.1763,  -0.0343},
  		{-0.0401,   0.2893,  -0.0093,   0.9552,  -0.0217,   0.0398,  -0.0122},
  		{ 0.1106,   0.1575,  -0.0382,  -0.0217,   0.9798,   0.0211,   0.0221},
 		{-0.0903,  -0.2620,  -0.1763,   0.0398,   0.0211,   0.9086,   0.2542},
 		{-0.4132,  -0.0143,  -0.0343,  -0.0122,   0.0221,   0.2542,   0.8733}};
		
    	// Transform Normal distributed random numbers to scenario specific distributions
    	double a = initModule->par("DS_mu_LOS"); 			// mean of delay spread
	double b = initModule->par("DS_eps_LOS"); 			// epsilon of delay spread
	double c = initModule->par("AoD_mu_LOS"); 			// mean of AoD
	double d = initModule->par("AoD_eps_LOS"); 		// epsilon of AoD
	double e = initModule->par("AoA_mu_LOS"); 			// mean of AoA
	double f = initModule->par("AoA_eps_LOS"); 		// epsilon of AoA
	double g = initModule->par("ZoA_mu_LOS"); 			// mean of ZoA
	double h = initModule->par("ZoA_eps_LOS"); 		// epsilon of ZoA
	double k = initModule->par("SF_sigma_LOS"); 		// sigma for shadow fading
	double l = initModule->par("K_mu"); 			// mean of Ricean K-factor
	double sigma_K = initModule->par("K_sigma");	// spread of K-factor
	
	// TODO: Take Square root
	// cross_matrix = sqrt(cross_matrix)
	for(size_t r = 0; r < receivers.size(); r++){
		for(size_t s = 0; s < senders.size(); s++){
			double ksi[7];
			for(int j = 0; j < 7;j++){
				ksi[j] = cross_matrix[j][0] * correlation[r][s][j] + cross_matrix[j][1] * correlation[r][s][j] + cross_matrix[j][2] * correlation[r][s][j] + cross_matrix[j][3] * correlation[r][s][j] + cross_matrix[j][4] * correlation[r][s][j] + cross_matrix[j][5] * correlation[r][s][j] + cross_matrix[j][6] * correlation[r][s][j];
			}
			sigma_ds_LOS[r][s]  = std::pow(10.0, (b*ksi[0] + a));      								// Log-Normal 
			sigma_asD_LOS[r][s] = std::min(104.0, std::pow(10.0, (d*ksi[1] + c)));      						// Log-Normal (maximum value should be 104 degrees) 
			sigma_asA_LOS[r][s] = std::min(104.0, std::pow(10.0, (f*ksi[2] + e)));      						// Log-Normal (maximum value should be 104 degrees) 
			sigma_zsA_LOS[r][s] = std::min(52.0, std::pow(10.0, (h*ksi[4] + g)));							// Log-Normal (maximum value should be 52 degrees) 
			dist2D = sqrt(pow((senders[s].x - receivers[r].x),2) + pow((senders[s].y - receivers[r].y),2));	
			sigma_zsD_LOS[r][s] = std::min(52.0, std::pow(10.0, (ksi[3] * sigma_ZSD(mean_ZSD(dist2D, heightUE, true), true))));	// Log-Normal (maximum value should be 52 degrees) 
			sigma_sf_LOS[r][s]  = std::pow(10.0, (0.1*k*ksi[5]));      								// Log-Normal dB
			sigma_kf_LOS[r][s]  = std::pow(10.0, (0.1*(sigma_K*ksi[6] + l)));	   						// Log-Normal dB
		}
	}

    	// Initialize large scale parameters (Non Line of Sight)
	std::cout << "Before Auto Correlation (NLOS)!" << std::endl;
	
	// Generate Autocorrelation
	correlation.clear();
	generateAutoCorrelation_NLOS(senders,receivers,correlation);
	
	std::cout << "After Auto Correlation (NLOS)!" << std::endl;
       
	//TODO: Remove hardcoded, values taken after taking the square-root of the cross_matrix using sqrtm() in MATLAB 
	double cross_matrix2[6][6] = {
		{ 0.8302,   0.0709,   0.1949,  -0.3252,  -0.0341,  -0.4011},
 		{ 0.0709,   0.9079,  -0.0277,   0.3008,   0.2808,   0.0259},
  		{ 0.1949,  -0.0277,   0.9580,   0.0352,   0.1131,  -0.1716},
  		{-0.3252,   0.3008,   0.0352,   0.8911,  -0.0542,  -0.0741},
  		{-0.0341,   0.2808,   0.1131,  -0.0542,   0.9509,  -0.0030},
  		{-0.4011,   0.0259,  -0.1716,  -0.0741,  -0.0030,   0.8964}};
		
    	// Transform Normal distributed random numbers to scenario specific distributions
    	a = initModule->par("DS_mu_NLOS"); 				// departure AS vs delay spread
	b = initModule->par("DS_eps_NLOS"); 			// arrival AS vs delay spread
	c = initModule->par("AoD_mu_NLOS"); 			// arrival AS vs shadowing std
	d = initModule->par("AoD_eps_NLOS"); 			// departure AS vs shadoving std
	e = initModule->par("AoA_mu_NLOS"); 			// delay spread vs shadoving std
	f = initModule->par("AoA_eps_NLOS"); 			// departure AS vs arrival AS
	g = initModule->par("SF_sigma_NLOS"); 			// departure AS vs k-factor
	h = initModule->par("ZoA_mu_NLOS"); 			// arrival AS vs k-factor
	k = initModule->par("ZoA_eps_NLOS"); 			// delay spread vs k-factor	
	
	// TODO: Take Square root
	// cross_matrix = sqrt(cross_matrix)
	for(size_t r = 0; r < receivers.size(); r++){
		for(size_t s = 0; s < senders.size(); s++){
				double ksi[6];
				for(int j = 0;j < 6;j++){
					ksi[j] = cross_matrix2[j][0] * correlation[r][s][j] + cross_matrix2[j][1] * correlation[r][s][j] + cross_matrix2[j][2] * correlation[r][s][j] + cross_matrix2[j][3] * correlation[r][s][j] + cross_matrix2[j][4] * correlation[r][s][j] + cross_matrix2[j][5] * correlation[r][s][j];
				}
				sigma_ds_NLOS[r][s]  = std::pow(10.0, (b*ksi[0] + a));      								// Log-Normal 
				sigma_asD_NLOS[r][s] = std::min(104.0, std::pow(10.0, (d*ksi[1] + c)));      						// Log-Normal (maximum value should be 104 degrees)  
				sigma_asA_NLOS[r][s] = std::min(104.0, std::pow(10.0, (f*ksi[2] + e)));      						// Log-Normal (maximum value should be 104 degrees) 
				sigma_zsA_NLOS[r][s] = std::min(52.0, std::pow(10.0, (k*ksi[4] + h)));							// Log-Normal (maximum value should be 52 degrees) 
				dist2D = sqrt(pow((senders[s].x - receivers[r].x),2) + pow((senders[s].y - receivers[r].y),2));
				sigma_zsD_NLOS[r][s] = std::min(52.0, std::pow(10.0, (ksi[3] * sigma_ZSD(mean_ZSD(dist2D, heightUE, false), false))));	// Log-Normal (maximum value should be 52 degrees) 
				sigma_sf_NLOS[r][s]  = std::pow(10.0, (0.1*g*ksi[5]));      								// Log-Normal dB
		}
	} 	

}

/**
 * Determines for each sender/reciever pair whether they have a line of sight
 */
vector<vector<bool>> METISChannel::genLosCond(const vector<Position>& sendPos,
		const vector<Position>& receivePos){
	vector<vector<bool>> losCond(receivePos.size(),
			vector<bool>(sendPos.size()));
	double dist;
	for(size_t i = 0; i < receivePos.size(); i++){
		for(size_t j=0; j<sendPos.size(); j++){
			dist = sqrt(pow(sendPos[j].x - receivePos[i].x,2) 
					+ pow(sendPos[j].y - receivePos[i].y,2));
			losCond[i][j] = LineOfSight(dist);
		}
	}
	return losCond;
}

tuple<vector<vector<vector<double>>>,vector<vector<vector<double>>>>
METISChannel::recomputeClusterDelays(const vector<vector<bool>>& LOSCondition,
		const vector<vector<double>>& sigmaDS_LOS,
		const vector<vector<double>>& sigmaDS_NLOS,
		const vector<vector<double>>& sigmaKF_LOS
		){
	vector<vector<vector<double>>> clusterDelays_NLOS(LOSCondition.size(),
			vector<vector<double>>(LOSCondition[0].size(),
				vector<double>())); 
	vector<vector<vector<double>>> clusterDelays_LOS(LOSCondition.size(),
			vector<vector<double>>(LOSCondition[0].size(),
				vector<double>())); 
    
	int n_clusters;
	double delayScaling;
	double min_delay;
	double sigma_ds;
    	for(size_t i = 0; i < LOSCondition.size(); i++){
		for(size_t j = 0; j<LOSCondition[i].size(); j++){
			if(LOSCondition[i][j]){
				// LOS Condition
				n_clusters = N_cluster_LOS;
				delayScaling = initModule->par("DelayScaling_LOS");
				sigma_ds = sigmaDS_LOS[i][j];
			}
			else{
				// NLOS Condition
				n_clusters = N_cluster_NLOS;
				delayScaling = initModule->par("DelayScaling_NLOS");
				sigma_ds = sigmaDS_NLOS[i][j];
			}
			clusterDelays_LOS[i][j].resize(n_clusters,0.0);
			clusterDelays_NLOS[i][j].resize(n_clusters,0.0);
			for(int k = 0; k < n_clusters; k++){
				clusterDelays_NLOS[i][j][k] = -1.0*delayScaling*sigma_ds*log(uniform(0,1));
			}
			min_delay = *std::min_element(clusterDelays_NLOS[i][j].cbegin(), clusterDelays_NLOS[i][j].cbegin() + n_clusters);
			// Normalize the delays (7.39)
			for(int k = 0; k < n_clusters; k++){
				clusterDelays_NLOS[i][j][k] = clusterDelays_NLOS[i][j][k] - min_delay;
			}
			// Sort the delays (7.39)
			std::sort(clusterDelays_NLOS[i][j].begin(), clusterDelays_NLOS[i][j].begin() + n_clusters, std::less<double>());
			if(LOSCondition[i][j]){
				// Apply LOS Peak compensation factor
				// Compute LOS Peak compensation factor (7.41)
				double K = 10.0 * log10(abs(sigmaKF_LOS[i][j]));
				double C_DS = 0.7705 - 0.0433 * K + 0.0002 * pow(K,2) + 0.000017 * pow(K,3);

				// Apply LOS compensation factor
				for(int k = 0; k < n_clusters; k++){
					clusterDelays_LOS[i][j][k] = clusterDelays_NLOS[i][j][k] / C_DS;
				}	
			}	
		}
	}
	return std::make_tuple(clusterDelays_LOS,clusterDelays_NLOS);
}

vector<vector<vector<double>>> METISChannel::genClusterPowers(const vector<vector<bool>>& LOSCondition,
		const vector<vector<vector<double>>>& clusterDelays,
		const vector<vector<double>>& sigmaDS_LOS,
		const vector<vector<double>>& sigmaDS_NLOS,
		const vector<vector<double>>& sigmaKF_LOS
		){
	vector<vector<vector<double>>> clusterPowers(clusterDelays.size(),
			vector<vector<double>>(clusterDelays[0].size(),
				vector<double>()));
	size_t numReceivers = clusterDelays.size();
	size_t numSenders = clusterDelays[0].size();
	double cluster_shadowing;
	double sum;
	double temp_CS;
	int n_clusters;
	double sigma_ds;
	double K;
	double P1_LOS;
	double delayScaling;
	
	// Generate cluster powers. 
	for(size_t i = 0; i < numReceivers; i++){
		for(size_t j = 0; j < numSenders; j++){
			if(LOSCondition[i][j]){
				temp_CS = initModule->par("PerClusterShadowing_LOS");
				delayScaling = initModule->par("DelayScaling_LOS");
				n_clusters = N_cluster_LOS;
				sigma_ds = sigmaDS_LOS[i][j];
			}
			else{
				temp_CS = initModule->par("PerClusterShadowing_NLOS");
				delayScaling = initModule->par("DelayScaling_NLOS");
				n_clusters = N_cluster_NLOS;
				sigma_ds = sigmaDS_NLOS[i][j];
			}
			clusterPowers[i][j].resize(n_clusters,0.0);
			sum = 0;
			cluster_shadowing = pow(temp_CS,2); //METIS
			for(int k = 0; k < n_clusters; k++){
				// Formula 7.42
				clusterPowers[i][j][k] = pow(10,-1.0*normal(0,cluster_shadowing)/10) * exp(-1.0 * clusterDelays[i][j][k] * ((delayScaling - 1) / (delayScaling * sigma_ds)) );
				sum += clusterPowers[i][j][k];
			}
			if(LOSCondition[i][j]){
				K = sigmaKF_LOS[i][j]; 			// K-factor in linear scale
				P1_LOS = K / (K + 1);
				for(int k = 0; k < n_clusters; k++){
					if (k==0){
						clusterPowers[i][j][k] = ((1/(K+1))*clusterPowers[i][j][k]/sum) + P1_LOS;
					}else{
						clusterPowers[i][j][k] = (1/(K+1))*clusterPowers[i][j][k]/sum;
					}
				}
			}else{
				for(int k = 0; k < n_clusters; k++){
					clusterPowers[i][j][k] /= sum;
				}
			}
		}
	}
	return clusterPowers;
}

vector<vector<vector<double>>> METISChannel::recomputeRayPowers(const vector<vector<bool>>& LOSCondition,
		vector<vector<vector<double>>>& clusterPowers
		){
	size_t numReceivers = LOSCondition.size();
	size_t numSenders = LOSCondition[0].size();
	vector<vector<vector<double>>> rayPowers(numReceivers,
			vector<vector<double>>(numSenders,vector<double>()));
	int n_clusters;
	int n_rays;
	for(size_t i = 0; i < numReceivers; i++){
		for(size_t j = 0; j<numSenders; j++){
			if(LOSCondition[i][j]){
				n_clusters = N_cluster_LOS;
				n_rays = numOfRays_LOS;
			}
			else{
				n_clusters = N_cluster_NLOS;
				n_rays = numOfRays_NLOS;
			}
			rayPowers[i][j].resize(n_clusters,0.0);
			for(int k = 0; k < n_clusters; k++){
				rayPowers[i][j][k] = clusterPowers[i][j][k] / n_rays;
			}
		}
	}
	return rayPowers;	
}

tuple<vector<vector<double>>,vector<vector<double>>,vector<vector<double>>,vector<vector<double>>> 
METISChannel::recomputeAngleDirection(
		const vector<Position>& receivers,
		const vector<Position>& senders,
		double heightSenders,
		double heightReceivers
		){
	size_t numReceivers = receivers.size();
	size_t numSenders = senders.size();
	vector<vector<double>> AoADir(numReceivers,
			vector<double>(numSenders));
	vector<vector<double>> ZoADir(numReceivers,
			vector<double>(numSenders));
	vector<vector<double>> AoDDir(numReceivers,
			vector<double>(numSenders));
	vector<vector<double>> ZoDDir(numReceivers,
			vector<double>(numSenders));
	double x_dir;
	double y_dir;
	vec cartLOS_RecToSend_Angle = zeros(3);
	vec cartLOS_SendToRec_Angle = zeros(3);
	vec sphLOSAngle;
		
	// Cycle through all mobile stations
    	for(size_t i = 0; i < numReceivers; i++){
		for(size_t j = 0; j < numSenders; j++){
				// Angles of Arrival
				x_dir = receivers[i].x - senders[j].x;
				y_dir = receivers[i].y - senders[j].y;

				cartLOS_RecToSend_Angle.set(0,x_dir);
				cartLOS_RecToSend_Angle.set(1,y_dir);
				cartLOS_RecToSend_Angle.set(2,heightReceivers);
				
				sphLOSAngle = Cart_to_Sph(cartLOS_RecToSend_Angle);
				AoADir[i][j] = sphLOSAngle.get(0);
				ZoADir[i][j] = sphLOSAngle.get(1);

				// Angles of Departure
				x_dir = senders[j].x - receivers[i].x;
				y_dir = senders[j].y - receivers[i].y;

				cartLOS_SendToRec_Angle.set(0,x_dir);
				cartLOS_SendToRec_Angle.set(1,y_dir);
				cartLOS_SendToRec_Angle.set(2,heightSenders);
				
				sphLOSAngle = Cart_to_Sph(cartLOS_SendToRec_Angle);
				AoDDir[i][j] = sphLOSAngle.get(0);
				ZoDDir[i][j] = sphLOSAngle.get(1);
		}
	}
	return std::make_tuple(AoADir,ZoADir,AoDDir,ZoDDir);
}

vector<vector<vector<vector<double>>>> METISChannel::recomputeAzimuthAngles(
		const vector<vector<bool>>& LOSCondition,
		const vector<vector<double>>& sigma_as_LOS,
		const vector<vector<double>>& sigma_as_NLOS,
		const vector<vector<double>>& sigma_kf,
		const vector<vector<vector<double>>>& clusterPowers,
		const vector<vector<double>>& angleDir,
		const bool arrival
		){
	size_t numReceivers = LOSCondition.size();
	size_t numSenders = LOSCondition[0].size();
	vector<vector<double>> azimuthCluster(numReceivers,
			vector<double>());
	vector<vector<vector<vector<double>>>> azimuthRays(numReceivers,
			vector<vector<vector<double>>>(numSenders,
				vector<vector<double>>()));
	int n_clusters;
	int n_rays;
	double X_1;
	double Y_1;
	double X_N;
	double Y_N;
	double sigma_angleSpread;
	double cluster_AS;
	for(size_t i = 0; i < numReceivers; i++){
		for(size_t j = 0; j < numSenders; j++){
			if(LOSCondition[i][j]){
				n_clusters = N_cluster_LOS;
				n_rays = numOfRays_LOS;
				sigma_angleSpread = sigma_as_LOS[i][j];
				if(arrival){
					cluster_AS = initModule->par("Cluster_ASA_LOS");
				}
				else{
					cluster_AS = initModule->par("Cluster_ASD_LOS");
				}
			}
			else{
				n_clusters = N_cluster_NLOS;
				n_rays = numOfRays_NLOS;
				sigma_angleSpread = sigma_as_NLOS[i][j];
				if(arrival){
					cluster_AS = initModule->par("Cluster_ASA_NLOS");
				}
				else{
					cluster_AS = initModule->par("Cluster_ASD_NLOS");
				}
			}
			azimuthRays[i][j].resize(n_clusters,vector<double>(
						n_rays));
			azimuthCluster[i].resize(n_clusters);
			X_1 = uniform(-1.0,1.0);
			Y_1 = normal(0.0, (pow(sigma_angleSpread / 7.0,2)));
			for(int k = 0; k < n_clusters; k++){
		
				//Formula 7.47
				azimuthCluster[i][k] = (sigma_angleSpread / (0.7 * C_AS(n_clusters, LOSCondition[i][j], i, sigma_kf))) * sqrt( -1.0 * log(clusterPowers[i][j][k] / *(std::max_element(clusterPowers[i][j].cbegin(), clusterPowers[i][j].cbegin() + n_clusters))) );
				X_N = uniform(-1.0,1.0);
				Y_N = normal(0.0, (pow(sigma_angleSpread / 7.0,2)));
		
				// First Ray has geometric LOS direction?!
				if(LOSCondition[i][j]){
					azimuthCluster[i][k] = azimuthCluster[i][k] * (X_N - X_1) + (Y_N - Y_1) + (angleDir[i][j]*180.0/pi);   // since AOA_LOS_dir is in radians
				}
				else{
					azimuthCluster[i][k] = azimuthCluster[i][k] * X_N + Y_N + (angleDir[i][j]*180.0/pi);   // since AOA_LOS_dir is in radians
				}
		
				for(int r = 0; r < n_rays; r++){
					// Final angle computation per ray 
					azimuthRays[i][j][k][r] = azimuthCluster[i][k] + cluster_AS * ray_offset[r];
				}
			}
		}
	}
	return azimuthRays;
}

/**
 * ATTENTION! This implementation simply sets all angles to 90°!
 * The computations specified in section 7.3.14.3 of METIS D1.2 are not yet 
 * implemented!
 */
vector<vector<vector<vector<double>>>> METISChannel::recomputeZenithAngles(
		const vector<vector<bool>>& LOSCondition,
		const vector<vector<double>>& sigma_zs_LOS,
		const vector<vector<double>>& sigma_zs_NLOS,
		const vector<vector<double>>& sigma_kf,
		const vector<vector<vector<double>>>& clusterPowers,
		const vector<vector<double>>& angleDir,
		const bool arrival
		){
	size_t numReceivers = LOSCondition.size();
	size_t numSenders = LOSCondition[0].size();
	vector<vector<vector<vector<double>>>> zenithRays(numReceivers,
			vector<vector<vector<double>>>(numSenders,
				vector<vector<double>>()));

	int n_clusters;
	int n_rays;
	for(size_t i = 0; i < numReceivers; i++){
		for(size_t j = 0; j<numSenders; j++){
			if(LOSCondition[i][j]){
				n_clusters = N_cluster_LOS;
				n_rays = numOfRays_LOS;
			}
			else{
				n_clusters = N_cluster_NLOS;
				n_rays = numOfRays_NLOS;
			}
			zenithRays[i][j].resize(n_clusters,vector<double>(n_rays));
			for(int k = 0; k < n_clusters; k++){
				for(int r = 0; r < n_rays; r++){
					zenithRays[i][j][k][r] = 90;
				}
			}
		}
	} 
	return zenithRays;
}

tuple<vector<vector<vector<vector<vector<double>>>>>,vector<vector<double>>>
METISChannel::genRandomPhases(
		const vector<vector<bool>>& LOSCondition
		){
	size_t numReceivers = LOSCondition.size();
	size_t numSenders = LOSCondition[0].size();
	vector<vector<vector<vector<vector<double>>>>> phases(numReceivers,
			vector<vector<vector<vector<double>>>>(numSenders,
				vector<vector<vector<double>>>()));
	vector<vector<double>> phases_LOS(numReceivers,vector<double>(
				numSenders));

	int n_clusters;
	int n_rays;
	for(size_t i = 0; i < numReceivers; i++){
		for(size_t j = 0; j < numSenders; j++){
			if(LOSCondition[i][j]){
				n_clusters = N_cluster_LOS;
				n_rays = numOfRays_LOS;
			}
			else{
				n_clusters = N_cluster_NLOS;
				n_rays = numOfRays_NLOS;
			}
			phases[i][j].resize(n_clusters,
					vector<vector<double>>(n_rays,
						vector<double>(4)));
			phases_LOS[i][j] = uniform(-1.0*pi, pi);
			for(int k = 0; k < n_clusters; k++){
				for(int r = 0; r < n_rays; r++){
					for(int l = 0; l < 4; l++){
						phases[i][j][k][r][l] = uniform(-1.0*pi, pi);			// for the random phases of NLOS component in equation 7-61 of METIS 1.2
					}
				}
			}
		}
	}
	return std::make_tuple(phases,phases_LOS);
}

vector<vector<vector<vector<double>>>> METISChannel::genCrossPolarization(
		vector<vector<bool>>& LOSCondition) {
	size_t numReceivers = LOSCondition.size();
	size_t numSenders = LOSCondition[0].size();
	vector<vector<vector<vector<double>>>> xpv(numReceivers,
			vector<vector<vector<double>>>(numSenders,
				vector<vector<double>>()));

	int n_clusters;
	int n_rays;
	for(size_t i = 0; i < numReceivers; i++){
		for(size_t j = 0; j < numSenders; j++){
			if(LOSCondition[i][j]){
				n_clusters = N_cluster_LOS;
				n_rays = numOfRays_LOS;
			}
			else{
				n_clusters = N_cluster_NLOS;
				n_rays = numOfRays_NLOS;
			}
			xpv[i][j].resize(n_clusters,vector<double>(n_rays));
			for(int k = 0; k < n_clusters; k++){
				for(int r = 0; r < n_rays; r++){
					xpv[i][j][k][r] = normal(XPR_Mean_LOS, pow(XPR_Std_LOS,2) );
				}
			}
		}
	}
	return xpv;
}

void METISChannel::computeRaySumCluster(
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
		vector<vector<vector<complex<double>>>>& raySum
		){
	double AoA[3];
	double AoD[3];
	complex<double> exp_arrival;
	complex<double> exp_departure;
	double receiverGain;
	double senderGain;
	complex<double> pol;
	complex<double> doppler;
	size_t rayIdx;
	// Cycle through all NLOS Rays
	for(size_t m = 0; m < numRays; m++){
		rayIdx = (subcluster!=nullptr) ? (*subcluster)[m] : m;
		AoA[0] = sin(zenithASA[rayIdx]*pi/180) * cos(azimuthASA[rayIdx]*pi/180);
		AoA[1] = sin(zenithASA[rayIdx]*pi/180) * sin(azimuthASA[rayIdx]*pi/180);
		AoA[2] = cos(zenithASA[rayIdx]*pi/180);

		AoD[0] = sin(zenithASD[rayIdx]*pi/180) * cos(azimuthASD[rayIdx]*pi/180);
		AoD[1] = sin(zenithASD[rayIdx]*pi/180) * sin(azimuthASD[rayIdx]*pi/180);
		AoD[2] = cos(zenithASD[rayIdx]*pi/180);

		exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * receiverAntennaPos[0] + AoA[1] * receiverAntennaPos[1] + AoA[2] * receiverAntennaPos[2])) );
		exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * senderAntennaPos[0] + AoD[1] * senderAntennaPos[1] + AoD[2] * senderAntennaPos[2])) );

		// TODO: Include Polarization
		// Replacement for polarization. (Below 4.14)
		//pol = exp( complex<double>(0,uniform(0,2*pi)));
		receiverGain = getMSGain(azimuthASA[rayIdx]*pi/180, zenithASA[rayIdx]*pi/180);
		senderGain = getBSGain(azimuthASD[rayIdx]*pi/180, zenithASD[rayIdx]*pi/180);
		pol = receiverGain * senderGain * exp(complex<double>(0, randomPhase[rayIdx][0]));

		// Calculate Doppler component of final formula
		// Formula: 
		// doppler = exp( j * k_0 * |v| * cos(AoA - AoMD) * t)
		// Where: 
		// k_0 = wavenumber
		// |v| = magnitude of velocity
		// AoA = Azimuth Angle of Arrival
		// AoMD = Azimuth Angle of MS Movement Direction
		// t = time vector

		// Cycle through the time axis (Formula 4.15)
		for(int t = 0; t < timeSamples; t++){
			// We do not have MS movement for the moment, so no doppler effect either
			//doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(azimuth_ASA[i][0][n][m]*pi/180 - AoMD) * timeVector[i][t] ) );
			doppler = 1.0;
			raySum[t][receiverAntennaIndex][senderAntennaIndex] += pol * doppler * exp_arrival * exp_departure;
		} // End time axis							
	} // End cycle Rays
	for(int t = 0; t < timeSamples; t++){
		raySum[t][receiverAntennaIndex][senderAntennaIndex] = prefactor * raySum[t][receiverAntennaIndex][senderAntennaIndex];
	} 							
}

tuple<vector<vector<vector<vector<vector<vector<complex<double>>>>>>>,
	vector<vector<vector<vector<vector<vector<complex<double>>>>>>>>
METISChannel::computeRaySums(vector<vector<bool>>& LOSCondition,
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
		){
	size_t numReceivers = LOSCondition.size();
	size_t numSenders = LOSCondition[0].size();
	vector<vector<vector<vector<vector<vector<complex<double>>>>>>> raySum(
			numReceivers,vector<vector<vector<vector<vector<complex<double>>>>>>(
				numSenders,vector<vector<vector<vector<complex<double>>>>>(
					N_cluster_NLOS + 4,vector<vector<vector<complex<double>>>>(
						timeSamples,vector<vector<complex<double>>>(
							numReceiverAntenna,vector<complex<double>>(
								numSenderAntenna))))));
	vector<vector<vector<vector<vector<vector<complex<double>>>>>>> raySum_LOS(
			numReceivers,vector<vector<vector<vector<vector<complex<double>>>>>>(
				numSenders,vector<vector<vector<vector<complex<double>>>>>(
					N_cluster_LOS + 4,vector<vector<vector<complex<double>>>>(
						timeSamples,vector<vector<complex<double>>>(
							numReceiverAntenna,vector<complex<double>>(
								numSenderAntenna))))));
	// Wavenumber k_0
	double k_0 = (2 * pi) / (speedOfLight / freq_c);
	
	// Create all local variables:
	complex<double> doppler, pol;
	double sq_P_over_M;
	double AoA[3];
	double AoD[3];
	double MSgain, BSgain;
	vector<vector<int>> subclusters{{0,1,2,3,4,5,6,7,18,19},
		{8,9,10,11,16,17},{12,13,14,15}};
	
	int n_clusters;
	// Cycle through all MS 
	for(size_t i = 0; i < numReceivers; i++){
		// Cycle through all base stations
		int clusterIdx;
		int clusterIdx_LOS;
		
		for(size_t j=0; j<numSenders; j++){
			double K = sigma_kf[i][j]; 
			// Cycle through all Receiver antennas (MS)
			for(int u = 0; u < numReceiverAntenna; u++){
				// Cycle through all Transmitter antennas (BS)
				for(int s = 0; s < numSenderAntenna; s++){
					// Cycle through all Paths/Clusters
					clusterIdx = 0;
					clusterIdx_LOS = 0;
					if(!LOSCondition[i][j]){
						n_clusters = N_cluster_NLOS;
					}
					else{
						n_clusters = N_cluster_LOS;
					}
					double P_n[3]; // Oversized, but 2 doubles really doesnt matter
					for(int n = 0; n < 2; n++){
						P_n[0] = 10.0/20.0 * clusterPowers[i][j][n];
						P_n[1] =  6.0/20.0 * clusterPowers[i][j][n];
						P_n[2] =  4.0/20.0 * clusterPowers[i][j][n];
						int id_subclust = 0;
						for(vector<int> &clust:subclusters){
							if(!LOSCondition[i][j]){
								sq_P_over_M = sqrt(P_n[id_subclust] / clust.size());
								computeRaySumCluster(clust.size(),
										sq_P_over_M,
										k_0,
										elevation_ASA[i][j][n],
										elevation_ASD[i][j][n],
										azimuth_ASA[i][j][n],
										azimuth_ASD[i][j][n],
										senderAntennaPos[j][s],
										receiverAntennaPos[i][u],
										u,
										s,
										randomPhase[i][j][n],
										&clust,
										raySum[i][j][clusterIdx]
										);
								clusterIdx++;
							}
							else{
								sq_P_over_M = (sqrt(1/(K + 1))) * sqrt(P_n[id_subclust] / clust.size());
								computeRaySumCluster(clust.size(),
										sq_P_over_M,
										k_0,
										elevation_ASA[i][j][n],
										elevation_ASD[i][j][n],
										azimuth_ASA[i][j][n],
										azimuth_ASD[i][j][n],
										senderAntennaPos[j][s],
										receiverAntennaPos[i][u],
										u,
										s,
										randomPhase[i][j][n],
										&clust,
										raySum_LOS[i][j][clusterIdx_LOS]
										);
								clusterIdx_LOS++;
								if(n == 0){ // for adding the additional LOS component, according to formula 7-61 in METIS 1.2
									AoA[0] = sin(ZoA_LOS_dir[i][j]) * cos(AoA_LOS_dir[i][j]);
									AoA[1] = sin(ZoA_LOS_dir[i][j]) * sin(AoA_LOS_dir[i][j]);
									AoA[2] = cos(ZoA_LOS_dir[i][j]);

									AoD[0] = sin(ZoD_LOS_dir[i][j]) * cos(AoD_LOS_dir[i][j]);
									AoD[1] = sin(ZoD_LOS_dir[i][j]) * sin(AoD_LOS_dir[i][j]);
									AoD[2] = cos(ZoD_LOS_dir[i][j]);

									complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * receiverAntennaPos[i][u][0] + AoA[1] * receiverAntennaPos[i][u][1] + AoA[2] * receiverAntennaPos[i][u][2])) );
									complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * senderAntennaPos[j][s][0] + AoD[1] * senderAntennaPos[j][s][1] + AoD[2] * senderAntennaPos[j][s][2])) );

									MSgain = getMSGain(AoA_LOS_dir[i][j], ZoA_LOS_dir[i][j]);
									BSgain = getBSGain(AoD_LOS_dir[i][j], ZoD_LOS_dir[i][j]);
									pol = MSgain * BSgain * exp(complex<double>(0, randomPhase_LOS[i][j]));

									// Cycle through the time axis (Formula 4.15)
									for(int t = 0; t < timeSamples; t++){
										//doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(AoA_LOS_dir[i][idIdx] - AoMD) * timeVector[i][t] ) );
										doppler = 1.0;
										raySum_LOS[i][j][0][t][u][s] += (sqrt(K / (K + 1))) * pol * doppler * exp_arrival * exp_departure;
									} // End time axis
								}
							}
							id_subclust++;
						}
					}
					double P;
					for(int n = 2; n < n_clusters; n++){
						P = clusterPowers[i][j][n];
						if(!LOSCondition[i][j]){
							sq_P_over_M = sqrt(P / numOfRays_NLOS);
							// Cycle through all NLOS Ray
							computeRaySumCluster(numOfRays_NLOS,
									sq_P_over_M,
									k_0,
									elevation_ASA[i][j][n],
									elevation_ASD[i][j][n],
									azimuth_ASA[i][j][n],
									azimuth_ASD[i][j][n],
									senderAntennaPos[j][s],
									receiverAntennaPos[i][u],
									u,
									s,
									randomPhase[i][j][n],
									nullptr,
									raySum[i][j][clusterIdx]
									);
							clusterIdx++;
						}
						else{
							sq_P_over_M = (sqrt(1/(K + 1))) * sqrt(P / numOfRays_LOS);
							// Cycle through all LOS Rays
							computeRaySumCluster(numOfRays_LOS,
									sq_P_over_M,
									k_0,
									elevation_ASA[i][j][n],
									elevation_ASD[i][j][n],
									azimuth_ASA[i][j][n],
									azimuth_ASD[i][j][n],
									senderAntennaPos[j][s],
									receiverAntennaPos[i][u],
									u,
									s,
									randomPhase[i][j][n],
									nullptr,
									raySum_LOS[i][j][clusterIdx_LOS]
									);
							clusterIdx_LOS++;
						}
					} // End cycle Clusters
//					if(i==2){
//						for(size_t clust=0; clust<raySum[i][j].size(); clust++){
//							for(size_t t=0; t<raySum[i][j][clust].size(); t++){
//								std::cout << "Ray Sum " << i << "," << j << "," << clust << "," << t << "," << u << "," << s << ": " << raySum[i][j][clust][t][u][s] << std::endl;
//							}
//							
//						}
//					}
				} // End Sender antenna
			} // End Receiver antenna
		} //End loop for all Senders
	} // End Links/Receivers

	return std::make_tuple(std::move(raySum),std::move(raySum_LOS));
}

vector<vector<vector<vector<double>>>> METISChannel::computeCoeffs(
		const vector<vector<bool>>& LOSCondition,
		const vector<Position>& receiverPos,
		const vector<Position>& senderPos,
		double heightReceivers,
		double heightSenders,
		bool up,
		int numRBs,
		int numReceiverAntenna,
		int numSenderAntenna,
		const vector<vector<vector<vector<vector<vector<complex<double>>>>>>>& raySum,
		const vector<vector<vector<vector<vector<vector<complex<double>>>>>>>& raySum_LOS,
		const vector<vector<vector<double>>>& clusterDelays,
		const vector<vector<vector<double>>>& clusterDelays_LOS
		){
	size_t numReceivers = receiverPos.size();
	size_t numSenders = senderPos.size();
	vector<vector<vector<vector<double>>>> coeffs(numReceivers,
			vector<vector<vector<double>>>(
				numSenders,vector<vector<double>>(
					timeSamples,vector<double>(numRBs))));
	double pathloss, dist3D, dist2D;
	int n_clusters;
	
	double delay_SC_1 = 5 * pow(10,-9); // delay for sub-cluster 1 (7-60)
	double delay_SC_2 = 10 * pow(10,-9); // delay for sub-cluster 2 (7-60)
	const vector<double> *delays;
	const vector<vector<vector<vector<complex<double>>>>> *sum;
	
	for(size_t i = 0; i < numReceivers; i++){ //TODO: Correct the implementation for each Tx-Rx antenna
		//Fourier Transform for interferers:
		for(size_t idIdx = 0; idIdx<numSenders; idIdx++){
			dist2D = sqrt(pow((senderPos[idIdx].x - receiverPos[i].x),2) + pow((senderPos[idIdx].y - receiverPos[i].y),2));
			dist3D = sqrt(pow(dist2D,2) + pow((heightSenders - heightReceivers),2));
			complex<double> res = complex<double>(0.0,0.0);
			pathloss = CalcPathloss(dist2D, dist3D, LOSCondition[i][idIdx]);
			for(int t = 0; t < timeSamples; t++){
				for(int f = 0; f < numRBs; f++){
					double freq_;
					if(up){
						// Computing values for UP
						// resource blocks
						freq_ = freq_c + (f+1)*180000;
					}
					else{
						// Computing values for DOWN
						// resource blocks
						freq_ = freq_c - f*180000;
					}
					res = complex<double>(0.0,0.0);
					if(LOSCondition[i][idIdx]){
						n_clusters = N_cluster_LOS;
						delays = &clusterDelays_LOS[i][idIdx];
						sum = &raySum_LOS[i][idIdx];
					}
					else{
						n_clusters = N_cluster_NLOS;
						delays = &clusterDelays[i][idIdx];
						sum = &raySum[i][idIdx];
					}
					for(int u = 0; u < numReceiverAntenna; u++){
						for(int s = 0; s < numSenderAntenna; s++){
							for(int n = 0; n < n_clusters; n++){
								if(n < 2){ // add additional sub-cluster delays at this stage (7-60 in METIS 1.2)
									res = res + (*sum)[n][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (*delays)[n]) );
									res = res + (*sum)[n+1][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * ((*delays)[n] + delay_SC_1)) );
									res = res + (*sum)[n+2][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * ((*delays)[n] + delay_SC_2)) );
									res = res + (*sum)[n+3][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (*delays)[n+1]) );
									res = res + (*sum)[n+4][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * ((*delays)[n+1] + delay_SC_1)) );
									res = res + (*sum)[n+5][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * ((*delays)[n+1] + delay_SC_2)) );
									n++;
								}else{
									res = res + (*sum)[n + 4][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (*delays)[n]) );
								}
							}
						}
					}
					double tempRes = pow(res.real(),2) + pow(res.imag(),2);
					coeffs[i][idIdx][t][f] = pathloss * tempRes;
				}
			}

		}
	}
	return coeffs;
}

void METISChannel::recomputeDownCoefficients(const vector<Position>& msPositions,
		const vector<Position>& bsPositions){
	int numReceiverAntenna = NumMsAntenna;
	int numSenderAntenna = NumBsAntenna;
    
    	// Copy MS Positions
	vector<Position> receiverPos(msPositions);	/*!< Position of the MS */
    	// Copy BS Positions
	vector<Position> senderPos(bsPositions);

	vector<vector<array<double,3>>> receiverAntennaPos(computeAntennaPos(
				receiverPos,numReceiverAntenna,heightUE));
	vector<vector<array<double,3>>>& senderAntennaPos = bsAntennaPositions;

	vector<vector<double>> AoA_LOS_dir;
	vector<vector<double>> ZoA_LOS_dir;
	vector<vector<double>> AoD_LOS_dir;
	vector<vector<double>> ZoD_LOS_dir;
	std::tie(AoA_LOS_dir,ZoA_LOS_dir,AoD_LOS_dir,ZoD_LOS_dir) = recomputeAngleDirection(
			receiverPos,
			senderPos,
			heightBS,
			heightUE
			);
	
	//Assign LOS Conditions:
	vector<vector<bool>> LOSCondition(genLosCond(senderPos,receiverPos));
	
	vector<vector<double>> sigma_ds_LOS(receiverPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_asD_LOS(receiverPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_asA_LOS(receiverPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_zsD_LOS(receiverPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_zsA_LOS(receiverPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_sf_LOS(receiverPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_kf_LOS(receiverPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_ds_NLOS(receiverPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_asD_NLOS(receiverPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_asA_NLOS(receiverPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_zsD_NLOS(receiverPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_zsA_NLOS(receiverPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_sf_NLOS(receiverPos.size(),
			vector<double>(neighbourPositions.size()));
	recomputeLargeScaleParameters(senderPos,receiverPos,
			sigma_ds_LOS,
			sigma_asD_LOS,
			sigma_asA_LOS,
			sigma_zsD_LOS,
			sigma_zsA_LOS,
			sigma_sf_LOS,
			sigma_kf_LOS,
			sigma_ds_NLOS,
			sigma_asD_NLOS,
			sigma_asA_NLOS,
			sigma_zsD_NLOS,
			sigma_zsA_NLOS,
			sigma_sf_NLOS
			);
	std::cout << "Finished Large Scale parameter.." << std::endl;

    	// Begin small scale parameter generation.
    
    	// Generate delays for each cluster according to Formula: 7:38 (METIS Document)
	vector<vector<vector<double>>> clusterDelays;
	vector<vector<vector<double>>> clusterDelays_LOS;
	std::tie(clusterDelays_LOS,clusterDelays) = recomputeClusterDelays(
			LOSCondition,
			sigma_ds_LOS,
			sigma_ds_NLOS,
			sigma_kf_LOS
			);

	vector<vector<vector<double>>> clusterPowers(genClusterPowers(LOSCondition,
				clusterDelays,
				sigma_ds_LOS,
				sigma_ds_NLOS,
				sigma_kf_LOS
				));

	// Precompute powers per ray (7.46)
	// While METIS D1.2 clearly states how to compute individual ray 
	// powers in section 7.3.13, only cluster powers are used 
	// for further computations.
//	vector<vector<vector<double>>> rayPowers(recomputeRayPowers(
//				LOSCondition,
//				clusterPowers)); /*!< The ray power for each cluster */
	
	
	
	// Generate azimuth angles of arrival
	vector<vector<vector<vector<double>>>> azimuth_ASA(recomputeAzimuthAngles(
			LOSCondition,
			sigma_asA_LOS,
			sigma_asA_NLOS,
			sigma_kf_LOS,
			clusterPowers,
			AoA_LOS_dir,
			true));

	// Generate azimuth angles of departure in the same way 
	vector<vector<vector<vector<double>>>> azimuth_ASD(recomputeAzimuthAngles(
			LOSCondition,
			sigma_asD_LOS,
			sigma_asD_NLOS,
			sigma_kf_LOS,
			clusterPowers,
			AoD_LOS_dir,
			false));

	// Generate Zenith angles 
	vector<vector<vector<vector<double>>>> elevation_ASA(recomputeZenithAngles(
				LOSCondition,
				sigma_zsA_LOS,
				sigma_zsA_NLOS,
				sigma_kf_LOS,
				clusterPowers,
				ZoA_LOS_dir,
				true));

	vector<vector<vector<vector<double>>>> elevation_ASD(recomputeZenithAngles(
				LOSCondition,
				sigma_zsD_LOS,
				sigma_zsD_NLOS,
				sigma_kf_LOS,
				clusterPowers,
				ZoD_LOS_dir,
				false));

	// Generate random phases (7.3.17)
	vector<vector<vector<vector<vector<double>>>>> randomPhase;
	vector<vector<double>> randomPhase_LOS;
	std::tie(randomPhase,randomPhase_LOS) = genRandomPhases(LOSCondition);
	
	// Generate cross polarization values
	vector<vector<vector<vector<double>>>> Xn_m(genCrossPolarization(
				LOSCondition));
	

	// initialize arrays for interferer ray sums
	vector<vector<vector<vector<vector<vector<complex<double>>>>>>> raySum_LOS;
	vector<vector<vector<vector<vector<vector<complex<double>>>>>>> raySum;

	std::cout << "START MAIN LOOP for BS: " << bsId << std::endl;

	std::tie(raySum,raySum_LOS) = computeRaySums(LOSCondition,
			sigma_kf_LOS,
			numReceiverAntenna,
			numSenderAntenna,
			clusterPowers,
			azimuth_ASA,
			azimuth_ASD,
			elevation_ASA,
			elevation_ASD,
			receiverAntennaPos,
			senderAntennaPos,
			randomPhase,
			randomPhase_LOS,
			AoA_LOS_dir,
			ZoA_LOS_dir,
			AoD_LOS_dir,
			ZoD_LOS_dir
			);
	
	std::cout << "FINISHED MAIN LOOP for BS: " << bsId << std::endl;

	//output2 << "Init 3 METIS at BS " << bsId << " with rand: " << normal(0,1) << std::endl;
		
	//------------------------------------------------------------------
	// Apply Fourier transform, to get time/frequency domain from time/delay
	
	std::cout << "START FOURIER TRANSFORM for BS: " << bsId << std::endl;

	coeffDownTable = std::move(computeCoeffs(
				LOSCondition,
				receiverPos,
				senderPos,
				heightUE,
				heightBS,
				false,
				downRBs,
				numReceiverAntenna,
				numSenderAntenna,
				raySum,
				raySum_LOS,
				clusterDelays,
				clusterDelays_LOS
				));
	
	std::cout << "FINISHED FOURIER TRANSFORM for BS: " << bsId << std::endl;
	
}

void METISChannel::recomputeUpCoefficients(const vector<vector<Position>>& msPositions,
		const vector<Position>& bsPositions){
	int numReceiverAntenna = NumBsAntenna;
	int numSenderAntenna = NumMsAntenna;
	// Copy BS Positions
	// This is only the position of the local base station, because 
	// that base station is the only receiver for the UP direction in any 
	// given cell.
	vector<Position> receiverPos{bsPositions[bsId]};
	vector<vector<array<double,3>>> receiverAntennaPos{bsAntennaPositions[bsId]};
	coeffUpTable = vector<vector<vector<vector<vector<double>>>>>(bsPositions.size(),
			vector<vector<vector<vector<double>>>>(1));
	for(size_t j=0; j<msPositions.size(); j++){
		// Copy MS Positions
		vector<Position> senderPos(msPositions[j]);

		vector<vector<array<double,3>>> senderAntennaPos(computeAntennaPos(
					senderPos,numSenderAntenna,heightUE));

		vector<vector<double>> AoA_LOS_dir;
		vector<vector<double>> ZoA_LOS_dir;
		vector<vector<double>> AoD_LOS_dir;
		vector<vector<double>> ZoD_LOS_dir;
		std::tie(AoA_LOS_dir,ZoA_LOS_dir,AoD_LOS_dir,ZoD_LOS_dir) = recomputeAngleDirection(
				receiverPos,
				senderPos,
				heightUE,
				heightBS
				);

		//Assign LOS Conditions:
		vector<vector<bool>> LOSCondition(genLosCond(senderPos,receiverPos));

		vector<vector<double>> sigma_ds_LOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_asD_LOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_asA_LOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_zsD_LOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_zsA_LOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_sf_LOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_kf_LOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_ds_NLOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_asD_NLOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_asA_NLOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_zsD_NLOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_zsA_NLOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_sf_NLOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		recomputeLargeScaleParameters(senderPos,receiverPos,
				sigma_ds_LOS,
				sigma_asD_LOS,
				sigma_asA_LOS,
				sigma_zsD_LOS,
				sigma_zsA_LOS,
				sigma_sf_LOS,
				sigma_kf_LOS,
				sigma_ds_NLOS,
				sigma_asD_NLOS,
				sigma_asA_NLOS,
				sigma_zsD_NLOS,
				sigma_zsA_NLOS,
				sigma_sf_NLOS
				);
		std::cout << "Finished Large Scale parameter.." << std::endl;

		// Begin small scale parameter generation.

		// Generate delays for each cluster according to Formula: 7:38 (METIS Document)
		vector<vector<vector<double>>> clusterDelays;
		vector<vector<vector<double>>> clusterDelays_LOS;
		std::tie(clusterDelays_LOS,clusterDelays) = recomputeClusterDelays(
				LOSCondition,
				sigma_ds_LOS,
				sigma_ds_NLOS,
				sigma_kf_LOS
				);

		vector<vector<vector<double>>> clusterPowers(genClusterPowers(LOSCondition,
					clusterDelays,
					sigma_ds_LOS,
					sigma_ds_NLOS,
					sigma_kf_LOS
					));

		// Precompute powers per ray (7.46)
		// While METIS D1.2 clearly states how to compute individual ray 
		// powers in section 7.3.13, only cluster powers are used 
		// for further computations.
		//	vector<vector<vector<double>>> rayPowers(recomputeRayPowers(
		//				LOSCondition,
		//				clusterPowers)); /*!< The ray power for each cluster */



		// Generate azimuth angles of arrival
		vector<vector<vector<vector<double>>>> azimuth_ASA(recomputeAzimuthAngles(
					LOSCondition,
					sigma_asA_LOS,
					sigma_asA_NLOS,
					sigma_kf_LOS,
					clusterPowers,
					AoA_LOS_dir,
					true));

		// Generate azimuth angles of departure in the same way 
		vector<vector<vector<vector<double>>>> azimuth_ASD(recomputeAzimuthAngles(
					LOSCondition,
					sigma_asD_LOS,
					sigma_asD_NLOS,
					sigma_kf_LOS,
					clusterPowers,
					AoD_LOS_dir,
					false));

		// Generate Zenith angles 
		vector<vector<vector<vector<double>>>> elevation_ASA(recomputeZenithAngles(
					LOSCondition,
					sigma_zsA_LOS,
					sigma_zsA_NLOS,
					sigma_kf_LOS,
					clusterPowers,
					ZoA_LOS_dir,
					true));

		vector<vector<vector<vector<double>>>> elevation_ASD(recomputeZenithAngles(
					LOSCondition,
					sigma_zsD_LOS,
					sigma_zsD_NLOS,
					sigma_kf_LOS,
					clusterPowers,
					ZoD_LOS_dir,
					false));

		// Generate random phases (7.3.17)
		vector<vector<vector<vector<vector<double>>>>> randomPhase;
		vector<vector<double>> randomPhase_LOS;
		std::tie(randomPhase,randomPhase_LOS) = genRandomPhases(LOSCondition);

		// Generate cross polarization values
		vector<vector<vector<vector<double>>>> Xn_m(genCrossPolarization(
					LOSCondition));


		// initialize arrays for interferer ray sums
		vector<vector<vector<vector<vector<vector<complex<double>>>>>>> raySum_LOS;
		vector<vector<vector<vector<vector<vector<complex<double>>>>>>> raySum;

		std::cout << "START MAIN LOOP for BS: " << bsId << std::endl;

		std::tie(raySum,raySum_LOS) = computeRaySums(LOSCondition,
				sigma_kf_LOS,
				numReceiverAntenna,
				numSenderAntenna,
				clusterPowers,
				azimuth_ASA,
				azimuth_ASD,
				elevation_ASA,
				elevation_ASD,
				receiverAntennaPos,
				senderAntennaPos,
				randomPhase,
				randomPhase_LOS,
				AoA_LOS_dir,
				ZoA_LOS_dir,
				AoD_LOS_dir,
				ZoD_LOS_dir
				);

		std::cout << "FINISHED MAIN LOOP for BS: " << bsId << std::endl;

		//output2 << "Init 3 METIS at BS " << bsId << " with rand: " << normal(0,1) << std::endl;

		//------------------------------------------------------------------
		// Apply Fourier transform, to get time/frequency domain from time/delay

		std::cout << "START FOURIER TRANSFORM for BS: " << bsId << std::endl;

		coeffUpTable[j] = std::move(computeCoeffs(
					LOSCondition,
					receiverPos,
					senderPos,
					heightBS,
					heightUE,
					true,
					upRBs,
					numReceiverAntenna,
					numSenderAntenna,
					raySum,
					raySum_LOS,
					clusterDelays,
					clusterDelays_LOS
					));

		std::cout << "FINISHED FOURIER TRANSFORM for BS: " << bsId << std::endl;
	
	}
    
	
}

void METISChannel::recomputeMETISParams(const vector<vector<Position>>& msPositions){
	int numMs;
	vector<vector<Position>> msPos(neighbourPositions.size(),
			vector<Position>());
	for(size_t j=0; j<neighbourPositions.size(); j++){
		numMs = neighbourIdMatching->getNumberOfMS(j);
		msPos[j].resize(numMs);
		for(size_t i=0; i<numMs; i++){
			msPos[j][i].x = msPositions[j][i].x;
			msPos[j][i].y = msPositions[j][i].y;
		}
	}
	vector<Position> bsPos(neighbourPositions.size());
	for(size_t j=0; j<neighbourPositions.size(); j++){
		bsPos[j] = neighbourPositions[j];
	}
	recomputeDownCoefficients(msPos[bsId],bsPos);
	recomputeUpCoefficients(msPos,bsPos);
	recomputeD2DCoefficients(msPos);
}

void METISChannel::recomputeD2DCoefficients(const vector<vector<Position>>& msPositions){
	int numReceiverAntenna = NumMsAntenna;
	int numSenderAntenna = NumMsAntenna;

	vector<Position> receiverPos{msPositions[bsId]};
	vector<vector<array<double,3>>> receiverAntennaPos(computeAntennaPos(
				receiverPos,numReceiverAntenna,heightUE));
	coeffUpD2DTable = vector<vector<vector<vector<vector<double>>>>>(msPositions.size(),
			vector<vector<vector<vector<double>>>>(numberOfMobileStations));
	coeffDownD2DTable = vector<vector<vector<vector<vector<double>>>>>(msPositions.size(),
			vector<vector<vector<vector<double>>>>(numberOfMobileStations));
	for(size_t j=0; j<msPositions.size(); j++){
		// Copy MS Positions
		vector<Position> senderPos(msPositions[j]);

		vector<vector<array<double,3>>> senderAntennaPos(computeAntennaPos(
					senderPos,numSenderAntenna,heightUE));

		vector<vector<double>> AoA_LOS_dir;
		vector<vector<double>> ZoA_LOS_dir;
		vector<vector<double>> AoD_LOS_dir;
		vector<vector<double>> ZoD_LOS_dir;
		std::tie(AoA_LOS_dir,ZoA_LOS_dir,AoD_LOS_dir,ZoD_LOS_dir) = recomputeAngleDirection(
				receiverPos,
				senderPos,
				heightUE,
				heightUE
				);

		//Assign LOS Conditions:
		vector<vector<bool>> LOSCondition(genLosCond(senderPos,receiverPos));

		vector<vector<double>> sigma_ds_LOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_asD_LOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_asA_LOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_zsD_LOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_zsA_LOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_sf_LOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_kf_LOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_ds_NLOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_asD_NLOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_asA_NLOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_zsD_NLOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_zsA_NLOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		vector<vector<double>> sigma_sf_NLOS(receiverPos.size(),
				vector<double>(senderPos.size()));
		recomputeLargeScaleParameters(senderPos,receiverPos,
				sigma_ds_LOS,
				sigma_asD_LOS,
				sigma_asA_LOS,
				sigma_zsD_LOS,
				sigma_zsA_LOS,
				sigma_sf_LOS,
				sigma_kf_LOS,
				sigma_ds_NLOS,
				sigma_asD_NLOS,
				sigma_asA_NLOS,
				sigma_zsD_NLOS,
				sigma_zsA_NLOS,
				sigma_sf_NLOS
				);
		std::cout << "Finished Large Scale parameter.." << std::endl;

		// Begin small scale parameter generation.

		// Generate delays for each cluster according to Formula: 7:38 (METIS Document)
		vector<vector<vector<double>>> clusterDelays;
		vector<vector<vector<double>>> clusterDelays_LOS;
		std::tie(clusterDelays_LOS,clusterDelays) = recomputeClusterDelays(
				LOSCondition,
				sigma_ds_LOS,
				sigma_ds_NLOS,
				sigma_kf_LOS
				);

		vector<vector<vector<double>>> clusterPowers(genClusterPowers(LOSCondition,
					clusterDelays,
					sigma_ds_LOS,
					sigma_ds_NLOS,
					sigma_kf_LOS
					));

		// Precompute powers per ray (7.46)
		// While METIS D1.2 clearly states how to compute individual ray 
		// powers in section 7.3.13, only cluster powers are used 
		// for further computations.
		//	vector<vector<vector<double>>> rayPowers(recomputeRayPowers(
		//				LOSCondition,
		//				clusterPowers)); /*!< The ray power for each cluster */



		// Generate azimuth angles of arrival
		vector<vector<vector<vector<double>>>> azimuth_ASA(recomputeAzimuthAngles(
					LOSCondition,
					sigma_asA_LOS,
					sigma_asA_NLOS,
					sigma_kf_LOS,
					clusterPowers,
					AoA_LOS_dir,
					true));

		// Generate azimuth angles of departure in the same way 
		vector<vector<vector<vector<double>>>> azimuth_ASD(recomputeAzimuthAngles(
					LOSCondition,
					sigma_asD_LOS,
					sigma_asD_NLOS,
					sigma_kf_LOS,
					clusterPowers,
					AoD_LOS_dir,
					false));

		// Generate Zenith angles 
		vector<vector<vector<vector<double>>>> elevation_ASA(recomputeZenithAngles(
					LOSCondition,
					sigma_zsA_LOS,
					sigma_zsA_NLOS,
					sigma_kf_LOS,
					clusterPowers,
					ZoA_LOS_dir,
					true));

		vector<vector<vector<vector<double>>>> elevation_ASD(recomputeZenithAngles(
					LOSCondition,
					sigma_zsD_LOS,
					sigma_zsD_NLOS,
					sigma_kf_LOS,
					clusterPowers,
					ZoD_LOS_dir,
					false));

		// Generate random phases (7.3.17)
		vector<vector<vector<vector<vector<double>>>>> randomPhase;
		vector<vector<double>> randomPhase_LOS;
		std::tie(randomPhase,randomPhase_LOS) = genRandomPhases(LOSCondition);

		// Generate cross polarization values
		vector<vector<vector<vector<double>>>> Xn_m(genCrossPolarization(
					LOSCondition));


		// initialize arrays for interferer ray sums
		vector<vector<vector<vector<vector<vector<complex<double>>>>>>> raySum_LOS;
		vector<vector<vector<vector<vector<vector<complex<double>>>>>>> raySum;

		std::cout << "START MAIN LOOP D2D for BS: " << bsId << std::endl;

		std::tie(raySum,raySum_LOS) = computeRaySums(LOSCondition,
				sigma_kf_LOS,
				numReceiverAntenna,
				numSenderAntenna,
				clusterPowers,
				azimuth_ASA,
				azimuth_ASD,
				elevation_ASA,
				elevation_ASD,
				receiverAntennaPos,
				senderAntennaPos,
				randomPhase,
				randomPhase_LOS,
				AoA_LOS_dir,
				ZoA_LOS_dir,
				AoD_LOS_dir,
				ZoD_LOS_dir
				);

		std::cout << "FINISHED MAIN LOOP D2D for BS: " << bsId << std::endl;

		//output2 << "Init 3 METIS at BS " << bsId << " with rand: " << normal(0,1) << std::endl;

		//------------------------------------------------------------------
		// Apply Fourier transform, to get time/frequency domain from time/delay

		std::cout << "START FOURIER TRANSFORM D2D for BS: " << bsId << std::endl;

		coeffUpD2DTable[j] = std::move(computeCoeffs(
					LOSCondition,
					receiverPos,
					senderPos,
					heightUE,
					heightUE,
					true,
					upRBs,
					numReceiverAntenna,
					numSenderAntenna,
					raySum,
					raySum_LOS,
					clusterDelays,
					clusterDelays_LOS
					));
		coeffDownD2DTable[j] = std::move(computeCoeffs(
					LOSCondition,
					receiverPos,
					senderPos,
					heightUE,
					heightUE,
					false,
					downRBs,
					numReceiverAntenna,
					numSenderAntenna,
					raySum,
					raySum_LOS,
					clusterDelays,
					clusterDelays_LOS
					));

		std::cout << "FINISHED FOURIER TRANSFORM D2D for BS: " << bsId << std::endl;
	
	}
    
	
}

/**
* Function that computes the scaling factor according to table 7.5 and formula 7.48
* @param numCluster The number of clusters per link
* @param LOS true iff LOS Condition, false otherwise
* @return Scaling factor
*/
double METISChannel::C_AS(int numCluster, bool LOS, int i,
		vector<vector<double>> sigma_kf_LOS){
	if(LOS){
		double K = 10.0 * log10(abs(sigma_kf_LOS[i][0]));
		double k_factor_LOS = 1.1035 - 0.028 * K - 0.002 * pow(K,2) + 0.0001 * pow(K,3);
		
		return C_AS(numCluster, false, i,sigma_kf_LOS) * k_factor_LOS;
	}else{
		switch(numCluster){
			case 4:
				return 0.779;
			case 5:
				return 0.860;
			case 8:
				return 1.018;
			case 10:
				return 1.090;
			case 11:
				return 1.123;
			case 12:
				return 1.146;
			case 14:
				return 1.190;
			case 15:
				return 1.211;
			case 16:
				return 1.226;
			case 19:
				return 1.273;
			case 20:
				return 1.289;
		};
	}
	
	// Wrong cluster number.
	return -1.0;
}

/**
* Function that computes the scaling factor according to table 7.7 and formula 7.52
* @param numCluster The number of clusters per link
* @param LOS true iff LOS Condition, false otherwise
* @return Scaling factor
*/
double METISChannel::C_ZS(int numCluster, bool LOS,
		vector<vector<double>> sigma_kf_LOS){
	int i = 0;
	if(LOS){
		double K = 10.0 * log10(abs(sigma_kf_LOS[i][0]));
		double k_factor_LOS = 1.3086 - 0.0339 * K - 0.0077 * pow(K,2) + 0.0002 * pow(K,3);
		
		return C_ZS(numCluster, false,sigma_kf_LOS) * k_factor_LOS;
	}else{
		switch(numCluster){
			case 12:
				return 1.104;
			case 19:
				return 1.184;
			case 20:
				return 1.178;
		};
	}
	
	// Wrong cluster number.
	return -1.0;
}

/**
* Function that computes the exponential auto correlation for LOS.
* @return True if succesful. False otherwise.
*/
void METISChannel::generateAutoCorrelation_LOS(const vector<Position>& senders,
		const vector<Position>& receivers,
		vector<vector<vector<double>>>& correlation){
	correlation.resize(receivers.size(),vector<vector<double>>(senders.size(),
				vector<double>(7)));
	
	// Read in LOS Decorrelation parameters.
	int decorr_LOS_DS = initModule->par("Decorr_LOS_DS");
	int decorr_LOS_ASD = initModule->par("Decorr_LOS_ASD");
	int decorr_LOS_ASA = initModule->par("Decorr_LOS_ASA");
	int decorr_LOS_ZSD = initModule->par("Decorr_LOS_ZSD");
	int decorr_LOS_ZSA = initModule->par("Decorr_LOS_ZSA");
	int decorr_LOS_SF = initModule->par("Decorr_LOS_SF");
	int decorr_LOS_K = initModule->par("Decorr_LOS_K");
	
	
	// Temporary variables
	double ***grid;
	double ***tmpX;
	double ***tmpY; 
	double filter[7][11]; 
	double sum[7][11]; 
	
	grid = new double**[7];
	tmpX = new double**[7];
	tmpY = new double**[7];
	
	for(int i = 0; i < 7; i++){
		for(int j = 0; j <= 10; j++){
			filter[i][j] = 0;
			sum[i][j] = 0;
		}
	}
	
	for(int i = 0; i < 7; i++){
		grid[i] = new double*[(int) sizeX+10];
		tmpX[i] = new double*[(int) sizeX+10];
		tmpY[i] = new double*[(int) sizeX+10];
		for(int j = 0; j < (int) sizeX+10; j++){
			grid[i][j] = new double[(int) sizeY+10];
			tmpX[i][j] = new double[(int) sizeY+10];
			tmpY[i][j] = new double[(int) sizeY+10];
		}
	}
		
	for(int i = 0; i <= 10; i++){
		filter[0][i] = exp((-1.0 * i) / decorr_LOS_DS);
		sum[0][0] += filter[0][i];
		filter[1][i] = exp((-1.0 * i) / decorr_LOS_ASD);
		sum[1][0] += filter[1][i];
		filter[2][i] = exp((-1.0 * i) / decorr_LOS_ASA);
		sum[2][0] += filter[2][i];
		filter[3][i] = exp((-1.0 * i) / decorr_LOS_ZSD);
		sum[3][0] += filter[3][i];
		filter[4][i] = exp((-1.0 * i) / decorr_LOS_ZSA);
		sum[4][0] += filter[4][i];
		filter[5][i] = exp((-1.0 * i) / decorr_LOS_SF);
		sum[5][0] += filter[5][i];
		filter[6][i] = exp((-1.0 * i) / decorr_LOS_K);
		sum[6][0] += filter[6][i];
	}
	
	for(int i = 0; i <= 10; i++){
		filter[0][i] = filter[0][i] / sum[0][0];
		filter[1][i] = filter[1][i] / sum[1][0];
		filter[2][i] = filter[2][i] / sum[2][0];
		filter[3][i] = filter[3][i] / sum[3][0];
		filter[4][i] = filter[4][i] / sum[4][0];
		filter[5][i] = filter[5][i] / sum[5][0];
		filter[6][i] = filter[6][i] / sum[6][0];
	}
		
	// Generate grid of (2*r+200)*(2*r+200) iid gaussian random numbers
	for(int i = 0; i < 7; i++){
		for(int j = 0; j < (int) sizeX+10; j++){
			for(int k = 0; k < (int) sizeY+10; k++){
				grid[i][j][k] = normal(0,1);
				tmpX[i][j][k] = 0;
				tmpY[i][j][k] = 0;
			}
		}
	}
		
	// Filter in X direction
	// For all Large scale parameters
	for(int i = 0; i < 7; i++){
		// Forall x points except offset
		for(int j = 10; j < (int) sizeX+10; j++){
			// Forall y points except offset
			for(int k = 10; k < (int) sizeY+10; k++){
				// Filter 100 points
				for(int l = 0; l < 10; l++){
					tmpX[i][j][k] += grid[i][j][k - l] * filter[i][l];
				}
			}
		}
	}
	
	// Filter in Y direction
	// For all Large scale parameters
	for(int i = 0; i < 7; i++){
		// For all x points except offset
		for(int j = 10; j < (int) sizeY+10; j++){
			// For all y points except offset
			for(int k = 10; k < (int) sizeX+10; k++){
				// Filter 100 points
				for(int l = 0; l < 10; l++){
					tmpY[i][k][j] += tmpX[i][k - l][j] * filter[i][l];
				}
			}
		}
	}
		
	// For each Link generate Large scale parameters
	for(size_t l = 0; l < receivers.size(); l++){
		for(size_t i=0; i<senders.size(); i++){
			for(int j=0; j<7; j++){
				correlation[l][i][j] = tmpY[j][(int) receivers[l].x ][(int) receivers[l].y];
			}
		}
	}
	// Clean up all local variables allocated on the heap
	// TODO Look into the possibility of actually putting those variables 
	// directly on the stack, instead of using entirely unnecessary 
	// dynamic memory allocation
	for(int i = 0; i < 7; i++){
		for(int j = 0; j < (int) sizeX+10; j++){
			delete[] grid[i][j];
			delete[] tmpX[i][j];
			delete[] tmpY[i][j];
		}
		delete[] grid[i];
		delete[] tmpX[i];
		delete[] tmpY[i];
	}
	delete[] grid;
	delete[] tmpX;
	delete[] tmpY;
}

/**
* Function that computes the exponential auto correlation for NLOS.
* @return True if succesful. False otherwise.
*/
void METISChannel::generateAutoCorrelation_NLOS(const vector<Position>& senders,
		const vector<Position>& receivers,
		vector<vector<vector<double>>>& correlation){
	correlation.resize(receivers.size(),vector<vector<double>>(senders.size(),
				vector<double>(6)));
	
	// Read in NLOS Decorrelation parameters.
	int decorr_NLOS_DS = initModule->par("Decorr_NLOS_DS");
	int decorr_NLOS_ASD = initModule->par("Decorr_NLOS_ASD");
	int decorr_NLOS_ASA = initModule->par("Decorr_NLOS_ASA");
	int decorr_NLOS_ZSD = initModule->par("Decorr_NLOS_ZSD");
	int decorr_NLOS_ZSA = initModule->par("Decorr_NLOS_ZSA");
	int decorr_NLOS_SF = initModule->par("Decorr_NLOS_SF");
	
	// Temporary variables
	double ***grid;
	double ***tmpX;
	double ***tmpY; 
	double **filter; 
	double **sum; 
	
	grid = new double**[6];
	tmpX = new double**[6];
	tmpY = new double**[6];
	filter = new double*[6];
	sum = new double*[6];
	
	for(int i = 0; i < 6; i++){
		filter[i] = new double[11];
		sum[i] = new double[11];
	}
	
	for(int i = 0; i < 6; i++){
		for(int j = 0; j <= 10; j++){
			filter[i][j] = 0;
			sum[i][j] = 0;
		}
	}
	
	for(int i = 0; i < 6; i++){
		grid[i] = new double*[(int) sizeX+10];
		tmpX[i] = new double*[(int) sizeX+10];
		tmpY[i] = new double*[(int) sizeX+10];
		for(int j = 0; j < (int) sizeX+10; j++){
			grid[i][j] = new double[(int) sizeY+10];
			tmpX[i][j] = new double[(int) sizeY+10];
			tmpY[i][j] = new double[(int) sizeY+10];
		}
	}
		
	for(int i = 0; i <= 10; i++){
		filter[0][i] = exp((-1.0 * i) / decorr_NLOS_DS);
		sum[0][0] += filter[0][i];
		filter[1][i] = exp((-1.0 * i) / decorr_NLOS_ASD);
		sum[1][0] += filter[1][i];
		filter[2][i] = exp((-1.0 * i) / decorr_NLOS_ASA);
		sum[2][0] += filter[2][i];
		filter[3][i] = exp((-1.0 * i) / decorr_NLOS_ZSD);
		sum[3][0] += filter[3][i];
		filter[4][i] = exp((-1.0 * i) / decorr_NLOS_ZSA);
		sum[4][0] += filter[4][i];
		filter[5][i] = exp((-1.0 * i) / decorr_NLOS_SF);
		sum[5][0] += filter[5][i];
	}
	
	for(int i = 0; i <= 10; i++){
		filter[0][i] = filter[0][i] / sum[0][0];
		filter[1][i] = filter[1][i] / sum[1][0];
		filter[2][i] = filter[2][i] / sum[2][0];
		filter[3][i] = filter[3][i] / sum[3][0];
		filter[4][i] = filter[4][i] / sum[4][0];
		filter[5][i] = filter[5][i] / sum[5][0];
	}
		
	// Generate grid of (2*r+200)*(2*r+200) iid gaussian random numbers
	for(int i = 0; i < 6; i++){
		for(int j = 0; j < (int) sizeX+10; j++){
			for(int k = 0; k < (int) sizeY+10; k++){
				grid[i][j][k] = normal(0,1);
				tmpX[i][j][k] = 0;
				tmpY[i][j][k] = 0;
			}
		}
	}
		
	// Filter in X direction
	// Forall Large scale parameters
	for(int i = 0; i < 6; i++){
		// Forall x points except offset
		for(int j = 10; j < (int) sizeX+10; j++){
			// Forall y points except offset
			for(int k = 10; k < (int) sizeY+10; k++){
				// Filter 100 points
				for(int l = 0; l < 10; l++){
					tmpX[i][j][k] += grid[i][j][k - l] * filter[i][l];
				}
			}
		}
	}
	
	// Filter in Y direction
	// Forall Large scale parameters
	for(int i = 0; i < 6; i++){
		// Forall x points except offset
		for(int j = 10; j < (int) sizeX+10; j++){
			// Forall y points except offset
			for(int k = 10; k < (int) sizeY+10; k++){
				// Filter 100 points
				for(int l = 0; l < 10; l++){
					tmpY[i][k][j] += tmpX[i][k - l][j] * filter[i][l];
				}
			}
		}
	}
		
	// For each Link generate Large scale parameters
	for(size_t l = 0; l < receivers.size(); l++){
		for(size_t i=0; i<senders.size(); i++){
			for(int j=0; j<6; j++){
				correlation[l][i][j] = tmpY[j][(int) receivers[l].x ][(int) receivers[l].y];
			}
		}
	}
	// Clean up all local variables allocated on the heap
	// TODO Look into the possibility of actually putting those variables 
	// directly on the stack, instead of using entirely unnecessary 
	// dynamic memory allocation
	for(int i = 0; i < 6; i++){
		for(int j = 0; j < (int) sizeX+10; j++){
			delete[] grid[i][j];
			delete[] tmpX[i][j];
			delete[] tmpY[i][j];
		}
		delete[] grid[i];
		delete[] tmpX[i];
		delete[] tmpY[i];
	}
	delete[] grid;
	delete[] tmpX;
	delete[] tmpY;
	for(int i=0; i<6; i++){
		delete[] filter[i];
		delete[] sum[i];
	}
	delete[] filter;
	delete[] sum;
}

/**
* Function that computes the pathloss for a given link.
* @param dist2D Distance in horizontal plane between BS and MS of the link.
* @param dist3D 3D distance between BS and MS of the link.
* @param LOS Boolean value which is true, iff the link is a LOS link.
* @return The pathloss for this link.
*/
double METISChannel::CalcPathloss(double dist2D, double dist3D, bool LOS){
	double pathloss;
	double distBP = 4 * (heightUE - 1.0) * (heightBS - 1.0) * (freq_c / speedOfLight); // Breakpoint Distance
	double pl_a;
	
	// TODO Reinsert minimum distance checks
	// The following METIS fomulas normally have a minimum distance 
	// requirement. To allow for easy random placement of mobile 
	// stations without having to account for that minimum distance
	// requirements, the checks were removed until another model 
	// can be implemented.
	if(LOS){
		if(heightUE >= 1.5 && ( dist2D < distBP )){
			// If the 2D distance is below the breakpoint distance (see Document 36873, page 23 for reference)
			pathloss = 22.0 * log10(dist3D) + 28.0 + 20.0 * log10(freq_c/1000000000);
		}else{
			pathloss = 22.0 * log10(dist3D) + 28.0 + 20.0 * log10(freq_c/1000000000) - 9.0 * log10( pow(distBP,2) + pow((heightBS - heightUE),2) );
		}
	}else{
		if(heightUE >= 1.5 && ( dist2D < distBP )){
			// If the 2D distance is below the breakpoint distance (see Document 36873, page 23 for reference)
			pathloss = 22.0 * log10(dist3D) + 28.0 + 20.0 * log10(freq_c/1000000000);
		}else{
			pathloss = 22.0 * log10(dist3D) + 28.0 + 20.0 * log10(freq_c/1000000000) - 9.0 * log10( pow(distBP,2) + pow((heightBS - heightUE),2) );
		}
		pathloss = std::max( pathloss, ( 36.7 * log10(dist3D) + 23.15 + 26.0 *log10(freq_c/1000000000) - 0.3 * (heightUE - 1.5) ) );
	}
	
	// Return path gain instead of loss
	//std::cout << "Distance: " << dist3D << std::endl;
	//std::cout << "Pathloss: " << 1/pl_a << std::endl;
	pl_a = pow(10,(pathloss/10));
	return 1/pl_a;
	//return 1/pow(10, (22*log10(dist2D) + 28 + 20*log10(freq_c / 1000000000)) /10);
}

/**
* Function that computes the line of sight probability.
* @param dist2D Distance in horizontal plane between BS and MS of the link.
* @return true iff LOS link, false otherwise
*/
bool METISChannel::LineOfSight(double dist2D){
	//double prob = std::min(18 / dist2D,1.0) * (1 - exp( -1.0*dist2D/36 )) + exp(-1.0 * dist2D/36);
	double prob;
	if (dist2D < 18.0){
		prob = 1.0;
	}else{
		prob = 	(18 / dist2D) + (1 - ( 18/dist2D )) * exp(-1.0 * dist2D/36);
	}
	double coin = uniform(0,1);
	return (coin <= prob); // check for the value of coin when running simulation
}

/**
* Function to compute the mean of ZSD
* @param dist2D Distance in horizontal plane between BS and MS of the link.
* @param heightUE Height of the UE in meters
* @param LOS Boolean value which is true, iff the link is a LOS link.
* @return Value of mean of ZSD
*/
double METISChannel::mean_ZSD(double dist2D, double heightUE, bool LOS){
	if (LOS){
		double mean = std::max(-1.0 * 0.5, (-1.0 * 2.1 * dist2D) - (0.01 * heightUE) + 0.765);
		return mean;
	}else{
		double mean = std::max(-1.0 * 0.5, (-1.0 * 2.1 * dist2D) - (0.01 * heightUE) + 0.915);
		return mean;
	}
}

/**
* Function to compute the sigma ZSD
* @param meanZSD The mean of ZSD
* @param LOS Boolean value which is true, iff the link is a LOS link.
* @return Value of sigma ZSD
*/
double METISChannel::sigma_ZSD(double meanZSD, bool LOS){
	if (LOS){
		double epsilon_ZSD = 0.4;
		double sigma;
		return sigma = normal(meanZSD, pow(epsilon_ZSD,2) );
	}else{
		double epsilon_ZSD = 0.49;
		double sigma;
		return sigma = normal(meanZSD, pow(epsilon_ZSD,2) );
	}
}

/**
* Function to allow the channel to receive message from other BS's
* @param msg Concrete message of potentially arbitrary subtype
*/
void METISChannel::handleMessage(cMessage* msg){
	if(msg->isName("COUNTER")){
		SINRcounter++;
		if(bsId == 3){
			ofstream myfile;
			myfile.open ("Counter.txt", std::ofstream::app);
			myfile << "SINR Counter METIS at BS " << bsId << " with rand: " << normal(0,1) << std::endl;
			//std::cout << "Counter: " << SINRcounter << std::endl;
		}
	}
	else if(msg->isName("DEBUG")){
		ofstream upTables;
		std::string fname("coeff_table_up_"+std::to_string(bsId)+".dat");
		upTables.open(fname,std::ofstream::trunc);
		printCoeffUpTables(upTables);
		upTables.close();
		ofstream downTables;
		fname = "coeff_table_down_"+std::to_string(bsId)+".dat";
		downTables.open(fname,std::ofstream::trunc);
		printCoeffDownTables(downTables);
		downTables.close();
	}
	delete msg;
}

double METISChannel::calcUpSINR(int RB, 
		std::forward_list<TransInfo*> &interferers,
		int msId,
		double transPower){
	int SINRCounter = 3; //originally set to std::round( simTime().dbl() * 1000.0* 4.0) - 1 
	double received = 0;
	double interference = 0;
	forward_list<TransInfo*>::iterator prev(interferers.before_begin());
	for(auto it = interferers.begin(); it!=interferers.end(); prev=it++){
		if((*it)->getCreationTime()>=simTime()-tti 
				&& ((*it)->getMessageDirection()==MessageDirection::up 
						|| (*it)->getMessageDirection()==MessageDirection::d2dUp)){
			// Sum interference arriving at the receiving base station
			interference += (*it)->getPower() * coeffUpTable[(*it)->getBsId()][0][(*it)->getMsId()][SINRCounter][RB];
		}
		else{
			delete *it;
			interferers.erase_after(prev);
			it=prev;
		}
	}
	received = transPower * coeffUpTable[bsId][0][msId][SINRcounter][RB];
	interference += getTermalNoise(300,180000);
	// Convert to db scale
	return 10 * log10( received / interference );
}

double METISChannel::calcDownSINR(int RB, 
		std::forward_list<TransInfo*> &interferers,
		int msId,
		double transPower){
	int SINRCounter = 3; //originally set to std::round( simTime().dbl() * 1000.0* 4.0) - 1 
	double received = 0;
	double interference = 0;
	forward_list<TransInfo*>::iterator prev(interferers.before_begin());
	for(auto it = interferers.begin(); it!=interferers.end(); prev=it++){
		if((*it)->getCreationTime()>=simTime()-tti){
			if((*it)->getMessageDirection()==MessageDirection::down){
				// Interference from a neighbouring BS
				interference += (*it)->getPower() * coeffDownTable[msId][(*it)->getBsId()][SINRCounter][RB];
			}
			else if((*it)->getMessageDirection()==MessageDirection::d2dDown){
				// Interfering from a MS transmitting D2D on the same down resource 
				// block
				interference += (*it)->getPower() * coeffDownD2DTable[(*it)->getBsId()][msId][(*it)->getMsId()][SINRCounter][RB];
			}
		}
		else{
			delete *it;
			interferers.erase_after(prev);
			it=prev;
		}
	}
	received = transPower * coeffDownTable[msId][bsId][SINRcounter][RB];
	interference += getTermalNoise(300,180000);
	// Convert to db scale
	return 10 * log10( received / interference );
}

double METISChannel::calcD2DSINR(int RB, 
		std::forward_list<TransInfo*> &interferers,
		int sendMsID,
		int receiveMsId,
		MessageDirection direction,
		double transPower){
	int SINRCounter = 3; //originally set to std::round( simTime().dbl() * 1000.0* 4.0) - 1 
	double received = 0;
	double interference = 0;
	forward_list<TransInfo*>::iterator prev(interferers.before_begin());
	for(auto it = interferers.begin(); it!=interferers.end(); prev=it++){
		if((*it)->getCreationTime()>=simTime()-tti){
			if(direction==MessageDirection::d2dDown){
				// Interferers when using the down frequency bands
				if((*it)->getMessageDirection()==MessageDirection::down){
					// The interferer is a neighbouring BS, link coefficients for those to all 
					// local MS are stored in the normal downlink table.
					interference += (*it)->getPower() * coeffDownTable[receiveMsId][(*it)->getBsId()][SINRCounter][RB];
				}
				else if((*it)->getMessageDirection()==MessageDirection::d2dDown){
					// Interference from a mobile station transmitting D2D on the same down 
					// resource block.
					interference += (*it)->getPower() * coeffDownD2DTable[(*it)->getBsId()][receiveMsId][(*it)->getMsId()][SINRCounter][RB];
				}
			}
			else if(direction==MessageDirection::d2dUp 
					&& ((*it)->getMessageDirection()==MessageDirection::up 
						|| (*it)->getMessageDirection()==MessageDirection::d2dUp)){
				interference += (*it)->getPower() * coeffUpD2DTable[(*it)->getBsId()][receiveMsId][(*it)->getMsId()][SINRCounter][RB];
			}
		}
		else{
			delete *it;
			interferers.erase_after(prev);
			it=prev;
		}
	}
	if(direction==MessageDirection::d2dDown){
		received = transPower * coeffDownD2DTable[bsId][receiveMsId][sendMsID][SINRCounter][RB];
	}
	else{
		received = transPower * coeffUpD2DTable[bsId][receiveMsId][sendMsID][SINRCounter][RB];
	}
	interference += getTermalNoise(300,180000);
	// Convert to db scale
	return 10 * log10( received / interference );
}

double METISChannel::calcAvgUpSINR(int RB, 
            std::forward_list<TransInfo*> &interferers,
            int msId,
            double transPower){
  // Compute SINR for MS->BS communication
  double res = calcUpSINR(RB,interferers,msId,transPower);
  // Add the average over all possible D2D connections to all local MS
  for(int i=0; i<numberOfMobileStations;++i){
    if(i!=msId){
      res += calcD2DSINR(RB,interferers,msId,i,MessageDirection::d2dUp,
          transPower);
    }
  }
  // Return the average over all SINR values
  // numberOfMobileStations is the number of values from the D2D computations,
  // +1 from the MS->BS computation.
  return res/(numberOfMobileStations+1);
}

double METISChannel::calcAvgDownSINR(int RB, 
            std::forward_list<TransInfo*> &interferers,
            double transPower){
  // Compute SINR for DOWN resource blocks for base stations.
  double res = 0.0;
  // Average over SINR values for all local mobile stations
  for(int i=0; i<numberOfMobileStations;++i){
    res += calcDownSINR(RB,interferers,i,transPower);
  }
  return res/numberOfMobileStations;
}

double METISChannel::calcAvgD2DDownSINR(int RB, 
            std::forward_list<TransInfo*> &interferers,
            int msId,
            double transPower){
  // Calculate average SINR for DOWN resource blocks when used for D2D by MS
  double res = 0.0;
  // Average over SINR values for all local mobile stations
  for(int i=0; i<numberOfMobileStations;++i){
    if(msId!=i){
      res += calcD2DSINR(RB,interferers,msId,i,MessageDirection::d2dDown,transPower);
    }
  }
  return res/numberOfMobileStations;
}

void METISChannel::updateChannel(const vector<vector<Position>>& msPos){
	recomputeMETISParams(msPos);
}

// Johnson Nyquist Noise
double METISChannel::getTermalNoise(double temp, double bandwidth){
	return (temp * bandwidth * 1.3806488e-23);
}

double METISChannel::calcPathloss(double dist){
	//Placeholder
	return 1.0;
}

std::ostream& METISChannel::printCoeffUpTables(std::ostream& out){
	out << "BS" << "\t" << "MS" << "\t" << "Time" << "\t" << "RB" << "\t" << "Coeff" << std::endl;
	for(size_t i=0; i<coeffUpTable.size(); i++){
		for(size_t j=0; j<coeffUpTable[i][0].size(); j++){
			for(size_t t=0; t<timeSamples; t++){
				for(size_t r=0; r<downRBs; r++){
					out << i << "\t" << j << "\t" << t << "\t" << r << "\t" << coeffUpTable[i][0][j][t][r] << std::endl;
				}

			}
		}
	}
	return out;
}

std::ostream& METISChannel::printCoeffDownTables(std::ostream& out){
	out << "MS" << "\t" << "BS" << "\t" << "Time" << "\t" << "RB" << "\t" << "Coeff" << std::endl;
	for(size_t i=0; i<coeffDownTable.size(); i++){
		for(size_t j=0; j<coeffDownTable[i].size(); j++){
			for(size_t t=0; t<timeSamples; t++){
				for(size_t r=0; r<downRBs; r++){
					out << i << "\t" << j << "\t" << t << "\t" << r << "\t" << coeffDownTable[i][j][t][r] << std::endl;
				}

			}
		}
	}
	return out;
}

METISChannel::~METISChannel(){
	//TODO: Delete all dynamic memory
	if(initialized){
		// All of the following member variables are only allocated 
		// in the init method, not the constructor. As a result,
		// the destructor might be called with all of the following 
		// variables uninitialized, leading to errors. Thus, they 
		// are only freed when init has been called, as indicated by 
		// the initialized variable.
		for(int i=0; i<numberOfMobileStations; i++){
			delete[] timeVector[i];
		}
		delete neighbourIdMatching; 
		delete[] timeVector;
	
	}

}
