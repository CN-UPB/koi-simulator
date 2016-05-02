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

/**
* Function, that initializes all large scale and small scale parameters according to METIS specifications.
* @param module OMNeT++ module, which calls this function to allow later .ini access.
* @param msPositions Positions of the mobile stations.
* @param neighbourPositions Positions of the Neighbour BS.
* @return True if initialization was successful, false otherwise
*/
bool METISChannel::init(cSimpleModule* module, Position** msPositions, std::map <int,Position> neighbourPositions){
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
    
    	OwnBsAntennaPosition = new double**[1];
    
    	// Find the neighbours and store the pair (bsId, position in data structures) in a map
    	cModule *cell = module->getParentModule()->getParentModule();
    	neighbourIdMatching = new NeighbourIdMatching(bsId, maxNumberOfNeighbours, cell);
    
    	// Actually, this counts the own BS as well, so substract 1 
    	numOfInterferers = neighbourIdMatching->numberOfNeighbours() - 1;
    
    
	// One Channel module per base station
    	// Half wavelength distance between antennas; give the position of Tx and Rx antennas in GCS
	// For even value of NumBsAntenna, the antenna elements will be equally spaced around the center of Tx
    	double wavelength = speedOfLight / freq_c;
    	for(int i = 0; i < 1; i++){
		OwnBsAntennaPosition[i] = new double*[NumBsAntenna];
		for(int j = 0; j < (NumBsAntenna/2); j++){
			OwnBsAntennaPosition[i][j] = new double[3];
			OwnBsAntennaPosition[i][j][0] = xPos - (0.25 * wavelength * (NumBsAntenna - 1 - (j*2.0)));
			OwnBsAntennaPosition[i][j][1] = yPos;
			OwnBsAntennaPosition[i][j][2] = heightBS;
		}

		for(int j = (NumBsAntenna/2); j < NumBsAntenna; j++){
			OwnBsAntennaPosition[i][j] = new double[3];
			OwnBsAntennaPosition[i][j][0] = OwnBsAntennaPosition[i][j-1][0] + (0.5 * wavelength);
			OwnBsAntennaPosition[i][j][1] = yPos;
			OwnBsAntennaPosition[i][j][2] = heightBS;
		}

	}
	
		
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
	
	// properly initialize the SINRtable with all zeros
	SINRtable = new double**[numberOfMobileStations];
	for(int i = 0; i < numberOfMobileStations; i++){
		SINRtable[i] = new double*[timeSamples];
		for(int j = 0; j < timeSamples; j++){
			SINRtable[i][j] = new double[downRBs];
			for(int k=0; k<downRBs; k++){
				SINRtable[i][j][k] = 0;
			}
		}
	}

	// properly initialize the SINRneighbour with all zeros
	SINRneighbour = new double***[numberOfMobileStations];
	for(int i = 0; i < numberOfMobileStations; i++){
		SINRneighbour[i] = new double**[numOfInterferers];
		for(int j = 0; j < numOfInterferers; j++){
			SINRneighbour[i][j] = new double*[timeSamples];
			for(int k = 0; k < timeSamples; k++){
				SINRneighbour[i][j][k] = new double[downRBs];
				for(int s=0; s<downRBs; s++){
					SINRneighbour[i][j][k][s] = 0;
				}
			}
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
					azimuthCluster[i][j] = azimuthCluster[i][k] * X_N + Y_N + (angleDir[i][j]*180.0/pi);   // since AOA_LOS_dir is in radians
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

void METISChannel::recomputeMETISParams(Position** msPositions){
    	int fromBsId = neighbourIdMatching->getDataStrId(bsId);
    	double wavelength = speedOfLight / freq_c;
	double dist2D;
    
    	// randomly generate direction and magnitude of movement vector
	double *MSVelMag = new double[numberOfMobileStations]; /*!< Magnitude of MS Velocity */
	double **MSVelDir = new double*[numberOfMobileStations]; /*!< Direction of MS Velocity */
    	// TODO: sync with OMNeT++ part
    	for(int i = 0; i < numberOfMobileStations; i++){
		MSVelDir[i] = new double[3];
		MSVelMag[i] = (vel*1000)/3600; // converting velocity from km/h to m/s 
		MSVelDir[i][0] = 1.0; // +1 for moving towards the BS, -1 for moving away from the BS
		MSVelDir[i][1] = 0.0; // set to zero b/c MS is not moving in y or z-directions
		MSVelDir[i][2] = 0.0; // same as above; originally all are set to uniform(-1.0,1.0)
		double dirMag = sqrt(pow(MSVelDir[i][0],2) + pow(MSVelDir[i][1],2) + pow(MSVelDir[i][2],2));
		MSVelDir[i][0] = MSVelDir[i][0] / dirMag;
		MSVelDir[i][1] = MSVelDir[i][1] / dirMag;
		MSVelDir[i][2] = MSVelDir[i][2] / dirMag;
	}

    	// Copy MS Positions
	vector<Position> MSPos(numberOfMobileStations);	/*!< Position of the MS */
    	for(int i = 0; i < numberOfMobileStations; i++){
		MSPos[i].x = msPositions[fromBsId][i].x;
		MSPos[i].y = msPositions[fromBsId][i].y;
	}
    	// Copy BS Positions
	vector<Position> BSPos(neighbourPositions.size());
	for(auto it = neighbourPositions.begin(); it!=neighbourPositions.end();
			it++){
		BSPos[it->first] = it->second;
	}

	// Mobile stations should be randomly rotated..?
	double ***MsAntennaPosition = new double**[numberOfMobileStations]; /*!< Position vector of Receiver antenna */
    	for(int i = 0; i < numberOfMobileStations; i++){
		MsAntennaPosition[i] = new double*[NumMsAntenna];
		for(int j = 0; j < (NumMsAntenna/2); j++){
			MsAntennaPosition[i][j] = new double[3];
			MsAntennaPosition[i][j][0] = MSPos[i].x - (0.25 * wavelength * (NumMsAntenna - 1 - (j*2.0)));
			MsAntennaPosition[i][j][1] = MSPos[i].y;
			MsAntennaPosition[i][j][2] = heightUE;
		}
		for(int j = (NumMsAntenna/2); j < NumMsAntenna ; j++){
			MsAntennaPosition[i][j] = new double[3];
			MsAntennaPosition[i][j][0] = MsAntennaPosition[i][j-1][0] + (0.5 * wavelength);
			MsAntennaPosition[i][j][1] = MSPos[i].y;
			MsAntennaPosition[i][j][2] = heightUE;
		}
	}

	vector<vector<double>> AoA_LOS_dir;
	vector<vector<double>> ZoA_LOS_dir;
	vector<vector<double>> AoD_LOS_dir;
	vector<vector<double>> ZoD_LOS_dir;
	std::tie(AoA_LOS_dir,ZoA_LOS_dir,AoD_LOS_dir,ZoD_LOS_dir) = recomputeAngleDirection(
			MSPos,
			BSPos,
			heightBS,
			heightUE
			);
	
	//Assign LOS Conditions:
	vector<vector<bool>> LOSCondition(genLosCond(BSPos,MSPos));
	
	vector<vector<double>> sigma_ds_LOS(MSPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_asD_LOS(MSPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_asA_LOS(MSPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_zsD_LOS(MSPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_zsA_LOS(MSPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_sf_LOS(MSPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_kf_LOS(MSPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_ds_NLOS(MSPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_asD_NLOS(MSPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_asA_NLOS(MSPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_zsD_NLOS(MSPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_zsA_NLOS(MSPos.size(),
			vector<double>(neighbourPositions.size()));
	vector<vector<double>> sigma_sf_NLOS(MSPos.size(),
			vector<double>(neighbourPositions.size()));
	recomputeLargeScaleParameters(BSPos,MSPos,
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
	
	// Main Loop:
	// TODO: allow Elevation Angle
	// TODO: adjust for METIS
	// TODO: add polarization case
	//------------------------------------------------------------------
		
	// Compute all constant values only once:
	
	// Wavenumber k_0
	double k_0 = (2 * pi) / (speedOfLight / freq_c);
	
	// Create all local variables:
	complex<double> doppler, pol;
	double sq_P_over_M;

	// For LOS case:
	complex<double> *****raySum_LOS = new complex<double>****[numberOfMobileStations];
	int clusterIdx_LOS;
	
	for(int m = 0; m < numberOfMobileStations; m++){
		// +4 Clusters, because 2 are subdivided into 6, resulting in 4 more.
		raySum_LOS[m] = new complex<double>***[(N_cluster_LOS + 4)];
		for(int i = 0; i < (N_cluster_LOS + 4); i++){
			raySum_LOS[m][i] = new complex<double>**[timeSamples];
			for(int j = 0; j < timeSamples; j++){
				raySum_LOS[m][i][j] = new complex<double>*[NumMsAntenna];
				for(int k = 0; k < NumMsAntenna; k++){
					raySum_LOS[m][i][j][k] = new complex<double>[NumBsAntenna]();
				}
			}
		}
	}
	// For NLOS case:
	complex<double> *****raySum = new complex<double>****[numberOfMobileStations];
	int clusterIdx;
	for(int m = 0; m < numberOfMobileStations; m++){
		// +4 Clusters, because 2 are subdivided into 6, resulting in 4 more.
		raySum[m] = new complex<double>***[(N_cluster_NLOS + 4)];
			for(int i = 0; i < (N_cluster_NLOS + 4); i++){
			raySum[m][i] = new complex<double>**[timeSamples];
			for(int j = 0; j < timeSamples; j++){
				raySum[m][i][j] = new complex<double>*[NumMsAntenna];
				for(int k = 0; k < NumMsAntenna; k++){
					raySum[m][i][j][k] = new complex<double>[NumBsAntenna]();
				}
			}
		}
	}

	// initialize arrays for interferer ray sums
	complex<double> ******raySumInterferer_LOS;
	complex<double> ******raySumInterferer;
	if(numOfInterferers>0){
		// Only compute values for interferers if there actually are any 
		// interfering neighbours
		raySumInterferer_LOS = new complex<double>*****[numberOfMobileStations];
		for (int m = 0; m < numberOfMobileStations; m++){
			raySumInterferer_LOS[m] = new complex<double>****[numOfInterferers];
			for(int i = 0; i < numOfInterferers; i++){
				raySumInterferer_LOS[m][i] = new complex<double>***[(N_cluster_LOS + 4)];
				for(int j = 0; j < (N_cluster_LOS + 4); j++){
					raySumInterferer_LOS[m][i][j] = new complex<double>**[timeSamples];
					for(int k = 0; k < timeSamples; k++){
						raySumInterferer_LOS[m][i][j][k] = new complex<double>*[NumMsAntenna];
						for(int l = 0; l < NumMsAntenna; l++){
							raySumInterferer_LOS[m][i][j][k][l] = new complex<double>[NumBsAntenna]();
						}
					}
				}
			}
		}

		raySumInterferer = new complex<double>*****[numberOfMobileStations];
		for(int m = 0; m < numberOfMobileStations; m++){
			// +4 Clusters, because 2 are subdivided into 6, resulting in 4 more.
			raySumInterferer[m] = new complex<double>****[numOfInterferers];
			for(int i = 0; i < numOfInterferers; i++){
				raySumInterferer[m][i] = new complex<double>***[(N_cluster_NLOS + 4)];
				for(int j = 0; j < (N_cluster_NLOS + 4); j++){
					raySumInterferer[m][i][j] = new complex<double>**[timeSamples];
					for(int k = 0; k < timeSamples; k++){
						raySumInterferer[m][i][j][k] = new complex<double>*[NumMsAntenna];
						for(int l = 0; l < NumMsAntenna; l++){
							raySumInterferer[m][i][j][k][l] = new complex<double>[NumBsAntenna]();
						}
					}
				}
			}
		}
	}

	std::cout << "START MAIN LOOP for BS: " << bsId << std::endl;
	
	double AoA[3];
	double AoD[3];
	double MSgain, BSgain;
	int SC_1[10] = {0,1,2,3,4,5,6,7,18,19}; //rays for sub-cluster 1 of first two clusters
	int SC_2[6] = {8,9,10,11,16,17}; //rays for sub-cluster 2 of first two clusters
	int SC_3[4] = {12,13,14,15}; //rays for sub-cluster 3 of first two clusters
	int size_SC_1 = sizeof(SC_1)/sizeof(*SC_1); //length of SC_1
	int size_SC_2 = sizeof(SC_2)/sizeof(*SC_2); //length of SC_2
	int size_SC_3 = sizeof(SC_3)/sizeof(*SC_3); //length of SC_3
	
	// Cycle through all MS, TODO: change the loop iteration for sub clusters as well as for proper link to link implementation
	for(int i = 0; i < numberOfMobileStations; i++){
		// Cycle through all interferer base stations
		int idIdx = 1;
		
		// Convert Cartesian Direction to spherical azimuth angle
		double AoMD = atan((MSVelDir[i][1] / MSVelDir[i][0]));
		
		for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
			if(it->first == bsId){
				// NLOS CASE:
				if(!LOSCondition[i][0]){
					// Cycle through all Receiver antennas (MS)
						std::cout << "NLOS" << std::endl;
					for(int u = 0; u < NumMsAntenna; u++){
						// Cycle through all Transmitter antennas (BS)
						for(int s = 0; s < NumBsAntenna; s++){
							// Cycle through all Paths/Clusters
							clusterIdx = 0;
							for(int n = 0; n < N_cluster_NLOS; n++){
								int L;
								double P_n[3]; // Oversized, but 2 doubles really doesnt matter
								// Two strongest clusters are divided into 3 subclusters!
								//std::cout << "C_Power: " << clusterPowers[i][0][n] << std::endl;
								if(n < 2){
									L = 3;
									P_n[0] = 10.0/20.0 * clusterPowers[i][0][n];
									P_n[1] =  6.0/20.0 * clusterPowers[i][0][n];
									P_n[2] =  4.0/20.0 * clusterPowers[i][0][n];
								}else{
									L = 1;
									P_n[0] = clusterPowers[i][0][n];
								}
								// Cycle through all subclusters (only a loop for cluster 1 and 2)
								for(int l = 0; l < L; l++){
									if (L == 3){
										if (l == 0){ // for sub-cluster 1 of first two clusters
											sq_P_over_M = sqrt(P_n[l] / size_SC_1);
											//std::cout << "P_n[l]: " << P_n[l] << std::endl;
											//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
											// Cycle through all NLOS Ray
											for(int m = 0; m < size_SC_1; m++){		
												AoA[0] = sin(elevation_ASA[i][0][n][SC_1[m]]*pi/180) * cos(azimuth_ASA[i][0][n][SC_1[m]]*pi/180);
												AoA[1] = sin(elevation_ASA[i][0][n][SC_1[m]]*pi/180) * sin(azimuth_ASA[i][0][n][SC_1[m]]*pi/180);
												AoA[2] = cos(elevation_ASA[i][0][n][SC_1[m]]*pi/180);
										
												AoD[0] = sin(elevation_ASD[i][0][n][SC_1[m]]*pi/180) * cos(azimuth_ASD[i][0][n][SC_1[m]]*pi/180);
												AoD[1] = sin(elevation_ASD[i][0][n][SC_1[m]]*pi/180) * sin(azimuth_ASD[i][0][n][SC_1[m]]*pi/180);
												AoD[2] = cos(elevation_ASD[i][0][n][SC_1[m]]*pi/180);
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2] )) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
												// TODO: Include Polarization
												// Replacement for polarization. (Below 4.14)
												//pol = exp( complex<double>(0,uniform(0,2*pi)));
												MSgain = getMSGain(azimuth_ASA[i][0][n][SC_1[m]]*pi/180, elevation_ASA[i][0][n][SC_1[m]]*pi/180);
												BSgain = getBSGain(azimuth_ASD[i][0][n][SC_1[m]]*pi/180, elevation_ASD[i][0][n][SC_1[m]]*pi/180);
												pol = MSgain * BSgain * exp(complex<double>(0, randomPhase[i][0][n][SC_1[m]][0]));
										
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
													//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
													doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(azimuth_ASA[i][0][n][SC_1[m]]*pi/180 - AoMD) * timeVector[i][t] ) );
													//std::cout << "Doppler: " << doppler << std::endl;
													raySum[i][clusterIdx][t][u][s] += pol * doppler * exp_arrival * exp_departure;
													if(m + 1 == size_SC_1){
														raySum[i][clusterIdx][t][u][s] = sq_P_over_M * raySum[i][clusterIdx][t][u][s];
													}
												} // End time axis							
											} // End cycle Rays
											clusterIdx++;
										} else if (l == 1){ // for sub-cluster 2 of first two clusters
											sq_P_over_M = sqrt(P_n[l] / size_SC_2);
											// Cycle through all NLOS Rays
											for(int m = 0; m < size_SC_2; m++){		
												AoA[0] = sin(elevation_ASA[i][0][n][SC_2[m]]*pi/180) * cos(azimuth_ASA[i][0][n][SC_2[m]]*pi/180);
												AoA[1] = sin(elevation_ASA[i][0][n][SC_2[m]]*pi/180) * sin(azimuth_ASA[i][0][n][SC_2[m]]*pi/180);
												AoA[2] = cos(elevation_ASA[i][0][n][SC_2[m]]*pi/180);
										
												AoD[0] = sin(elevation_ASD[i][0][n][SC_2[m]]*pi/180) * cos(azimuth_ASD[i][0][n][SC_2[m]]*pi/180);
												AoD[1] = sin(elevation_ASD[i][0][n][SC_2[m]]*pi/180) * sin(azimuth_ASD[i][0][n][SC_2[m]]*pi/180);
												AoD[2] = cos(elevation_ASD[i][0][n][SC_2[m]]*pi/180);
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
												// TODO: Include Polarization
												// Replacement for polarization. (Below 4.14)
												//pol = exp( complex<double>(0,uniform(0,2*pi)));
												MSgain = getMSGain(azimuth_ASA[i][0][n][SC_2[m]]*pi/180, elevation_ASA[i][0][n][SC_2[m]]*pi/180);
												BSgain = getBSGain(azimuth_ASD[i][0][n][SC_2[m]]*pi/180, elevation_ASD[i][0][n][SC_2[m]]*pi/180);
												pol = MSgain * BSgain * exp(complex<double>(0, randomPhase[i][0][n][SC_2[m]][0]));
										
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
													//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
													doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(azimuth_ASA[i][0][n][SC_2[m]]*pi/180 - AoMD) * timeVector[i][t] ) );
													//std::cout << "Doppler: " << doppler << std::endl;
													raySum[i][clusterIdx][t][u][s] += pol * doppler * exp_arrival * exp_departure;
													if(m + 1 == size_SC_2){
														//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
														raySum[i][clusterIdx][t][u][s] = sq_P_over_M * raySum[i][clusterIdx][t][u][s];
													}
												} // End time axis							
											} // End cycle Rays
											clusterIdx++;
										} else { //for sub-cluster 3 of first two clusters
											sq_P_over_M = sqrt(P_n[l] / size_SC_3);
											// Cycle through all NLOS Ray
											for(int m = 0; m < size_SC_3; m++){		
												AoA[0] = sin(elevation_ASA[i][0][n][SC_3[m]]*pi/180) * cos(azimuth_ASA[i][0][n][SC_3[m]]*pi/180);
												AoA[1] = sin(elevation_ASA[i][0][n][SC_3[m]]*pi/180) * sin(azimuth_ASA[i][0][n][SC_3[m]]*pi/180);
												AoA[2] = cos(elevation_ASA[i][0][n][SC_3[m]]*pi/180);
										
												AoD[0] = sin(elevation_ASD[i][0][n][SC_3[m]]*pi/180) * cos(azimuth_ASD[i][0][n][SC_3[m]]*pi/180);
												AoD[1] = sin(elevation_ASD[i][0][n][SC_3[m]]*pi/180) * sin(azimuth_ASD[i][0][n][SC_3[m]]*pi/180);
												AoD[2] = cos(elevation_ASD[i][0][n][SC_3[m]]*pi/180);
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
												// TODO: Include Polarization
												// Replacement for polarization. (Below 4.14)
												//pol = exp( complex<double>(0,uniform(0,2*pi)));
												MSgain = getMSGain(azimuth_ASA[i][0][n][SC_3[m]]*pi/180, elevation_ASA[i][0][n][SC_3[m]]*pi/180);
												BSgain = getBSGain(azimuth_ASD[i][0][n][SC_3[m]]*pi/180, elevation_ASD[i][0][n][SC_3[m]]*pi/180);
												pol = MSgain * BSgain * exp(complex<double>(0, randomPhase[i][0][n][SC_3[m]][0]));
										
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
													//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
													doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(azimuth_ASA[i][0][n][SC_3[m]]*pi/180 - AoMD) * timeVector[i][t] ) );
													//std::cout << "Doppler: " << doppler << std::endl;
													raySum[i][clusterIdx][t][u][s] += pol * doppler * exp_arrival * exp_departure;
													if(m + 1 == size_SC_3){
														//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
														raySum[i][clusterIdx][t][u][s] = sq_P_over_M * raySum[i][clusterIdx][t][u][s];
													}
												} // End time axis							
											} // End cycle Rays
											clusterIdx++;
										} //end if (for all sub-clusters)
									}else{ // for all clusters from N = 3 onwards
										sq_P_over_M = sqrt(P_n[l] / numOfRays_NLOS);
										//std::cout << "P_n[i]: " << P_n[i] << std::endl;
										//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
										// Cycle through all NLOS Rays
										for(int m = 0; m < numOfRays_NLOS; m++){		
											AoA[0] = sin(elevation_ASA[i][0][n][m]*pi/180) * cos(azimuth_ASA[i][0][n][m]*pi/180);
											AoA[1] = sin(elevation_ASA[i][0][n][m]*pi/180) * sin(azimuth_ASA[i][0][n][m]*pi/180);
											AoA[2] = cos(elevation_ASA[i][0][n][m]*pi/180);
										
											AoD[0] = sin(elevation_ASD[i][0][n][m]*pi/180) * cos(azimuth_ASD[i][0][n][m]*pi/180);
											AoD[1] = sin(elevation_ASD[i][0][n][m]*pi/180) * sin(azimuth_ASD[i][0][n][m]*pi/180);
											AoD[2] = cos(elevation_ASD[i][0][n][m]*pi/180);
										
											complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2])) );
											complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
											// TODO: Include Polarization
											// Replacement for polarization. (Below 4.14)
											//pol = exp( complex<double>(0,uniform(0,2*pi)));
											MSgain = getMSGain(azimuth_ASA[i][0][n][m]*pi/180, elevation_ASA[i][0][n][m]*pi/180);
											BSgain = getBSGain(azimuth_ASD[i][0][n][m]*pi/180, elevation_ASD[i][0][n][m]*pi/180);
											pol = MSgain * BSgain * exp(complex<double>(0, randomPhase[i][0][n][m][0]));
										
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
												//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
												doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(azimuth_ASA[i][0][n][m]*pi/180 - AoMD) * timeVector[i][t] ) );
												//std::cout << "Doppler: " << doppler << std::endl;
												raySum[i][clusterIdx][t][u][s] += pol * doppler * exp_arrival * exp_departure;
												if(m + 1 == numOfRays_NLOS){
													//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
													raySum[i][clusterIdx][t][u][s] = sq_P_over_M * raySum[i][clusterIdx][t][u][s];
												}
											} // End time axis							
										} // End cycle Rays
										clusterIdx++;
									} //End if
								} // End cycle Subclusters
							} // End cycle Clusters
						} // End BS antenna
					} // End MS antenna
				}else{
					// LOS case, TODO: correct the computation of raySum_LOS
					// Cycle through all Receiver antennas (MS)
					std::cout<< "LOS" << std::endl;
					for(int u = 0; u < NumMsAntenna; u++){
						// Cycle through all Transmitter antennas (BS)
						for(int s = 0; s < NumBsAntenna; s++){
							// Cycle through all Paths/Clusters
							clusterIdx_LOS = 0;
							for(int n = 0; n < N_cluster_LOS; n++){
								int L;
								double P_n[3]; // Oversized, but 2 doubles really doesnt matter
								// Two strongest clusters are divided into 3 subclusters!
								//std::cout << "C_Power: " << clusterPowers[i][0][n] << std::endl;
								if(n < 2){
									L = 3;
									P_n[0] = 10.0/20.0 * clusterPowers[i][0][n];
									P_n[1] =  6.0/20.0 * clusterPowers[i][0][n];
									P_n[2] =  4.0/20.0 * clusterPowers[i][0][n];
								}else{
									L = 1;
									P_n[0] = clusterPowers[i][0][n];
								}
								// Cycle through all subclusters (only a loop for cluster 1 and 2)
								for(int l = 0; l < L; l++){
									if (L == 3){
										if (l == 0){ // for sub-cluster 1 of first two clusters
											sq_P_over_M = sqrt(P_n[l] / size_SC_1);
											//std::cout << "P_n[l]: " << P_n[l] << std::endl;
											//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
											// Cycle through all LOS Rays
											for(int m = 0; m < size_SC_1; m++){		
												AoA[0] = sin(elevation_ASA[i][0][n][SC_1[m]]*pi/180) * cos(azimuth_ASA[i][0][n][SC_1[m]]*pi/180);
												AoA[1] = sin(elevation_ASA[i][0][n][SC_1[m]]*pi/180) * sin(azimuth_ASA[i][0][n][SC_1[m]]*pi/180);
												AoA[2] = cos(elevation_ASA[i][0][n][SC_1[m]]*pi/180);
										
												AoD[0] = sin(elevation_ASD[i][0][n][SC_1[m]]*pi/180) * cos(azimuth_ASD[i][0][n][SC_1[m]]*pi/180);
												AoD[1] = sin(elevation_ASD[i][0][n][SC_1[m]]*pi/180) * sin(azimuth_ASD[i][0][n][SC_1[m]]*pi/180);
												AoD[2] = cos(elevation_ASD[i][0][n][SC_1[m]]*pi/180);
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
												// TODO: Include Polarization
												// Replacement for polarization. (Below 4.14)
												//pol = exp( complex<double>(0,uniform(0,2*pi)));
												MSgain = getMSGain(azimuth_ASA[i][0][n][SC_1[m]]*pi/180, elevation_ASA[i][0][n][SC_1[m]]*pi/180);
												BSgain = getBSGain(azimuth_ASD[i][0][n][SC_1[m]]*pi/180, elevation_ASD[i][0][n][SC_1[m]]*pi/180);
												pol = MSgain * BSgain * exp(complex<double>(0, randomPhase[i][0][n][SC_1[m]][0]));
										
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
													//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
													doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(azimuth_ASA[i][0][n][SC_1[m]]*pi/180 - AoMD) * timeVector[i][t] ) );
													//std::cout << "Doppler: " << doppler << std::endl;
													raySum_LOS[i][clusterIdx_LOS][t][u][s] += pol * doppler * exp_arrival * exp_departure;
													if(m + 1 == size_SC_1){
														//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
														double K = sigma_kf_LOS[i][0]; 
														raySum_LOS[i][clusterIdx_LOS][t][u][s] = (sqrt(1/(K + 1))) * sq_P_over_M * raySum_LOS[i][clusterIdx_LOS][t][u][s];
													}
												} // End time axis							
											} // End cycle Rays
											clusterIdx_LOS++;
										} else if (l == 1){ // for sub-cluster 2 of first two clusters
											sq_P_over_M = sqrt(P_n[l] / size_SC_2);
											//std::cout << "P_n[l]: " << P_n[l] << std::endl;
											//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
											// Cycle through all LOS Rays
											for(int m = 0; m < size_SC_2; m++){		
												AoA[0] = sin(elevation_ASA[i][0][n][SC_2[m]]*pi/180) * cos(azimuth_ASA[i][0][n][SC_2[m]]*pi/180);
												AoA[1] = sin(elevation_ASA[i][0][n][SC_2[m]]*pi/180) * sin(azimuth_ASA[i][0][n][SC_2[m]]*pi/180);
												AoA[2] = cos(elevation_ASA[i][0][n][SC_2[m]]*pi/180);
										
												AoD[0] = sin(elevation_ASD[i][0][n][SC_2[m]]*pi/180) * cos(azimuth_ASD[i][0][n][SC_2[m]]*pi/180);
												AoD[1] = sin(elevation_ASD[i][0][n][SC_2[m]]*pi/180) * sin(azimuth_ASD[i][0][n][SC_2[m]]*pi/180);
												AoD[2] = cos(elevation_ASD[i][0][n][SC_2[m]]*pi/180);
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
												// TODO: Include Polarization
												// Replacement for polarization. (Below 4.14)
												//pol = exp( complex<double>(0,uniform(0,2*pi)));
												MSgain = getMSGain(azimuth_ASA[i][0][n][SC_2[m]]*pi/180, elevation_ASA[i][0][n][SC_2[m]]*pi/180);
												BSgain = getBSGain(azimuth_ASD[i][0][n][SC_2[m]]*pi/180, elevation_ASD[i][0][n][SC_2[m]]*pi/180);
												pol = MSgain * BSgain * exp(complex<double>(0, randomPhase[i][0][n][SC_2[m]][0]));
										
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
													//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
													doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(azimuth_ASA[i][0][n][SC_2[m]]*pi/180 - AoMD) * timeVector[i][t] ) );
													//std::cout << "Doppler: " << doppler << std::endl;
													raySum_LOS[i][clusterIdx_LOS][t][u][s] += pol * doppler * exp_arrival * exp_departure;
													if(m + 1 == size_SC_2){
														//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
														double K = sigma_kf_LOS[i][0]; 
														raySum_LOS[i][clusterIdx_LOS][t][u][s] = (sqrt(1/(K + 1))) * sq_P_over_M * raySum_LOS[i][clusterIdx_LOS][t][u][s];
													}
												} // End time axis							
											} // End cycle Rays
											clusterIdx_LOS++;
										} else { //for sub-cluster 3 of first two clusters
											sq_P_over_M = sqrt(P_n[l] / size_SC_3);
											//std::cout << "P_n[l]: " << P_n[l] << std::endl;
											//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
											// Cycle through all LOS Rays
											for(int m = 0; m < size_SC_3; m++){		
												AoA[0] = sin(elevation_ASA[i][0][n][SC_3[m]]*pi/180) * cos(azimuth_ASA[i][0][n][SC_3[m]]*pi/180);
												AoA[1] = sin(elevation_ASA[i][0][n][SC_3[m]]*pi/180) * sin(azimuth_ASA[i][0][n][SC_3[m]]*pi/180);
												AoA[2] = cos(elevation_ASA[i][0][n][SC_3[m]]*pi/180);
										
												AoD[0] = sin(elevation_ASD[i][0][n][SC_3[m]]*pi/180) * cos(azimuth_ASD[i][0][n][SC_3[m]]*pi/180);
												AoD[1] = sin(elevation_ASD[i][0][n][SC_3[m]]*pi/180) * sin(azimuth_ASD[i][0][n][SC_3[m]]*pi/180);
												AoD[2] = cos(elevation_ASD[i][0][n][SC_3[m]]*pi/180);
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
												// TODO: Include Polarization
												// Replacement for polarization. (Below 4.14)
												//pol = exp( complex<double>(0,uniform(0,2*pi)));
												MSgain = getMSGain(azimuth_ASA[i][0][n][SC_3[m]]*pi/180, elevation_ASA[i][0][n][SC_3[m]]*pi/180);
												BSgain = getBSGain(azimuth_ASD[i][0][n][SC_3[m]]*pi/180, elevation_ASD[i][0][n][SC_3[m]]*pi/180);
												pol = MSgain * BSgain * exp(complex<double>(0, randomPhase[i][0][n][SC_3[m]][0]));
										
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
													//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
													doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(azimuth_ASA[i][0][n][SC_3[m]]*pi/180 - AoMD) * timeVector[i][t] ) );
													//std::cout << "Doppler: " << doppler << std::endl;
													raySum_LOS[i][clusterIdx_LOS][t][u][s] += pol * doppler * exp_arrival * exp_departure;
													if(m + 1 == size_SC_3){
														//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
														double K = sigma_kf_LOS[i][0]; 
														raySum_LOS[i][clusterIdx_LOS][t][u][s] = (sqrt(1/(K + 1))) * sq_P_over_M * raySum_LOS[i][clusterIdx_LOS][t][u][s];
													}
												} // End time axis							
											} // End cycle Rays
											clusterIdx_LOS++;
										} //end if (for all sub-clusters)
									}else{ // for all clusters from N = 3 onwards
										sq_P_over_M = sqrt(P_n[l] / numOfRays_LOS);
										//std::cout << "P_n[l]: " << P_n[l] << std::endl;
										//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
										// Cycle through all LOS Rays
										for(int m = 0; m < numOfRays_LOS; m++){		
											AoA[0] = sin(elevation_ASA[i][0][n][m]*pi/180) * cos(azimuth_ASA[i][0][n][m]*pi/180);
											AoA[1] = sin(elevation_ASA[i][0][n][m]*pi/180) * sin(azimuth_ASA[i][0][n][m]*pi/180);
											AoA[2] = cos(elevation_ASA[i][0][n][m]*pi/180);
										
											AoD[0] = sin(elevation_ASD[i][0][n][m]*pi/180) * cos(azimuth_ASD[i][0][n][m]*pi/180);
											AoD[1] = sin(elevation_ASD[i][0][n][m]*pi/180) * sin(azimuth_ASD[i][0][n][m]*pi/180);
											AoD[2] = cos(elevation_ASD[i][0][n][m]*pi/180);
										
											complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2])) );
											complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
											// TODO: Include Polarization
											// Replacement for polarization. (Below 4.14)
											//pol = exp( complex<double>(0,uniform(0,2*pi)));
											MSgain = getMSGain(azimuth_ASA[i][0][n][m]*pi/180, elevation_ASA[i][0][n][m]*pi/180);
											BSgain = getBSGain(azimuth_ASD[i][0][n][m]*pi/180, elevation_ASD[i][0][n][m]*pi/180);
											pol = MSgain * BSgain * exp(complex<double>(0, randomPhase[i][0][n][m][0]));
										
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
												//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
												doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(azimuth_ASA[i][0][n][m]*pi/180 - AoMD) * timeVector[i][t] ) );
												//std::cout << "Doppler: " << doppler << std::endl;
												raySum_LOS[i][clusterIdx_LOS][t][u][s] += pol * doppler * exp_arrival * exp_departure;
												if(m + 1 == numOfRays_LOS){
													//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
													double K = sigma_kf_LOS[i][0]; 
													raySum_LOS[i][clusterIdx_LOS][t][u][s] = (sqrt(1/(K + 1))) * sq_P_over_M * raySum_LOS[i][clusterIdx_LOS][t][u][s];
												}
											} // End time axis							
										} // End cycle Rays
										clusterIdx_LOS++;
									}//end if
								} // End cycle Subclusters
								if(n == 0){ // for adding the additional LOS component, according to formula 7-61 in METIS 1.2
									std::cout << "n: " << n << " LOS Computation for MS " << i << " and BS 0 " << std::endl;
									double K = sigma_kf_LOS[i][0]; 						// K-factor in linear scale 
									AoA[0] = sin(ZoA_LOS_dir[i][0]) * cos(AoA_LOS_dir[i][0]);
									AoA[1] = sin(ZoA_LOS_dir[i][0]) * sin(AoA_LOS_dir[i][0]);
									AoA[2] = cos(ZoA_LOS_dir[i][0]);
										
									AoD[0] = sin(ZoD_LOS_dir[i][0]) * cos(AoD_LOS_dir[i][0]);
									AoD[1] = sin(ZoD_LOS_dir[i][0]) * sin(AoD_LOS_dir[i][0]);
									AoD[2] = cos(ZoD_LOS_dir[i][0]);
										
									complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2])) );
									complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
									// TODO: Include Polarization in generic form
									// Replacement for polarization. (Below 4.14)
									//pol = exp( complex<double>(0,uniform(0,2*pi)));
									MSgain = getMSGain(AoA_LOS_dir[i][0], ZoA_LOS_dir[i][0]);
									BSgain = getBSGain(AoD_LOS_dir[i][0], ZoD_LOS_dir[i][0]);
									pol = MSgain * BSgain * exp(complex<double>(0, randomPhase_LOS[i][0]));
										
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
										//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
										doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(AoA_LOS_dir[i][0] - AoMD) * timeVector[i][t] ) );
										//std::cout << "Doppler: " << doppler << std::endl;
										raySum_LOS[i][0][t][u][s] += (sqrt(K / (K + 1))) * pol * doppler * exp_arrival * exp_departure;
									} // End time axis
								}
							} // End cycle Clusters
						} // End BS antenna
					} // End MS antenna
				} //End if for own BS
			}else{
				// calculate interferer
				// NLOS CASE:
				if(!LOSCondition[i][idIdx]){
					// Cycle through all Receiver antennas (MS)
					for(int u = 0; u < NumMsAntenna; u++){
						// Cycle through all Transmitter antennas (BS)
						for(int s = 0; s < NumBsAntenna; s++){
							// Cycle through all Paths/Clusters
							clusterIdx = 0;
							for(int n = 0; n < N_cluster_NLOS; n++){
								int L;
								double P_n[3]; // Oversized, but 2 doubles really doesnt matter
								// Two strongest clusters are divided into 3 subclusters!
								if(n < 2){
									L = 3;
									P_n[0] = 10.0/20.0 * clusterPowers[i][idIdx][n];
									P_n[1] =  6.0/20.0 * clusterPowers[i][idIdx][n];
									P_n[2] =  4.0/20.0 * clusterPowers[i][idIdx][n];
								}else{
									L = 1;
									P_n[0] = clusterPowers[i][idIdx][n];
								}
								// Cycle through all subclusters (only a loop for cluster 1 and 2)
								for(int l = 0; l < L; l++){
									if (L == 3){
										if (l == 0){ // for sub-cluster 1 of first two clusters
											sq_P_over_M = sqrt(P_n[l] / size_SC_1);
											//std::cout << "P_n[l]: " << P_n[l] << std::endl;
											//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
											// Cycle through all NLOS Ray
											for(int m = 0; m < size_SC_1; m++){		
												AoA[0] = sin(elevation_ASA[i][idIdx][n][SC_1[m]]*pi/180) * cos(azimuth_ASA[i][idIdx][n][SC_1[m]]*pi/180);
												AoA[1] = sin(elevation_ASA[i][idIdx][n][SC_1[m]]*pi/180) * sin(azimuth_ASA[i][idIdx][n][SC_1[m]]*pi/180);
												AoA[2] = cos(elevation_ASA[i][idIdx][n][SC_1[m]]*pi/180);
										
												AoD[0] = sin(elevation_ASD[i][idIdx][n][SC_1[m]]*pi/180) * cos(azimuth_ASD[i][idIdx][n][SC_1[m]]*pi/180);
												AoD[1] = sin(elevation_ASD[i][idIdx][n][SC_1[m]]*pi/180) * sin(azimuth_ASD[i][idIdx][n][SC_1[m]]*pi/180);
												AoD[2] = cos(elevation_ASD[i][idIdx][n][SC_1[m]]*pi/180);
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
												// TODO: Include Polarization
												// Replacement for polarization. (Below 4.14)
												//pol = exp( complex<double>(0,uniform(0,2*pi)));
												MSgain = getMSGain(azimuth_ASA[i][idIdx][n][SC_1[m]]*pi/180, elevation_ASA[i][idIdx][n][SC_1[m]]*pi/180);
												BSgain = getBSGain(azimuth_ASD[i][idIdx][n][SC_1[m]]*pi/180, elevation_ASD[i][idIdx][n][SC_1[m]]*pi/180);
												pol = MSgain * BSgain * exp(complex<double>(0, randomPhase[i][idIdx][n][SC_1[m]][0]));
										
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
													//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
													doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(azimuth_ASA[i][idIdx][n][SC_1[m]]*pi/180 - AoMD) * timeVector[i][t] ) );
													//std::cout << "Doppler: " << doppler << std::endl;
													raySumInterferer[i][idIdx-1][clusterIdx][t][u][s] += pol * doppler * exp_arrival * exp_departure;
													if(m + 1 == size_SC_1){
														//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
														raySumInterferer[i][idIdx-1][clusterIdx][t][u][s] = sq_P_over_M * raySumInterferer[i][idIdx-1][clusterIdx][t][u][s];
														//std::cout << "raySumInterferer[i][clusterIdx][t][u][s]: " << raySumInterferer[i][clusterIdx][t][u][s] << std::endl;
													}
												} // End time axis							
											} // End cycle Rays
											clusterIdx++;
										} else if (l == 1){ // for sub-cluster 2 of first two clusters
											sq_P_over_M = sqrt(P_n[l] / size_SC_2);
											//std::cout << "P_n[l]: " << P_n[l] << std::endl;
											//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
											// Cycle through all NLOS Ray
											for(int m = 0; m < size_SC_2; m++){		
												AoA[0] = sin(elevation_ASA[i][idIdx][n][SC_2[m]]*pi/180) * cos(azimuth_ASA[i][idIdx][n][SC_2[m]]*pi/180);
												AoA[1] = sin(elevation_ASA[i][idIdx][n][SC_2[m]]*pi/180) * sin(azimuth_ASA[i][idIdx][n][SC_2[m]]*pi/180);
												AoA[2] = cos(elevation_ASA[i][idIdx][n][SC_2[m]]*pi/180);
										
												AoD[0] = sin(elevation_ASD[i][idIdx][n][SC_2[m]]*pi/180) * cos(azimuth_ASD[i][idIdx][n][SC_2[m]]*pi/180);
												AoD[1] = sin(elevation_ASD[i][idIdx][n][SC_2[m]]*pi/180) * sin(azimuth_ASD[i][idIdx][n][SC_2[m]]*pi/180);
												AoD[2] = cos(elevation_ASD[i][idIdx][n][SC_2[m]]*pi/180);
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
												// TODO: Include Polarization
												// Replacement for polarization. (Below 4.14)
												//pol = exp( complex<double>(0,uniform(0,2*pi)));
												MSgain = getMSGain(azimuth_ASA[i][idIdx][n][SC_2[m]]*pi/180, elevation_ASA[i][idIdx][n][SC_2[m]]*pi/180);
												BSgain = getBSGain(azimuth_ASD[i][idIdx][n][SC_2[m]]*pi/180, elevation_ASD[i][idIdx][n][SC_2[m]]*pi/180);
												pol = MSgain * BSgain * exp(complex<double>(0, randomPhase[i][idIdx][n][SC_2[m]][0]));
										
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
													//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
													doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(azimuth_ASA[i][idIdx][n][SC_2[m]]*pi/180 - AoMD) * timeVector[i][t] ) );
													//std::cout << "Doppler: " << doppler << std::endl;
													raySumInterferer[i][idIdx-1][clusterIdx][t][u][s] += pol * doppler * exp_arrival * exp_departure;
													if(m + 1 == size_SC_2){
														//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
														raySumInterferer[i][idIdx-1][clusterIdx][t][u][s] = sq_P_over_M * raySumInterferer[i][idIdx-1][clusterIdx][t][u][s];
														//std::cout << "raySumInterferer[i][clusterIdx][t][u][s]: " << raySumInterferer[i][clusterIdx][t][u][s] << std::endl;
													}
												} // End time axis							
											} // End cycle Rays
											clusterIdx++;
										} else { //for sub-cluster 3 of first two clusters
											sq_P_over_M = sqrt(P_n[l] / size_SC_3);
											//std::cout << "P_n[l]: " << P_n[l] << std::endl;
											//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
											// Cycle through all NLOS Ray
											for(int m = 0; m < size_SC_3; m++){		
												AoA[0] = sin(elevation_ASA[i][idIdx][n][SC_3[m]]*pi/180) * cos(azimuth_ASA[i][idIdx][n][SC_3[m]]*pi/180);
												AoA[1] = sin(elevation_ASA[i][idIdx][n][SC_3[m]]*pi/180) * sin(azimuth_ASA[i][idIdx][n][SC_3[m]]*pi/180);
												AoA[2] = cos(elevation_ASA[i][idIdx][n][SC_3[m]]*pi/180);
										
												AoD[0] = sin(elevation_ASD[i][idIdx][n][SC_3[m]]*pi/180) * cos(azimuth_ASD[i][idIdx][n][SC_3[m]]*pi/180);
												AoD[1] = sin(elevation_ASD[i][idIdx][n][SC_3[m]]*pi/180) * sin(azimuth_ASD[i][idIdx][n][SC_3[m]]*pi/180);
												AoD[2] = cos(elevation_ASD[i][idIdx][n][SC_3[m]]*pi/180);
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
												// TODO: Include Polarization
												// Replacement for polarization. (Below 4.14)
												//pol = exp( complex<double>(0,uniform(0,2*pi)));
												MSgain = getMSGain(azimuth_ASA[i][idIdx][n][SC_3[m]]*pi/180, elevation_ASA[i][idIdx][n][SC_3[m]]*pi/180);
												BSgain = getBSGain(azimuth_ASD[i][idIdx][n][SC_3[m]]*pi/180, elevation_ASD[i][idIdx][n][SC_3[m]]*pi/180);
												pol = MSgain * BSgain * exp(complex<double>(0, randomPhase[i][idIdx][n][SC_3[m]][0]));
										
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
													//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
													doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(azimuth_ASA[i][idIdx][n][SC_3[m]]*pi/180 - AoMD) * timeVector[i][t] ) );
													//std::cout << "Doppler: " << doppler << std::endl;
													raySumInterferer[i][idIdx-1][clusterIdx][t][u][s] += pol * doppler * exp_arrival * exp_departure;
													if(m + 1 == size_SC_3){
														//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
														raySumInterferer[i][idIdx-1][clusterIdx][t][u][s] = sq_P_over_M * raySumInterferer[i][idIdx-1][clusterIdx][t][u][s];
														//std::cout << "raySumInterferer[i][clusterIdx][t][u][s]: " << raySumInterferer[i][clusterIdx][t][u][s] << std::endl;
													}
												} // End time axis							
											} // End cycle Rays
											clusterIdx++;
										} //end if (for all sub-clusters)
									}else{ // for all clusters from N = 3 onwards
										sq_P_over_M = sqrt(P_n[l] / numOfRays_NLOS);
										//std::cout << "P_n[l]: " << P_n[l] << std::endl;
										//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
										// Cycle through all NLOS Ray
										for(int m = 0; m < numOfRays_NLOS; m++){		
											AoA[0] = sin(elevation_ASA[i][idIdx][n][m]*pi/180) * cos(azimuth_ASA[i][idIdx][n][m]*pi/180);
											AoA[1] = sin(elevation_ASA[i][idIdx][n][m]*pi/180) * sin(azimuth_ASA[i][idIdx][n][m]*pi/180);
											AoA[2] = cos(elevation_ASA[i][idIdx][n][m]*pi/180);
										
											AoD[0] = sin(elevation_ASD[i][idIdx][n][m]*pi/180) * cos(azimuth_ASD[i][idIdx][n][m]*pi/180);
											AoD[1] = sin(elevation_ASD[i][idIdx][n][m]*pi/180) * sin(azimuth_ASD[i][idIdx][n][m]*pi/180);
											AoD[2] = cos(elevation_ASD[i][idIdx][n][m]*pi/180);
										
											complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2])) );
											complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
											// TODO: Include Polarization
											// Replacement for polarization. (Below 4.14)
											//pol = exp( complex<double>(0,uniform(0,2*pi)));
											MSgain = getMSGain(azimuth_ASA[i][idIdx][n][m]*pi/180, elevation_ASA[i][idIdx][n][m]*pi/180);
											BSgain = getBSGain(azimuth_ASD[i][idIdx][n][m]*pi/180, elevation_ASD[i][idIdx][n][m]*pi/180);
											pol = MSgain * BSgain * exp(complex<double>(0, randomPhase[i][idIdx][n][m][0]));
										
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
												//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
												doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(azimuth_ASA[i][idIdx][n][m]*pi/180 - AoMD) * timeVector[i][t] ) );
												//std::cout << "Doppler: " << doppler << std::endl;
												raySumInterferer[i][idIdx-1][clusterIdx][t][u][s] += pol * doppler * exp_arrival * exp_departure;
												if(m + 1 == numOfRays_NLOS){
													//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
													raySumInterferer[i][idIdx-1][clusterIdx][t][u][s] = sq_P_over_M * raySumInterferer[i][idIdx-1][clusterIdx][t][u][s];
													//std::cout << "raySumInterferer[i][idIdx-1][clusterIdx][t][u][s]: " << raySumInterferer[i][idIdx-1][clusterIdx][t][u][s] << std::endl;
												}
											} // End time axis							
										} // End cycle Rays
										clusterIdx++;
									}//end if
								} // End cycle Subclusters
							} // End cycle Clusters
						} // End BS antenna
					} // End MS antenna
				}else{
					// LOS case,
					// Cycle through all Receiver antennas (MS)
					for(int u = 0; u < NumMsAntenna; u++){
						// Cycle through all Transmitter antennas (BS)
						for(int s = 0; s < NumBsAntenna; s++){
							// Cycle through all Paths/Clusters
							clusterIdx_LOS = 0;
							for(int n = 0; n < N_cluster_LOS; n++){
								int L;
								double P_n[3]; // Oversized, but 2 doubles really doesnt matter
								// Two strongest clusters are divided into 3 subclusters!
								if(n < 2){
									L = 3;
									P_n[0] = 10.0/20.0 * clusterPowers[i][idIdx][n];
									P_n[1] =  6.0/20.0 * clusterPowers[i][idIdx][n];
									P_n[2] =  4.0/20.0 * clusterPowers[i][idIdx][n];
								}else{
									L = 1;
									P_n[0] = clusterPowers[i][idIdx][n];
								}
								// Cycle through all subclusters (only a loop for cluster 1 and 2)
								for(int l = 0; l < L; l++){
									if (L == 3){
										if (l == 0){ // for sub-cluster 1 of first two clusters
											sq_P_over_M = sqrt(P_n[l] / size_SC_1);
											//std::cout << "P_n[l]: " << P_n[l] << std::endl;
											//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
											// Cycle through all LOS Rays
											for(int m = 0; m < size_SC_1; m++){		
												AoA[0] = sin(elevation_ASA[i][idIdx][n][SC_1[m]]*pi/180) * cos(azimuth_ASA[i][idIdx][n][SC_1[m]]*pi/180);
												AoA[1] = sin(elevation_ASA[i][idIdx][n][SC_1[m]]*pi/180) * sin(azimuth_ASA[i][idIdx][n][SC_1[m]]*pi/180);
												AoA[2] = cos(elevation_ASA[i][idIdx][n][SC_1[m]]*pi/180);
										
												AoD[0] = sin(elevation_ASD[i][idIdx][n][SC_1[m]]*pi/180) * cos(azimuth_ASD[i][idIdx][n][SC_1[m]]*pi/180);
												AoD[1] = sin(elevation_ASD[i][idIdx][n][SC_1[m]]*pi/180) * sin(azimuth_ASD[i][idIdx][n][SC_1[m]]*pi/180);
												AoD[2] = cos(elevation_ASD[i][idIdx][n][SC_1[m]]*pi/180);
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
												// TODO: Include Polarization
												// Replacement for polarization. (Below 4.14)
												//pol = exp( complex<double>(0,uniform(0,2*pi)));
												MSgain = getMSGain(azimuth_ASA[i][idIdx][n][SC_1[m]]*pi/180, elevation_ASA[i][idIdx][n][SC_1[m]]*pi/180);
												BSgain = getBSGain(azimuth_ASD[i][idIdx][n][SC_1[m]]*pi/180, elevation_ASD[i][idIdx][n][SC_1[m]]*pi/180);
												pol = MSgain * BSgain * exp(complex<double>(0, randomPhase[i][idIdx][n][SC_1[m]][0]));
										
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
													//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
													doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(azimuth_ASA[i][idIdx][n][SC_1[m]]*pi/180 - AoMD) * timeVector[i][t] ) );
													//std::cout << "Doppler: " << doppler << std::endl;
													raySumInterferer_LOS[i][idIdx-1][clusterIdx_LOS][t][u][s] += pol * doppler * exp_arrival * exp_departure;
													if(m + 1 == size_SC_1){
														//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
														double K = sigma_kf_LOS[i][idIdx]; 
														raySumInterferer_LOS[i][idIdx-1][clusterIdx_LOS][t][u][s] = (sqrt(1/(K + 1))) * sq_P_over_M * raySumInterferer_LOS[i][idIdx-1][clusterIdx_LOS][t][u][s];
														//std::cout << "raySumInterferer_LOS[i][idIdx-1][clusterIdx][t][u][s]: " << raySumInterferer_LOS[i][idIdx-1][clusterIdx_LOS][t][u][s] << std::endl;
													}
												} // End time axis							
											} // End cycle Rays
											clusterIdx_LOS++;
										} else if (l == 1){ // for sub-cluster 2 of first two clusters
											sq_P_over_M = sqrt(P_n[l] / size_SC_2);
											//std::cout << "P_n[l]: " << P_n[l] << std::endl;
											//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
											// Cycle through all LOS Rays
											for(int m = 0; m < size_SC_2; m++){		
												AoA[0] = sin(elevation_ASA[i][idIdx][n][SC_2[m]]*pi/180) * cos(azimuth_ASA[i][idIdx][n][SC_2[m]]*pi/180);
												AoA[1] = sin(elevation_ASA[i][idIdx][n][SC_2[m]]*pi/180) * sin(azimuth_ASA[i][idIdx][n][SC_2[m]]*pi/180);
												AoA[2] = cos(elevation_ASA[i][idIdx][n][SC_2[m]]*pi/180);
										
												AoD[0] = sin(elevation_ASD[i][idIdx][n][SC_2[m]]*pi/180) * cos(azimuth_ASD[i][idIdx][n][SC_2[m]]*pi/180);
												AoD[1] = sin(elevation_ASD[i][idIdx][n][SC_2[m]]*pi/180) * sin(azimuth_ASD[i][idIdx][n][SC_2[m]]*pi/180);
												AoD[2] = cos(elevation_ASD[i][idIdx][n][SC_2[m]]*pi/180);
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
												// TODO: Include Polarization
												// Replacement for polarization. (Below 4.14)
												//pol = exp( complex<double>(0,uniform(0,2*pi)));
												MSgain = getMSGain(azimuth_ASA[i][idIdx][n][SC_2[m]]*pi/180, elevation_ASA[i][idIdx][n][SC_2[m]]*pi/180);
												BSgain = getBSGain(azimuth_ASD[i][idIdx][n][SC_2[m]]*pi/180, elevation_ASD[i][idIdx][n][SC_2[m]]*pi/180);
												pol = MSgain * BSgain * exp(complex<double>(0, randomPhase[i][idIdx][n][SC_2[m]][0]));
										
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
													//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
													doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(azimuth_ASA[i][idIdx][n][SC_2[m]]*pi/180 - AoMD) * timeVector[i][t] ) );
													//std::cout << "Doppler: " << doppler << std::endl;
													raySumInterferer_LOS[i][idIdx-1][clusterIdx_LOS][t][u][s] += pol * doppler * exp_arrival * exp_departure;
													if(m + 1 == size_SC_2){
														//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
														double K = sigma_kf_LOS[i][idIdx]; 
														raySumInterferer_LOS[i][idIdx-1][clusterIdx_LOS][t][u][s] = (sqrt(1/(K + 1))) * sq_P_over_M * raySumInterferer_LOS[i][idIdx-1][clusterIdx_LOS][t][u][s];
														//std::cout << "raySumInterferer_LOS[i][idIdx-1][clusterIdx_LOS][t][u][s]: " << raySumInterferer[i][idIdx-1][clusterIdx_LOS][t][u][s] << std::endl;
													}
												} // End time axis							
											} // End cycle Rays
											clusterIdx_LOS++;
										} else { //for sub-cluster 3 of first two clusters
											sq_P_over_M = sqrt(P_n[l] / size_SC_3);
											//std::cout << "P_n[l]: " << P_n[l] << std::endl;
											//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
											// Cycle through all LOS Rays
											for(int m = 0; m < size_SC_3; m++){		
												AoA[0] = sin(elevation_ASA[i][idIdx][n][SC_3[m]]*pi/180) * cos(azimuth_ASA[i][idIdx][n][SC_3[m]]*pi/180);
												AoA[1] = sin(elevation_ASA[i][idIdx][n][SC_3[m]]*pi/180) * sin(azimuth_ASA[i][idIdx][n][SC_3[m]]*pi/180);
												AoA[2] = cos(elevation_ASA[i][idIdx][n][SC_3[m]]*pi/180);
										
												AoD[0] = sin(elevation_ASD[i][idIdx][n][SC_3[m]]*pi/180) * cos(azimuth_ASD[i][idIdx][n][SC_3[m]]*pi/180);
												AoD[1] = sin(elevation_ASD[i][idIdx][n][SC_3[m]]*pi/180) * sin(azimuth_ASD[i][idIdx][n][SC_3[m]]*pi/180);
												AoD[2] = cos(elevation_ASD[i][idIdx][n][SC_3[m]]*pi/180);
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
												// TODO: Include Polarization
												// Replacement for polarization. (Below 4.14)
												//pol = exp( complex<double>(0,uniform(0,2*pi)));
												MSgain = getMSGain(azimuth_ASA[i][idIdx][n][SC_3[m]]*pi/180, elevation_ASA[i][idIdx][n][SC_3[m]]*pi/180);
												BSgain = getBSGain(azimuth_ASD[i][idIdx][n][SC_3[m]]*pi/180, elevation_ASD[i][idIdx][n][SC_3[m]]*pi/180);
												pol = MSgain * BSgain * exp(complex<double>(0, randomPhase[i][idIdx][n][SC_3[m]][0]));
										
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
													//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
													doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(azimuth_ASA[i][idIdx][n][SC_3[m]]*pi/180 - AoMD) * timeVector[i][t] ) );
													//std::cout << "Doppler: " << doppler << std::endl;
													raySumInterferer_LOS[i][idIdx-1][clusterIdx_LOS][t][u][s] += pol * doppler * exp_arrival * exp_departure;
													if(m + 1 == size_SC_3){
														//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
														double K = sigma_kf_LOS[i][idIdx]; 
														raySumInterferer_LOS[i][idIdx-1][clusterIdx_LOS][t][u][s] = (sqrt(1/(K + 1))) * sq_P_over_M * raySumInterferer_LOS[i][idIdx-1][clusterIdx_LOS][t][u][s];
														//std::cout << "raySumInterferer_LOS[i][idIdx-1][clusterIdx_LOS][t][u][s]: " << raySumInterferer[i][idIdx-1][clusterIdx_LOS][t][u][s] << std::endl;
													}
												} // End time axis							
											} // End cycle Rays
											clusterIdx_LOS++;
										} //end if (for all sub-clusters)
									}else{ // for all clusters from N = 3 onwards
										sq_P_over_M = sqrt(P_n[l] / numOfRays_LOS);
										//std::cout << "P_n[i]: " << P_n[i] << std::endl;
										//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
										// Cycle through all LOS Rays
										for(int m = 0; m < numOfRays_LOS; m++){		
											AoA[0] = sin(elevation_ASA[i][idIdx][n][m]*pi/180) * cos(azimuth_ASA[i][idIdx][n][m]*pi/180);
											AoA[1] = sin(elevation_ASA[i][idIdx][n][m]*pi/180) * sin(azimuth_ASA[i][idIdx][n][m]*pi/180);
											AoA[2] = cos(elevation_ASA[i][idIdx][n][m]*pi/180);
										
											AoD[0] = sin(elevation_ASD[i][idIdx][n][m]*pi/180) * cos(azimuth_ASD[i][idIdx][n][m]*pi/180);
											AoD[1] = sin(elevation_ASD[i][idIdx][n][m]*pi/180) * sin(azimuth_ASD[i][idIdx][n][m]*pi/180);
											AoD[2] = cos(elevation_ASD[i][idIdx][n][m]*pi/180);
										
											complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2])) );
											complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
											// TODO: Include Polarization
											// Replacement for polarization. (Below 4.14)
											//pol = exp( complex<double>(0,uniform(0,2*pi)));
											MSgain = getMSGain(azimuth_ASA[i][idIdx][n][m]*pi/180, elevation_ASA[i][idIdx][n][m]*pi/180);
											BSgain = getBSGain(azimuth_ASD[i][idIdx][n][m]*pi/180, elevation_ASD[i][idIdx][n][m]*pi/180);
											pol = MSgain * BSgain * exp(complex<double>(0, randomPhase[i][idIdx][n][m][0]));
										
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
												//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
												doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(azimuth_ASA[i][idIdx][n][m]*pi/180 - AoMD) * timeVector[i][t] ) );
												//std::cout << "Doppler: " << doppler << std::endl;
												raySumInterferer_LOS[i][idIdx-1][clusterIdx_LOS][t][u][s] += pol * doppler * exp_arrival * exp_departure;
												if(m + 1 == numOfRays_LOS){
													//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
													double K = sigma_kf_LOS[i][idIdx]; 
													raySumInterferer_LOS[i][idIdx-1][clusterIdx_LOS][t][u][s] = (sqrt(1/(K + 1))) * sq_P_over_M * raySumInterferer_LOS[i][idIdx-1][clusterIdx_LOS][t][u][s];
													//std::cout << "raySumInterferer_LOS[i][idIdx-1][clusterIdx_LOS][t][u][s]: " << raySumInterferer_LOS[i][idIdx-1][clusterIdx_LOS][t][u][s] << std::endl;
												}
											} // End time axis							
										} // End cycle Rays
										clusterIdx_LOS++;
									}//end if
								} // End cycle Subclusters
								if(n == 0){ // for adding the additional LOS component, according to formula 7-61 in METIS 1.2
									std::cout << "n: " << n << " LOS Computation for MS " << i << " and BS " << idIdx << std::endl;
									double K = sigma_kf_LOS[i][idIdx]; 						// K-factor in linear scale 
									AoA[0] = sin(ZoA_LOS_dir[i][idIdx]) * cos(AoA_LOS_dir[i][idIdx]);
									AoA[1] = sin(ZoA_LOS_dir[i][idIdx]) * sin(AoA_LOS_dir[i][idIdx]);
									AoA[2] = cos(ZoA_LOS_dir[i][idIdx]);
										
									AoD[0] = sin(ZoD_LOS_dir[i][idIdx]) * cos(AoD_LOS_dir[i][idIdx]);
									AoD[1] = sin(ZoD_LOS_dir[i][idIdx]) * sin(AoD_LOS_dir[i][idIdx]);
									AoD[2] = cos(ZoD_LOS_dir[i][idIdx]);
										
									complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * MsAntennaPosition[i][u][0] + AoA[1] * MsAntennaPosition[i][u][1] + AoA[2] * MsAntennaPosition[i][u][2])) );
									complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * OwnBsAntennaPosition[0][s][0] + AoD[1] * OwnBsAntennaPosition[0][s][1] + AoD[2] * OwnBsAntennaPosition[0][s][2])) );
										
									// TODO: Include Polarization in generic form
									// Replacement for polarization. (Below 4.14)
									//pol = exp( complex<double>(0,uniform(0,2*pi)));
									MSgain = getMSGain(AoA_LOS_dir[i][idIdx], ZoA_LOS_dir[i][idIdx]);
									BSgain = getBSGain(AoD_LOS_dir[i][idIdx], ZoD_LOS_dir[i][idIdx]);
									pol = MSgain * BSgain * exp(complex<double>(0, randomPhase_LOS[i][idIdx]));
										
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
										//std::cout << "n: " << n << " t: " << t << " Idx: " << clusterIdx << std::endl;
										doppler = exp( complex<double>(0,k_0 * MSVelMag[i] * cos(AoA_LOS_dir[i][idIdx] - AoMD) * timeVector[i][t] ) );
										//std::cout << "Doppler: " << doppler << std::endl;
										raySumInterferer_LOS[i][idIdx-1][0][t][u][s] += (sqrt(K / (K + 1))) * pol * doppler * exp_arrival * exp_departure;
									} // End time axis
								}
							} // End cycle Clusters
						} // End BS antenna
					} // End MS antenna
				} //End if condition for interferer BS
				idIdx++;
			} //End if condition for own or interferer BS
		} //End loop for all BSs
	} // End Links/MS
	
	std::cout << "FINISHED MAIN LOOP for BS: " << bsId << std::endl;

	//output2 << "Init 3 METIS at BS " << bsId << " with rand: " << normal(0,1) << std::endl;
		
	//------------------------------------------------------------------
	// Apply Fourier transform, to get time/frequency domain from time/delay
	
	std::cout << "START FOURIER TRANSFORM for BS: " << bsId << std::endl;
	
	double pathloss, dist3D;
	ofstream TimeFrequency;
	//ofstream TimeFrequencyInteferer;
	
	double delay_SC_1 = 5 * pow(10,-9); // delay for sub-cluster 1 (7-60)
	double delay_SC_2 = 10 * pow(10,-9); // delay for sub-cluster 2 (7-60)
	
	for(int i = 0; i < numberOfMobileStations; i++){
		//TimeFrequency.open ("/home/sahari/GSCM_Channel_model/abstractLTEChannelModel/results/TimeFrequency_BS_" + std::to_string( (long long) bsId ) + "_" + std::to_string((MSPos[i].x * 100)) + ".txt");
		dist2D = sqrt(pow((xPos - MSPos[i].x),2) + pow((yPos - MSPos[i].y),2));
		dist3D = sqrt(pow(dist2D,2) + pow((heightBS - heightUE),2));
		complex<double> res = complex<double>(0.0,0.0);
		for(int t = 0; t < timeSamples; t++){
			for(int f = 0; f < downRBs; f++){
				double freq_ = freq_c + f*180000;
				res = complex<double>(0.0,0.0);
				if(LOSCondition[i][0]){
					pathloss = CalcPathloss(dist2D, dist3D, true);
					for(int u = 0; u < NumMsAntenna; u++){
						for(int s = 0; s < NumBsAntenna; s++){
							for(int n = 0; n < N_cluster_LOS; n++){
								if(n < 2){ // add additional sub-cluster delays at this stage (7-60 in METIS 1.2)
									res = res + raySum_LOS[i][n][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * clusterDelays_LOS[i][0][n]) );
									res = res + raySum_LOS[i][n+1][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (clusterDelays_LOS[i][0][n] + delay_SC_1)) );
									res = res + raySum_LOS[i][n+2][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (clusterDelays_LOS[i][0][n] + delay_SC_2)) );
									res = res + raySum_LOS[i][n+3][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * clusterDelays_LOS[i][0][n+1]) );
									res = res + raySum_LOS[i][n+4][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (clusterDelays_LOS[i][0][n+1] + delay_SC_1)) );
									res = res + raySum_LOS[i][n+5][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (clusterDelays_LOS[i][0][n+1] + delay_SC_2)) );
									n++;
								}else{
									res = res + raySum_LOS[i][n + 4][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * clusterDelays_LOS[i][0][n]) );
								}
							}
						}
					}
					double tempRes = pow(res.real(),2) + pow(res.imag(),2);
					SINRtable[i][t][f] = pathloss * tempRes;
					//TimeFrequency << pathloss * tempRes << " ";
				}else{
					pathloss = CalcPathloss(dist2D, dist3D, false);
					for(int u = 0; u < NumMsAntenna; u++){
						for(int s = 0; s < NumBsAntenna; s++){
							for(int n = 0; n < N_cluster_NLOS; n++){
								if(n < 2){ // add additional sub-cluster delays at this stage (7-60 in METIS 1.2)
									res = res + raySum[i][n][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * clusterDelays[i][0][n]) );
									res = res + raySum[i][n+1][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (clusterDelays[i][0][n] + delay_SC_1)) );
									res = res + raySum[i][n+2][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (clusterDelays[i][0][n] + delay_SC_2)) );
									res = res + raySum[i][n+3][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * clusterDelays[i][0][n+1]) );
									res = res + raySum[i][n+4][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (clusterDelays[i][0][n+1] + delay_SC_1)) );
									res = res + raySum[i][n+5][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (clusterDelays[i][0][n+1] + delay_SC_2)) );
									n++;
								}else{
									res = res + raySum[i][n + 4][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * clusterDelays[i][0][n]) );
								}
							} //End Tx antenna
						} //End Rx antenna
					} //End clusters
					double tempRes = pow(res.real(),2) + pow(res.imag(),2);
					SINRtable[i][t][f] = pathloss * tempRes;
					//TimeFrequency << pathloss * tempRes << " ";
				} //End if condition for LOS/NLOS
			} //End loop for RBs
			TimeFrequency << "\n\n";
		}
		//TimeFrequency << pathloss << "\n";
		//TimeFrequency.close();
	}
	
	if(numOfInterferers>0){
		//Only compute SINRneighbour values if there actually are any 
		// neighbours to influence transmissions
		for(int i = 0; i < numberOfMobileStations; i++){ //TODO: Correct the implementation for each Tx-Rx antenna
			int idIdx = 0;
			//Fourier Transform for interferers:
			for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
				if(it->first == bsId){
					// Skip own BS
					continue;
				}else{
					//TimeFrequencyInteferer.open ("/home/sahari/GSCM_Channel_model/abstractLTEChannelModel/results/TimeFrequency_Interferer_BS_" + std::to_string( (long long) bsId ) + "_" + std::to_string((MSPos[i].x * 100)) + ".txt");
					dist2D = sqrt(pow((it->second.x - MSPos[i].x),2) + pow((it->second.y - MSPos[i].y),2));
					dist3D = sqrt(pow(dist2D,2) + pow((heightBS - heightUE),2));
					complex<double> res = complex<double>(0.0,0.0);
					for(int t = 0; t < timeSamples; t++){
						for(int f = 0; f < downRBs; f++){
							double freq_ = freq_c + f*180000;
							res = complex<double>(0.0,0.0);
							if(LOSCondition[i][idIdx+1]){
								pathloss = CalcPathloss(dist2D, dist3D, true);
								for(int u = 0; u < NumMsAntenna; u++){
									for(int s = 0; s < NumBsAntenna; s++){
										for(int n = 0; n < N_cluster_LOS; n++){
											if(n < 2){ // add additional sub-cluster delays at this stage (7-60 in METIS 1.2)
												res = res + raySumInterferer_LOS[i][idIdx][n][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * clusterDelays_LOS[i][idIdx+1][n]) );
												res = res + raySumInterferer_LOS[i][idIdx][n+1][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (clusterDelays_LOS[i][idIdx+1][n] + delay_SC_1)) );
												res = res + raySumInterferer_LOS[i][idIdx][n+2][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (clusterDelays_LOS[i][idIdx+1][n] + delay_SC_2)) );
												res = res + raySumInterferer_LOS[i][idIdx][n+3][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * clusterDelays_LOS[i][idIdx+1][n+1]) );
												res = res + raySumInterferer_LOS[i][idIdx][n+4][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (clusterDelays_LOS[i][idIdx+1][n+1] + delay_SC_1)) );
												res = res + raySumInterferer_LOS[i][idIdx][n+5][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (clusterDelays_LOS[i][idIdx+1][n+1] + delay_SC_2)) );
												n++;
											}else{
												res = res + raySumInterferer_LOS[i][idIdx][n + 4][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * clusterDelays_LOS[i][idIdx+1][n]) );
											}
										}
									}
								}
								double tempRes = pow(res.real(),2) + pow(res.imag(),2);
								SINRneighbour[i][idIdx][t][f] = pathloss * tempRes;
								//TimeFrequencyInteferer << pathloss * tempRes << " ";
							}else{
								pathloss = CalcPathloss(dist2D, dist3D, false);
								for(int u = 0; u < NumMsAntenna; u++){
									for(int s = 0; s < NumBsAntenna; s++){
										for(int n = 0; n < N_cluster_NLOS; n++){
											if(n < 2){ // add additional sub-cluster delays at this stage (7-60 in METIS 1.2)
												res = res + raySumInterferer[i][idIdx][n][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * clusterDelays[i][idIdx+1][n]) );
												res = res + raySumInterferer[i][idIdx][n+1][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (clusterDelays[i][idIdx+1][n] + delay_SC_1)) );
												res = res + raySumInterferer[i][idIdx][n+2][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (clusterDelays[i][idIdx+1][n] + delay_SC_2)) );
												res = res + raySumInterferer[i][idIdx][n+3][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * clusterDelays[i][idIdx+1][n+1]) );
												res = res + raySumInterferer[i][idIdx][n+4][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (clusterDelays[i][idIdx+1][n+1] + delay_SC_1)) );
												res = res + raySumInterferer[i][idIdx][n+5][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * (clusterDelays[i][idIdx+1][n+1] + delay_SC_2)) );
												n++;
											}else{
												res = res + raySumInterferer[i][idIdx][n + 4][t][u][s] * exp( complex<double>(0.0, -2.0 * pi * freq_ * clusterDelays[i][idIdx+1][n]) );
											}
										}
									}
								}
								double tempRes = pow(res.real(),2) + pow(res.imag(),2);
								SINRneighbour[i][idIdx][t][f] = pathloss * tempRes;
								//TimeFrequencyInteferer << pathloss * tempRes << " ";
							}
						}
						//TimeFrequencyInteferer << "\n";
					}
					idIdx++;
				}
			}
			//TimeFrequencyInteferer.close();
		}
	
	}
	
	std::cout << "FINISHED FOURIER TRANSFORM for BS: " << bsId << std::endl;
	
	//------------------------------------------------------------------
	
	//---------------------------------------------------------------------
	// Delete all heap allocated local variables

	for(int i = 0; i < numberOfMobileStations; i++){
		delete[] MSVelDir[i];
		for(int s=0; s<NumMsAntenna; s++){
			delete[] MsAntennaPosition[i][s];
		}
		delete[] MsAntennaPosition[i];
	}
	delete[] MSVelDir;
	delete[] MSVelMag;
	for(int m = 0; m < numberOfMobileStations; m++){
		for(int i = 0; i < (N_cluster_LOS + 4); i++){
			for(int j = 0; j < timeSamples; j++){
				for(int k = 0; k < NumMsAntenna; k++){
					delete[] raySum_LOS[m][i][j][k];
				}
				delete[] raySum_LOS[m][i][j];
			}
			delete[] raySum_LOS[m][i];
		}
		delete[] raySum_LOS[m];
		for(int i = 0; i < (N_cluster_NLOS + 4); i++){
			for(int j = 0; j < timeSamples; j++){
				for(int k = 0; k < NumMsAntenna; k++){
					delete[] raySum[m][i][j][k];
				}
				delete[] raySum[m][i][j];
			}
			delete[] raySum[m][i];
		}
		delete[] raySum[m];
	}
	delete[] raySum;
	delete[] raySum_LOS;
	if(numOfInterferers>0){
		for (int m = 0; m < numberOfMobileStations; m++){
			for(int i = 0; i < numOfInterferers; i++){
				for(int j = 0; j < (N_cluster_LOS + 4); j++){
					for(int k = 0; k < timeSamples; k++){
						for(int l = 0; l < NumMsAntenna; l++){
							delete[] raySumInterferer_LOS[m][i][j][k][l];
							delete[] raySumInterferer[m][i][j][k][l];
						}
						delete[] raySumInterferer_LOS[m][i][j][k];
						delete[] raySumInterferer[m][i][j][k];
					}
					delete[] raySumInterferer_LOS[m][i][j];
					delete[] raySumInterferer[m][i][j];
				}
				delete[] raySumInterferer_LOS[m][i];
				delete[] raySumInterferer[m][i];
			}
			delete[] raySumInterferer_LOS[m];
			delete[] raySumInterferer[m];
		}
		delete[] raySumInterferer_LOS;
		delete[] raySumInterferer;
	}
	delete[] MsAntennaPosition;
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
	int r = initModule->par("cellRadiusMETIS");
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
		grid[i] = new double*[2*r + 20];
		tmpX[i] = new double*[2*r + 20];
		tmpY[i] = new double*[2*r + 20];
		for(int j = 0; j < 2*r + 20; j++){
			grid[i][j] = new double[2*r + 20];
			tmpX[i][j] = new double[2*r + 20];
			tmpY[i][j] = new double[2*r + 20];
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
		for(int j = 0; j < 2*r + 20; j++){
			for(int k = 0; k < 2*r + 20; k++){
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
		for(int j = 10; j < 2*r + 10; j++){
			// Forall y points except offset
			for(int k = 10; k < 2*r + 10; k++){
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
		for(int j = 10; j < 2*r + 10; j++){
			// For all y points except offset
			for(int k = 10; k < 2*r + 10; k++){
				// Filter 100 points
				for(int l = 0; l < 10; l++){
					tmpY[i][k][j] += tmpX[i][k - l][j] * filter[i][l];
				}
			}
		}
	}
		
	// For each Link generate Large scale parameters
	double middle = (2*r + 10)/2;
	for(size_t l = 0; l < receivers.size(); l++){
		for(size_t i=0; i<senders.size(); i++){
			for(int j=0; j<7; j++){
				correlation[l][i][j] = tmpY[j][(int) (receivers[l].x + middle - senders[i].x)][(int) (receivers[l].y + middle - senders[i].y)];
			}
		}
	}
	// Clean up all local variables allocated on the heap
	// TODO Look into the possibility of actually putting those variables 
	// directly on the stack, instead of using entirely unnecessary 
	// dynamic memory allocation
	for(int i = 0; i < 7; i++){
		for(int j = 0; j < 2*r + 20; j++){
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
	int r = initModule->par("cellRadiusMETIS");
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
		grid[i] = new double*[2*r + 20];
		tmpX[i] = new double*[2*r + 20];
		tmpY[i] = new double*[2*r + 20];
		for(int j = 0; j < 2*r + 20; j++){
			grid[i][j] = new double[2*r + 20];
			tmpX[i][j] = new double[2*r + 20];
			tmpY[i][j] = new double[2*r + 20];
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
		for(int j = 0; j < 2*r + 20; j++){
			for(int k = 0; k < 2*r + 20; k++){
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
		for(int j = 10; j < 2*r + 10; j++){
			// Forall y points except offset
			for(int k = 10; k < 2*r + 10; k++){
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
		for(int j = 10; j < 2*r + 10; j++){
			// Forall y points except offset
			for(int k = 10; k < 2*r + 10; k++){
				// Filter 100 points
				for(int l = 0; l < 10; l++){
					tmpY[i][k][j] += tmpX[i][k - l][j] * filter[i][l];
				}
			}
		}
	}
		
	// For each Link generate Large scale parameters
	double middle = (2*r + 10)/2;
	for(size_t l = 0; l < receivers.size(); l++){
		for(size_t i=0; i<senders.size(); i++){
			for(int j=0; j<6; j++){
				std::cout << l << " " << i << " " << j << std::endl; 
				correlation[l][i][j] = tmpY[j][(int) (receivers[l].x + middle - senders[i].x)][(int) (receivers[l].y + middle - senders[i].y)];
			}
		}
	}
	// Clean up all local variables allocated on the heap
	// TODO Look into the possibility of actually putting those variables 
	// directly on the stack, instead of using entirely unnecessary 
	// dynamic memory allocation
	for(int i = 0; i < 6; i++){
		for(int j = 0; j < 2*r + 20; j++){
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
	
	if(LOS){
		if(dist2D < 10.0 || dist2D > 5000.0){
			// LOS: 2D distance must be between 10 and 5000 meters.
			std::cout << "WARNING: LOS Distance outside of allowed boundaries!";
		}else if(heightUE >= 1.5 && ( dist2D < distBP )){
			// If the 2D distance is below the breakpoint distance (see Document 36873, page 23 for reference)
			pathloss = 22.0 * log10(dist3D) + 28.0 + 20.0 * log10(freq_c/1000000000);
		}else{
			pathloss = 22.0 * log10(dist3D) + 28.0 + 20.0 * log10(freq_c/1000000000) - 9.0 * log10( pow(distBP,2) + pow((heightBS - heightUE),2) );
		}
	}else{
		//pl_a = 36.7 * log10(dist3D) + 23.15 + 26 * log10(freq_c / 1000000000) - 0.3 * (heightUE);
		//pl_a = pow(10,pl_a/10);
		
		if(dist2D < 10.0 || dist2D > 2000.0 || heightUE > 22.5){
			// NLOS: 2D distance must be between 10 and 2000 meters.
			std::cout << "WARNING: NLOS Distance outside of allowed boundaries!";
		}else if(heightUE >= 1.5 && ( dist2D < distBP )){
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
	delete msg;
}

double METISChannel::calcSINR(int RB, vector<double> &power, vector<Position> &pos, vector<int> &bsId_, bool up, int msId){
	NeighbourIdMatching *neighbourIdMatching;
	cModule *cell = initModule->getParentModule()->getParentModule();
	neighbourIdMatching = new NeighbourIdMatching(bsId, maxNumberOfNeighbours, cell);
	
	int SINRCounter = 3; //originally set to std::round( simTime().dbl() * 1000.0* 4.0) - 1 

	double interference = 0;
	for(int i = 0; i < neighbourIdMatching->numberOfNeighbours() - 1; i++){
		interference += SINRneighbour[msId][i][SINRCounter][RB];
	}
	delete neighbourIdMatching;
	interference += getTermalNoise(300,180000);
	// Convert to db scale
	//std::cout << 10 * log10( SINRtable[msId][SINRcounter][RB] / interference ) << std::endl;
	//std::cout << "Gain: " << 10 * log10( SINRtable[msId][SINRcounter][RB] ) << "  MS: " << msId << std::endl;
	//std::cout << "Interference: " << interference << "  MS: " << msId << std::endl;
	//std::cout << "Simtime: " << simTime() << " Counter: " << std::round( simTime().dbl() *1000.0*4.0) << std::endl;
	return 10 * log10( SINRtable[msId][SINRCounter][RB] / interference );
	//return 10 * log10( SINRtable[msId][SINRcounter][RB] / getTermalNoise(300,180000) );
}

vec METISChannel::calcSINR(vector<double> &power, vector<Position> &pos, vector<int> &bsId_, bool up, int msId){
	vec result(downRBs);
	for(int i = 0; i < downRBs; i++){
		result.set(i,calcSINR(i, power, pos, bsId_, up, msId));
	}
	//Placeholder
	return result;
}

void METISChannel::updateChannel(Position** msPos){
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
			for(int j=0; j<numOfInterferers; j++){
				for(int k=0; k<timeSamples; k++){
					delete[] SINRneighbour[i][j][k];
				}
				delete[] SINRneighbour[i][j];
			}
			delete[] SINRneighbour[i];
			for(int j=0; j<timeSamples; j++){
				delete[] SINRtable[i][j];
			}
			delete[] SINRtable[i];
			delete[] timeVector[i];
		}
		delete neighbourIdMatching; 
		delete[] SINRneighbour;
		delete[] SINRtable;
		for(int i=0; i<NumBsAntenna; i++){
			delete[] OwnBsAntennaPosition[0][i];
		}
		delete[] timeVector;
		delete[] OwnBsAntennaPosition[0];
		delete[] OwnBsAntennaPosition;
	
	}

}
