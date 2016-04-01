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

/**
* Function, that initializes all large scale and small scale parameters according to METIS specifications.
* @param module OMNeT++ module, which calls this function to allow later .ini access.
* @param msPositions Positions of the mobile stations.
* @param neighbourPositions Positions of the Neighbour BS.
* @return True if initialization was successful, false otherwise
*/
bool METISChannel::init(cSimpleModule* module, Position** msPositions, std::map <int,Position> neighbourPositions){
	// Basic Initialization
	neighbourPos = neighbourPositions;
	maxNumberOfNeighbours = module->par("maxNumberOfNeighbours");
	bsId = module->par("bsId");
	tti = module->par("tti");
    	numberOfMobileStations = module->par("numberOfMobileStations");
    	xPos = module->par("xPos");
	yPos = module->par("yPos");
	upRBs = module->par("upResourceBlocks");
	downRBs = module->par("downResourceBlocks");
    	bs_antenna_bearing = new double[3];
    	bs_antenna_downtilt = new double[3];
    	bs_antenna_slant = new double[3];
    	ms_antenna_bearing = new double[numberOfMobileStations];
    	ms_antenna_downtilt = new double[numberOfMobileStations];
    	ms_antenna_slant = new double[numberOfMobileStations];
    	LOSCondition = new bool*[numberOfMobileStations];
    	MSPos = new Position[numberOfMobileStations];
    	speedOfLight = 299792458.0;
   	initModule = module;
    	freq_c = module->par("CarrierFrequency");
    	SINRcounter = 0;
    	NumTxAntenna = module->par("NumTxAntenna");
    	NumRxAntenna = module->par("NumRxAntenna");
    	heightUE = module->par("OutdoorHeightUE");
   	heightBS = module->par("BsHeight");
    
    	MSVelMag = new double[numberOfMobileStations];
    	MSVelDir = new double*[numberOfMobileStations];
    
    	TxAntennaPosition = new double**[1];
    	RxAntennaPosition = new double**[numberOfMobileStations];
    
    	double wavelength = speedOfLight / freq_c;
	double dist2D;
    
    	//string out2 = std::to_string( (long long) bsId ) + "valid.txt";
    	//ofstream output2(out2);
    	//output2 << "Init METIS at BS " << bsId << " with rand: " << normal(0,1) << std::endl;
    
    	// randomly generate direction and magnitude of movement vector
    	// TODO: sync with OMNeT++ part
    	for(int i = 0; i < numberOfMobileStations; i++){
		MSVelDir[i] = new double[3];
		double vel = module->par("Velocity");
		MSVelMag[i] = (vel*1000)/3600; // converting velocity from km/h to m/s 
		MSVelDir[i][0] = 1.0; // +1 for moving towards the BS, -1 for moving away from the BS
		MSVelDir[i][1] = 0.0; // set to zero b/c MS is not moving in y or z-directions
		MSVelDir[i][2] = 0.0; // same as above; originally all are set to uniform(-1.0,1.0)
		double dirMag = sqrt(pow(MSVelDir[i][0],2) + pow(MSVelDir[i][1],2) + pow(MSVelDir[i][2],2));
		MSVelDir[i][0] = MSVelDir[i][0] / dirMag;
		MSVelDir[i][1] = MSVelDir[i][1] / dirMag;
		MSVelDir[i][2] = MSVelDir[i][2] / dirMag;
	}
	std::cout << "Length NeighbourPositions: " << neighbourPositions.size() << std::endl;
	    
    	// Find the neighbours and store the pair (bsId, position in data structures) in a map
	NeighbourIdMatching *neighbourIdMatching;
    	cModule *cell = module->getParentModule()->getParentModule();
    	neighbourIdMatching = new NeighbourIdMatching(bsId, maxNumberOfNeighbours, cell);
    
    	// Actually, this counts the own BS as well, so substract 1 
    	numOfInterferers = neighbourIdMatching->numberOfNeighbours() - 1;
    
    	int fromBsId = neighbourIdMatching->getDataStrId(bsId);
    
    	// Copy MS Positions
    	for(int i = 0; i < numberOfMobileStations; i++){
		MSPos[i].x = msPositions[fromBsId][i].x;
		MSPos[i].y = msPositions[fromBsId][i].y;
	}
	
	// One Channel module per base station
    	// Half wavelength distance between antennas; give the position of Tx and Rx antennas in GCS
	// For even value of NumTxAntenna, the antenna elements will be equally spaced around the center of Tx
    	for(int i = 0; i < 1; i++){
		TxAntennaPosition[i] = new double*[NumTxAntenna];
		for(int j = 0; j < (NumTxAntenna/2); j++){
			TxAntennaPosition[i][j] = new double[3];
			TxAntennaPosition[i][j][0] = xPos - (0.25 * wavelength * (NumTxAntenna - 1 - (j*2.0)));
			TxAntennaPosition[i][j][1] = yPos;
			TxAntennaPosition[i][j][2] = heightBS;
		}

		for(int j = (NumTxAntenna/2); j < NumTxAntenna; j++){
			TxAntennaPosition[i][j] = new double[3];
			TxAntennaPosition[i][j][0] = TxAntennaPosition[i][j-1][0] + (0.5 * wavelength);
			TxAntennaPosition[i][j][1] = yPos;
			TxAntennaPosition[i][j][2] = heightBS;
		}

	}
	
	// Do the same for Rx antennas
	// Mobile stations should be randomly rotated..?
    	for(int i = 0; i < numberOfMobileStations; i++){
		RxAntennaPosition[i] = new double*[NumRxAntenna];
		for(int j = 0; j < (NumRxAntenna/2); j++){
			RxAntennaPosition[i][j] = new double[3];
			RxAntennaPosition[i][j][0] = MSPos[i].x - (0.25 * wavelength * (NumRxAntenna - 1 - (j*2.0)));
			RxAntennaPosition[i][j][1] = MSPos[i].y;
			RxAntennaPosition[i][j][2] = heightUE;
		}
		for(int j = (NumRxAntenna/2); j < NumRxAntenna ; j++){
			RxAntennaPosition[i][j] = new double[3];
			RxAntennaPosition[i][j][0] = RxAntennaPosition[i][j-1][0] + (0.5 * wavelength);
			RxAntennaPosition[i][j][1] = MSPos[i].y;
			RxAntennaPosition[i][j][2] = heightUE;
		}
	}

	AoA_LOS_dir = new double*[numberOfMobileStations];
	AoD_LOS_dir = new double*[numberOfMobileStations];
	ZoA_LOS_dir = new double*[numberOfMobileStations];
	ZoD_LOS_dir = new double*[numberOfMobileStations];
	double x_dir;
	double y_dir;
	vec cartLOS_BStoMS_Angle = zeros(3);
	vec cartLOS_MStoBS_Angle = zeros(3);
	vec sphLOSAngle;
		
	// Cycle through all mobile stations
    	for(int i = 0; i < numberOfMobileStations; i++){
		AoA_LOS_dir[i] = new double[neighbourPositions.size()];
		AoD_LOS_dir[i] = new double[neighbourPositions.size()];
		ZoA_LOS_dir[i] = new double[neighbourPositions.size()];
		ZoD_LOS_dir[i] = new double[neighbourPositions.size()];
		
		// Cycle through all base stations
		int PosIt = 1;
		for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
			if(it->first == bsId){
				x_dir = MSPos[i].x - xPos;
				y_dir = MSPos[i].y - yPos;

				cartLOS_MStoBS_Angle.set(0,x_dir);
				cartLOS_MStoBS_Angle.set(1,y_dir);
				cartLOS_MStoBS_Angle.set(2,heightUE);
				
				x_dir = xPos - MSPos[i].x;
				y_dir = yPos - MSPos[i].y;
				
				cartLOS_BStoMS_Angle.set(0,x_dir);
				cartLOS_BStoMS_Angle.set(1,y_dir);
				cartLOS_BStoMS_Angle.set(2,heightBS);
				
				sphLOSAngle = Cart_to_Sph(cartLOS_MStoBS_Angle);
				AoA_LOS_dir[i][0] = sphLOSAngle.get(0);
				ZoA_LOS_dir[i][0] = pi/2;//sphLOSAngle.get(1);
				
				sphLOSAngle = Cart_to_Sph(cartLOS_BStoMS_Angle);
				AoD_LOS_dir[i][0] = sphLOSAngle.get(0);
				ZoD_LOS_dir[i][0] = pi/2;//sphLOSAngle.get(1);
			}else{
				x_dir = MSPos[i].x - it->second.x;
				y_dir = MSPos[i].y - it->second.y;

				cartLOS_MStoBS_Angle.set(0,x_dir);
				cartLOS_MStoBS_Angle.set(1,y_dir);
				cartLOS_MStoBS_Angle.set(2,heightUE);
				
				x_dir = it->second.x - MSPos[i].x;
				y_dir = it->second.y - MSPos[i].y;
				
				cartLOS_BStoMS_Angle.set(0,x_dir);
				cartLOS_BStoMS_Angle.set(1,y_dir);
				cartLOS_BStoMS_Angle.set(2,heightBS);
				
				sphLOSAngle = Cart_to_Sph(cartLOS_MStoBS_Angle);
				AoA_LOS_dir[i][PosIt] = sphLOSAngle.get(0);
				ZoA_LOS_dir[i][PosIt] = pi/2;//sphLOSAngle.get(1);
				
				sphLOSAngle = Cart_to_Sph(cartLOS_BStoMS_Angle);
				AoD_LOS_dir[i][PosIt] = sphLOSAngle.get(0);
				ZoD_LOS_dir[i][PosIt] = pi/2;//sphLOSAngle.get(1);
				
				PosIt++;
			}
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
	
	SINRtable = new double**[numberOfMobileStations];
	for(int i = 0; i < numberOfMobileStations; i++){
		SINRtable[i] = new double*[timeSamples];
		for(int j = 0; j < timeSamples; j++){
			SINRtable[i][j] = new double[downRBs];
		}
	}
	
    	// Initialize BS antenna parameters
    	for(int i = 0; i < 3; i++){
		bs_antenna_downtilt[i] = 12.0;			// "..commonly assumed to be 12 Degrees.."; page 57 METIS Document
		bs_antenna_slant[i] = 0.0;				// "..usually Zero Degrees.."; page 57 METIS Document
	}
	bs_antenna_bearing[0] = 0.0;
	bs_antenna_bearing[1] = 120.0;
	bs_antenna_bearing[2] = 240.0;
	
	// Initialize MS antenna parameters
	// (The orientation may be completely randomly generated or based on more realistic distributions.)
    	for(int i = 0; i < numberOfMobileStations; i++){
		ms_antenna_bearing[i] = uniform(0.0,360.0);
		ms_antenna_downtilt[i] = uniform(0.0,360.0);
		ms_antenna_slant[i] = uniform(0.0,360.0);
	}
	
	//Assign LOS Conditions:
	double MS_BS_dist;
	for(int i = 0; i < numberOfMobileStations; i++){
		LOSCondition[i] = new bool[neighbourPositions.size()];
		int k = 1;
		for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
			if(it->first == bsId){
				MS_BS_dist = sqrt((xPos - MSPos[i].x)*(xPos - MSPos[i].x) + (yPos - MSPos[i].y)*(yPos - MSPos[i].y));
				LOSCondition[i][0] = LineOfSight(MS_BS_dist);
			}else{
				MS_BS_dist = sqrt((it->second.x - MSPos[i].x)*(it->second.x - MSPos[i].x) + (it->second.y - MSPos[i].y)*(it->second.y - MSPos[i].y));
				LOSCondition[i][k] = LineOfSight(MS_BS_dist);
				k++;
			}
		}
	}
	
	// Initialize large scale parameters (Line of Sight)
	// Generate Cross-correlation Matrix
	double a = module->par("cross_a_LOS"); // departure AS vs delay spread
	double b = module->par("cross_b_LOS"); // arrival AS vs delay spread
	double c = module->par("cross_c_LOS"); // arrival AS vs shadowing std
	double d = module->par("cross_d_LOS"); // departure AS vs shadoving std
	double e = module->par("cross_e_LOS"); // delay spread vs shadoving std
	double f = module->par("cross_f_LOS"); // departure AS vs arrival AS
	double g = module->par("cross_g_LOS"); // departure AS vs k-factor
	double h = module->par("cross_h_LOS"); // arrival AS vs k-factor
	double k = module->par("cross_k_LOS"); // delay spread vs k-factor
	double l = module->par("cross_l_LOS"); // shadowing std vs k-factor
	double m = module->par("cross_m_LOS"); // departure ZS vs shadowing std
	double n = module->par("cross_n_LOS"); // arrival ZS vs shadowing std
	double o = module->par("cross_o_LOS"); // departure ZS vs k-factor
	double p = module->par("cross_p_LOS"); // arrival ZS vs k-factor
	double q = module->par("cross_q_LOS"); // departure ZS vs delay spread
	double r = module->par("cross_r_LOS"); // arrival ZS vs delay spread
	double s = module->par("cross_s_LOS"); // departure ZS vs departure AS
	double t = module->par("cross_t_LOS"); // arrival ZS vs departure AS
	double u = module->par("cross_u_LOS"); // departure ZS vs arrival AS
	double v = module->par("cross_v_LOS"); // arrival ZS vs arrival AS
	double w = module->par("cross_w_LOS"); // arrival ZS vs departure ZS
		
	std::cout << "Before Auto Correlation!" << std::endl;
	
	// Generate Autocorrelation
	generateAutoCorrelation_LOS();
	
	std::cout << "After Auto Correlation!" << std::endl;
       
	//Cross-correlation matrix 
	/*
    	double cross_matrix[7][7] = {
		{1,a,b,q,r,e,k},
		{a,1,f,s,t,d,g},
		{b,f,1,u,v,c,h},
		{q,s,u,1,w,m,o},
		{r,t,v,w,1,n,p},
		{e,d,c,m,n,1,l},
		{k,g,h,o,p,l,1}};
	*/
	
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
    	a = module->par("DS_mu_LOS"); 			// mean of delay spread
	b = module->par("DS_eps_LOS"); 			// epsilon of delay spread
	c = module->par("AoD_mu_LOS"); 			// mean of AoD
	d = module->par("AoD_eps_LOS"); 		// epsilon of AoD
	e = module->par("AoA_mu_LOS"); 			// mean of AoA
	f = module->par("AoA_eps_LOS"); 		// epsilon of AoA
	g = module->par("ZoA_mu_LOS"); 			// mean of ZoA
	h = module->par("ZoA_eps_LOS"); 		// epsilon of ZoA
	k = module->par("SF_sigma_LOS"); 		// sigma for shadow fading
	l = module->par("K_mu"); 			// mean of Ricean K-factor
	double sigma_K = module->par("K_sigma");	// spread of K-factor
	
	sigma_ds_LOS = new double*[numberOfMobileStations];
	sigma_asD_LOS = new double*[numberOfMobileStations];
	sigma_asA_LOS = new double*[numberOfMobileStations];
	sigma_zsD_LOS = new double*[numberOfMobileStations];
	sigma_zsA_LOS = new double*[numberOfMobileStations];
	sigma_sf_LOS = new double*[numberOfMobileStations];
	sigma_kf_LOS = new double*[numberOfMobileStations];
		
	// TODO: Take Square root
	// cross_matrix = sqrt(cross_matrix)
	for(int i = 0; i < numberOfMobileStations; i++){
		sigma_ds_LOS[i] = new double[neighbourPositions.size()];
		sigma_asD_LOS[i] = new double[neighbourPositions.size()];
		sigma_asA_LOS[i] = new double[neighbourPositions.size()];
		sigma_zsD_LOS[i] = new double[neighbourPositions.size()];
		sigma_zsA_LOS[i] = new double[neighbourPositions.size()];
		sigma_sf_LOS[i] = new double[neighbourPositions.size()];
		sigma_kf_LOS[i] = new double[neighbourPositions.size()];
		int id_BS = 1;
		for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
			if(it->first == bsId){
				double ksi[7];
				for(int j = 0; j < 7;j++){
					ksi[j] = cross_matrix[j][0] * autoCorrelation_LOS[j][i] + cross_matrix[j][1] * autoCorrelation_LOS[j][i] + cross_matrix[j][2] * autoCorrelation_LOS[j][i] + cross_matrix[j][3] * autoCorrelation_LOS[j][i] + cross_matrix[j][4] * autoCorrelation_LOS[j][i] + cross_matrix[j][5] * autoCorrelation_LOS[j][i] + cross_matrix[j][6] * autoCorrelation_LOS[j][i];
				}
				sigma_ds_LOS[i][0]  = std::pow(10.0, (b*ksi[0] + a));      								// Log-Normal 
				sigma_asD_LOS[i][0] = std::min(104.0, std::pow(10.0, (d*ksi[1] + c)));      						// Log-Normal (maximum value should be 104 degrees) 
				sigma_asA_LOS[i][0] = std::min(104.0, std::pow(10.0, (f*ksi[2] + e)));      						// Log-Normal (maximum value should be 104 degrees) 
				sigma_zsA_LOS[i][0] = std::min(52.0, std::pow(10.0, (h*ksi[4] + g)));							// Log-Normal (maximum value should be 52 degrees) 
				dist2D = sqrt(pow((xPos - MSPos[i].x),2) + pow((yPos - MSPos[i].y),2));	
				sigma_zsD_LOS[i][0] = std::min(52.0, std::pow(10.0, (ksi[3] * sigma_ZSD(mean_ZSD(dist2D, heightUE, true), true))));	// Log-Normal (maximum value should be 52 degrees) 
				sigma_sf_LOS[i][0]  = std::pow(10.0, (0.1*k*ksi[5]));      								// Log-Normal dB
				sigma_kf_LOS[i][0]  = std::pow(10.0, (0.1*(sigma_K*ksi[6] + l)));	   						// Log-Normal dB
			}else{
				double ksi[7];
				for(int j = 0; j < 7;j++){
					ksi[j] = cross_matrix[j][0] * autoCorrelation_LOS[j][i] + cross_matrix[j][1] * autoCorrelation_LOS[j][i] + cross_matrix[j][2] * autoCorrelation_LOS[j][i] + cross_matrix[j][3] * autoCorrelation_LOS[j][i] + cross_matrix[j][4] * autoCorrelation_LOS[j][i] + cross_matrix[j][5] * autoCorrelation_LOS[j][i] + cross_matrix[j][6] * autoCorrelation_LOS[j][i];
				}
				sigma_ds_LOS[i][id_BS]  = std::pow(10.0, (b*ksi[0] + a));      								// Log-Normal 
				sigma_asD_LOS[i][id_BS] = std::min(104.0, std::pow(10.0, (d*ksi[1] + c)));      					// Log-Normal (maximum value should be 104 degrees) 
				sigma_asA_LOS[i][id_BS] = std::min(104.0, std::pow(10.0, (f*ksi[2] + e)));      					// Log-Normal (maximum value should be 104 degrees) 
				sigma_zsA_LOS[i][id_BS] = std::min(52.0, std::pow(10.0, (h*ksi[4] + g)));						// Log-Normal (maximum value should be 52 degrees) 
				dist2D = sqrt(pow((it->second.x - MSPos[i].x),2) + pow((it->second.y - MSPos[i].y),2));
				sigma_zsD_LOS[i][id_BS] = std::min(52.0, std::pow(10.0, (ksi[3] * sigma_ZSD(mean_ZSD(dist2D, heightUE, true), true))));	// Log-Normal (maximum value should be 52 degrees) 
				sigma_sf_LOS[i][id_BS]  = std::pow(10.0, (0.1*k*ksi[5]));      								// Log-Normal dB
				sigma_kf_LOS[i][id_BS]  = std::pow(10.0, (0.1*(sigma_K*ksi[6] + l)));	   						// Log-Normal dB
				id_BS++;
			}
		}
	}
    
    	// Initialize large scale parameters (Non Line of Sight)
    	// Generate Cross-correlation Matrix
	a = module->par("cross_a"); // departure AS vs delay spread
	b = module->par("cross_b"); // arrival AS vs delay spread
	c = module->par("cross_c"); // arrival AS vs shadowing std
	d = module->par("cross_d"); // departure AS vs shadoving std
	e = module->par("cross_e"); // delay spread vs shadoving std
	f = module->par("cross_f"); // departure AS vs arrival AS
	m = module->par("cross_m"); // departure ZS vs shadowing std
	n = module->par("cross_n"); // arrival ZS vs shadowing std
	q = module->par("cross_q"); // departure ZS vs delay spread
	r = module->par("cross_r"); // arrival ZS vs delay spread
	s = module->par("cross_s"); // departure ZS vs departure AS
	t = module->par("cross_t"); // arrival ZS vs departure AS
	u = module->par("cross_u"); // departure ZS vs arrival AS
	v = module->par("cross_v"); // arrival ZS vs arrival AS
	w = module->par("cross_w"); // arrival ZS vs departure ZS

	std::cout << "Before Auto Correlation (NLOS)!" << std::endl;
	
	// Generate Autocorrelation
	generateAutoCorrelation_NLOS();
	
	std::cout << "After Auto Correlation (NLOS)!" << std::endl;
       
	//Cross-correlation matrix 
	/*
    	double cross_matrix2[6][6] = {
		{1,a,b,q,r,e},
		{a,1,f,s,t,d},
		{b,f,1,u,v,c},
		{q,s,u,1,w,m},
		{r,t,v,w,1,n},
		{e,d,c,m,n,1}};
	*/
		
	//TODO: Remove hardcoded, values taken after taking the square-root of the cross_matrix using sqrtm() in MATLAB 
	double cross_matrix2[6][6] = {
		{ 0.8302,   0.0709,   0.1949,  -0.3252,  -0.0341,  -0.4011},
 		{ 0.0709,   0.9079,  -0.0277,   0.3008,   0.2808,   0.0259},
  		{ 0.1949,  -0.0277,   0.9580,   0.0352,   0.1131,  -0.1716},
  		{-0.3252,   0.3008,   0.0352,   0.8911,  -0.0542,  -0.0741},
  		{-0.0341,   0.2808,   0.1131,  -0.0542,   0.9509,  -0.0030},
  		{-0.4011,   0.0259,  -0.1716,  -0.0741,  -0.0030,   0.8964}};
		
    	// Transform Normal distributed random numbers to scenario specific distributions
    	a = module->par("DS_mu_NLOS"); 				// departure AS vs delay spread
	b = module->par("DS_eps_NLOS"); 			// arrival AS vs delay spread
	c = module->par("AoD_mu_NLOS"); 			// arrival AS vs shadowing std
	d = module->par("AoD_eps_NLOS"); 			// departure AS vs shadoving std
	e = module->par("AoA_mu_NLOS"); 			// delay spread vs shadoving std
	f = module->par("AoA_eps_NLOS"); 			// departure AS vs arrival AS
	g = module->par("SF_sigma_NLOS"); 			// departure AS vs k-factor
	h = module->par("ZoA_mu_NLOS"); 			// arrival AS vs k-factor
	k = module->par("ZoA_eps_NLOS"); 			// delay spread vs k-factor	
	
	sigma_ds_NLOS = new double*[numberOfMobileStations];
	sigma_asD_NLOS = new double*[numberOfMobileStations];
	sigma_asA_NLOS = new double*[numberOfMobileStations];
	sigma_zsD_NLOS = new double*[numberOfMobileStations];
	sigma_zsA_NLOS = new double*[numberOfMobileStations];
	sigma_sf_NLOS = new double*[numberOfMobileStations];

	// TODO: Take Square root
	// cross_matrix = sqrt(cross_matrix)
	for(int i = 0; i < numberOfMobileStations; i++){
		sigma_ds_NLOS[i] = new double[neighbourPositions.size()];
		sigma_asD_NLOS[i] = new double[neighbourPositions.size()];
		sigma_asA_NLOS[i] = new double[neighbourPositions.size()];
		sigma_zsD_NLOS[i] = new double[neighbourPositions.size()];
		sigma_zsA_NLOS[i] = new double[neighbourPositions.size()];
		sigma_sf_NLOS[i] = new double[neighbourPositions.size()];
		int it_BS = 1;
		for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
			if(it->first == bsId){
				double ksi[6];
				for(int j = 0;j < 6;j++){
					ksi[j] = cross_matrix2[j][0] * autoCorrelation_NLOS[j][i] + cross_matrix2[j][1] * autoCorrelation_NLOS[j][i] + cross_matrix2[j][2] * autoCorrelation_NLOS[j][i] + cross_matrix2[j][3] * autoCorrelation_NLOS[j][i] + cross_matrix2[j][4] * autoCorrelation_NLOS[j][i] + cross_matrix2[j][5] * autoCorrelation_NLOS[j][i];
				}
				sigma_ds_NLOS[i][0]  = std::pow(10.0, (b*ksi[0] + a));      								// Log-Normal 
				sigma_asD_NLOS[i][0] = std::min(104.0, std::pow(10.0, (d*ksi[1] + c)));      						// Log-Normal (maximum value should be 104 degrees)  
				sigma_asA_NLOS[i][0] = std::min(104.0, std::pow(10.0, (f*ksi[2] + e)));      						// Log-Normal (maximum value should be 104 degrees) 
				sigma_zsA_NLOS[i][0] = std::min(52.0, std::pow(10.0, (k*ksi[4] + h)));							// Log-Normal (maximum value should be 52 degrees) 
				dist2D = sqrt(pow((xPos - MSPos[i].x),2) + pow((yPos - MSPos[i].y),2));
				sigma_zsD_NLOS[i][0] = std::min(52.0, std::pow(10.0, (ksi[3] * sigma_ZSD(mean_ZSD(dist2D, heightUE, false), false))));	// Log-Normal (maximum value should be 52 degrees) 
				sigma_sf_NLOS[i][0]  = std::pow(10.0, (0.1*g*ksi[5]));      								// Log-Normal dB

				//sigma_ds_NLOS[i]  = 7.6e-8;
				//sigma_asD_NLOS[i] = 15;
				//sigma_asA_NLOS[i] = 35;
				//sigma_sf_NLOS[i]  = 1.5;
			}else{
				double ksi[6];
				for(int j = 0;j < 6;j++){
					ksi[j] = cross_matrix2[j][0] * autoCorrelation_NLOS[j][i] + cross_matrix2[j][1] * autoCorrelation_NLOS[j][i] + cross_matrix2[j][2] * autoCorrelation_NLOS[j][i] + cross_matrix2[j][3] * autoCorrelation_NLOS[j][i] + cross_matrix2[j][4] * autoCorrelation_NLOS[j][i] + cross_matrix2[j][5] * autoCorrelation_NLOS[j][i];
				}
				sigma_ds_NLOS[i][it_BS]  = std::pow(10.0, (b*ksi[0] + a));      							// Log-Normal 
				sigma_asD_NLOS[i][it_BS] = std::min(104.0, std::pow(10.0, (d*ksi[1] + c)));      					// Log-Normal (maximum value should be 104 degrees) 
				sigma_asA_NLOS[i][it_BS] = std::min(104.0, std::pow(10.0, (f*ksi[2] + e)));      					// Log-Normal (maximum value should be 104 degrees) 
				sigma_zsA_NLOS[i][it_BS] = std::min(52.0, std::pow(10.0, (k*ksi[4] + h)));						// Log-Normal (maximum value should be 52 degrees) 
				dist2D = sqrt(pow((it->second.x - MSPos[i].x),2) + pow((it->second.y - MSPos[i].y),2));
				sigma_zsD_NLOS[i][it_BS] = std::min(52.0, std::pow(10.0, (ksi[3] * sigma_ZSD(mean_ZSD(dist2D, heightUE, false), false))));// Log-Normal (maximum value should be 52 degrees) 
				sigma_sf_NLOS[i][it_BS]  = std::pow(10.0, (0.1*g*ksi[5]));      							// Log-Normal dB
				
				//sigma_ds_NLOS[i]  = 7.6e-8;
				//sigma_asD_NLOS[i] = 15;
				//sigma_asA_NLOS[i] = 35;
				//sigma_sf_NLOS[i]  = 1.5;
				//sigma_kf_NLOS[i]  = 1.0;
				it_BS++;
			}
		}
	} 	
	
	std::cout << "Finished Large Scale parameter.." << std::endl;

    	// Begin small scale parameter generation.
    
    	// Generate delays for each cluster according to Formula: 7:38 (METIS Document)
	clusterDelays = new double**[numberOfMobileStations];
	clusterDelays_LOS = new double**[numberOfMobileStations];
	double delayScaling;
	int N_cluster_LOS = module->par("NumberOfClusters_LOS");
	int N_cluster_NLOS = module->par("NumberOfClusters_NLOS");
	std::cout << "N_cluster_LOS: " << N_cluster_LOS << std::endl;
    	std::cout << "N_cluster_NLOS: " << N_cluster_NLOS << std::endl;
    
    	for(int i = 0; i < numberOfMobileStations; i++){
		clusterDelays[i] = new double*[neighbourPositions.size()];
		clusterDelays_LOS[i] = new double*[neighbourPositions.size()];
		int k = 1;
		for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
			if(it->first == bsId){
				if(LOSCondition[i][0]){
					// LOS Condition
					clusterDelays[i][0] = new double[N_cluster_LOS];
					clusterDelays_LOS[i][0] = new double[N_cluster_LOS];
					delayScaling = module->par("DelayScaling_LOS");
					
					for(int j = 0; j < N_cluster_LOS; j++){
						clusterDelays[i][0][j] = -1.0*delayScaling*sigma_ds_LOS[i][0]*log(uniform(0,1));
					}
					
					double min_delay = *std::min_element(clusterDelays[i][0], clusterDelays[i][0] + N_cluster_LOS);
					// Normalize the delays (7.39)
					for(int j = 0; j < N_cluster_LOS; j++){
						clusterDelays[i][0][j] = clusterDelays[i][0][j] - min_delay;
					}

					// Sort the delays (7.39)
					std::sort(clusterDelays[i][0], clusterDelays[i][0] + N_cluster_LOS, std::less<double>());
					
					// Compute LOS Peak compensation factor (7.41)
					double K = 10.0 * log10(abs(sigma_kf_LOS[i][0]));
					double C_DS = 0.7705 - 0.0433 * K + 0.0002 * pow(K,2) + 0.000017 * pow(K,3);
					
					// Apply LOS compensation factor
					for(int j = 0; j < N_cluster_LOS; j++){
						clusterDelays_LOS[i][0][j] = clusterDelays[i][0][j] / C_DS;
					}
					
				}else{
					// NLOS Condition
					clusterDelays[i][0] = new double[N_cluster_NLOS];
					delayScaling = module->par("DelayScaling_NLOS");
					
					for(int j = 0; j < N_cluster_NLOS; j++){
						clusterDelays[i][0][j] = -1.0*delayScaling*sigma_ds_NLOS[i][0]*log(uniform(0,1));
					}
					
					double min_delay = *std::min_element(clusterDelays[i][0], clusterDelays[i][0] + N_cluster_NLOS);
					// Normalize the delays (7.39)
					for(int j = 0; j < N_cluster_NLOS; j++){
						clusterDelays[i][0][j] = clusterDelays[i][0][j] - min_delay;
						//std::cout << "Delay: " << clusterDelays[i][0][j] << std::endl;
					}

					// Sort the delays (7.39)
					std::sort(clusterDelays[i][0], clusterDelays[i][0] + N_cluster_NLOS,std::less<double>());
					
				}
			}else{
				if(LOSCondition[i][k]){
					// LOS Condition
					clusterDelays[i][k] = new double[N_cluster_LOS];
					clusterDelays_LOS[i][k] = new double[N_cluster_LOS];
					delayScaling = module->par("DelayScaling_LOS");
					
					for(int j = 0; j < N_cluster_LOS; j++){
						clusterDelays[i][k][j] = -1.0*delayScaling*sigma_ds_LOS[i][k]*log(uniform(0,1));
					}
					
					double min_delay = *std::min_element(clusterDelays[i][k], clusterDelays[i][k] + N_cluster_LOS);
					// Normalize the delays (7.39)
					for(int j = 0; j < N_cluster_LOS; j++){
						clusterDelays[i][k][j] = clusterDelays[i][k][j] - min_delay;
					}
					
					// Sort the delays (7.39)
					std::sort(clusterDelays[i][k], clusterDelays[i][k] + N_cluster_LOS, std::less<double>());
					
					// Compute LOS Peak compensation factor (7.41)
					double K = 10.0 * log10(abs(sigma_kf_LOS[i][k]));
					double C_DS = 0.7705 - 0.0433 * K + 0.0002 * pow(K,2) + 0.000017 * pow(K,3);
					
					// Apply LOS compensation factor
					for(int j = 0; j < N_cluster_LOS; j++){
						clusterDelays_LOS[i][k][j] = clusterDelays[i][k][j] / C_DS;
					}	
				}else{
					// NLOS Condition
					clusterDelays[i][k] = new double[N_cluster_NLOS];
					delayScaling = module->par("DelayScaling_NLOS");
					
					for(int j = 0; j < N_cluster_NLOS; j++){
						clusterDelays[i][k][j] = -1.0*delayScaling*sigma_ds_NLOS[i][k]*log(uniform(0,1));
					}
					
					double min_delay = *std::min_element(clusterDelays[i][k], clusterDelays[i][k] + N_cluster_NLOS);
					// Normalize the delays (7.39)
					for(int j = 0; j < N_cluster_NLOS; j++){
						clusterDelays[i][k][j] = clusterDelays[i][k][j] - min_delay;
						//std::cout << "Delay: " << clusterDelays[i][k][j] << std::endl;
					}

					// Sort the delays (7.39)
					std::sort(clusterDelays[i][k], clusterDelays[i][k] + N_cluster_NLOS,std::less<double>());
				}
				k++;
			}
		}
	}
		
	clusterPowers = new double**[numberOfMobileStations];
	double cluster_shadowing;
	double sum;
	
	// Generate cluster powers. 
	for(int i = 0; i < numberOfMobileStations; i++){
		clusterPowers[i] = new double*[neighbourPositions.size()];
		int itIdx = 1;
		for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
			if(it->first == bsId){
				if(LOSCondition[i][0]){
					//std::cout << "LOS Case!" << std::endl;
					clusterPowers[i][0] = new double[N_cluster_LOS];
					sum = 0;
					double K = sigma_kf_LOS[i][0];				// K-factor in linear scale
					
					double P1_LOS = K / (K + 1);
					
					double temp_CS = module->par("PerClusterShadowing_LOS");
					cluster_shadowing = pow(temp_CS,2);
					delayScaling = module->par("DelayScaling_LOS");
					for(int j = 0; j < N_cluster_LOS; j++){
						// Formula 7.42
						clusterPowers[i][0][j] = pow(10,-1.0*normal(0,cluster_shadowing)/10) * exp(-1.0 * clusterDelays[i][0][j] * ((delayScaling - 1) / (delayScaling * sigma_ds_LOS[i][0])) );
						sum += clusterPowers[i][0][j];
					}
					
					for(int j = 0; j < N_cluster_LOS; j++){
						if (j==0){
                                                        clusterPowers[i][0][j] = ((1/(K+1))*clusterPowers[i][0][j]/sum) + P1_LOS;
                                                }else{
						        clusterPowers[i][0][j] = (1/(K+1))*clusterPowers[i][0][j]/sum;
						}
					}
				}else{
					//std::cout << "NLOS Case!" << std::endl;
					clusterPowers[i][0] = new double[N_cluster_NLOS];
					sum = 0;
					
					double temp_CS = module->par("PerClusterShadowing_NLOS");
					cluster_shadowing = pow(temp_CS,2); //METIS
					delayScaling = module->par("DelayScaling_NLOS");
					for(int j = 0; j < N_cluster_NLOS; j++){
						// Formula 7.42
						clusterPowers[i][0][j] = pow(10,-1.0*normal(0,cluster_shadowing)/10) * exp(-1.0 * clusterDelays[i][0][j] * ((delayScaling - 1) / (delayScaling * sigma_ds_NLOS[i][0])) ); //METIS
						//cluster_shadowing = normal(0,1) * temp_CS; //WINNER II
						//cluster_shadowing = 2.0;
						//clusterDelays[i][j] = 7.6e-8;
						//clusterPowers[i][j] = pow(10,-1.0*cluster_shadowing/10) * exp(-1.0 * clusterDelays[i][j] * ((delayScaling - 1) / (delayScaling * sigma_ds_NLOS[i])) ); //WINNER II
						sum += clusterPowers[i][0][j];
					}
					
					for(int j = 0; j < N_cluster_NLOS; j++){
						clusterPowers[i][0][j] /= sum;
					}
				}
			}else{
				if(LOSCondition[i][itIdx]){
					//std::cout << "LOS Case!" << std::endl;
					clusterPowers[i][itIdx] = new double[N_cluster_LOS];
					sum = 0;
					double K = sigma_kf_LOS[i][itIdx]; 			// K-factor in linear scale
					
					double P1_LOS = K / (K + 1);
					
					double temp_CS = module->par("PerClusterShadowing_LOS");
					cluster_shadowing = pow(temp_CS,2);
					delayScaling = module->par("DelayScaling_LOS");
					for(int j = 0; j < N_cluster_LOS; j++){
						// Formula 7.42
						clusterPowers[i][itIdx][j] = pow(10,-1.0*normal(0,cluster_shadowing)/10) * exp(-1.0 * clusterDelays[i][itIdx][j] * ((delayScaling - 1) / (delayScaling * sigma_ds_LOS[i][itIdx])) );
						sum += clusterPowers[i][itIdx][j];
					}
					
					for(int j = 0; j < N_cluster_LOS; j++){
						if (j==0){
                                                        clusterPowers[i][itIdx][j] = ((1/(K+1))*clusterPowers[i][itIdx][j]/sum) + P1_LOS;
                                                }else{
						        clusterPowers[i][itIdx][j] = (1/(K+1))*clusterPowers[i][itIdx][j]/sum;
						}
					}
				}else{
					//std::cout << "NLOS Case!" << std::endl;
					clusterPowers[i][itIdx] = new double[N_cluster_NLOS];
					sum = 0;
					
					double temp_CS = module->par("PerClusterShadowing_NLOS");
					cluster_shadowing = pow(temp_CS,2); //METIS
					delayScaling = module->par("DelayScaling_NLOS");
					for(int j = 0; j < N_cluster_NLOS; j++){
						// Formula 7.42
						clusterPowers[i][itIdx][j] = pow(10,-1.0*normal(0,cluster_shadowing)/10) * exp(-1.0 * clusterDelays[i][itIdx][j] * ((delayScaling - 1) / (delayScaling * sigma_ds_NLOS[i][itIdx])) ); //METIS
						//cluster_shadowing = normal(0,1) * temp_CS; //WINNER II
						//cluster_shadowing = 2.0;
						//clusterDelays[i][j] = 7.6e-8;
						//clusterPowers[i][j] = pow(10,-1.0*cluster_shadowing/10) * exp(-1.0 * clusterDelays[i][j] * ((delayScaling - 1) / (delayScaling * sigma_ds_NLOS[i])) ); //WINNER II
						sum += clusterPowers[i][itIdx][j];
					}
					
					for(int j = 0; j < N_cluster_NLOS; j++){
						clusterPowers[i][itIdx][j] /= sum;
					}
				}
				itIdx++;
			}
		}
	}
		
	int numOfRays_LOS = module->par("NumberOfRays_LOS");
	int numOfRays_NLOS = module->par("NumberOfRays_NLOS");
	
	// Precompute powers per ray (7.46)
	rayPowers = new double**[numberOfMobileStations];
	for(int i = 0; i < numberOfMobileStations; i++){
		rayPowers[i] = new double*[neighbourPositions.size()];
		int itIdx = 1;
		for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
			if(it->first == bsId){
				if(LOSCondition[i][0]){
					rayPowers[i][0] = new double[N_cluster_LOS];
					for(int j = 0; j < N_cluster_LOS; j++){
						rayPowers[i][0][j] = clusterPowers[i][0][j] / numOfRays_LOS;
					}
				}else{
					rayPowers[i][0] = new double[N_cluster_NLOS];
					for(int j = 0; j < N_cluster_NLOS; j++){
						rayPowers[i][0][j] = clusterPowers[i][0][j] / numOfRays_NLOS;
					}
				}
			}else{
				if(LOSCondition[i][itIdx]){
					rayPowers[i][itIdx] = new double[N_cluster_LOS];
					for(int j = 0; j < N_cluster_LOS; j++){
						rayPowers[i][itIdx][j] = clusterPowers[i][itIdx][j] / numOfRays_LOS;
					}
				}else{
					rayPowers[i][itIdx] = new double[N_cluster_NLOS];
					for(int j = 0; j < N_cluster_NLOS; j++){
						rayPowers[i][itIdx][j] = clusterPowers[i][itIdx][j] / numOfRays_NLOS;
					}
				}
				itIdx++;
			}
		}
	}
	
	// Ray offset. Table 7.6
	double ray_offset[20] = {	
							0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492, 
							0.3715, -0.3715, 0.5129, -0.5129, 0.6797, -0.6797,
							0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195,
							2.1551, -2.1551
							};
	
	double **azimuth_cluster_ASA;
	azimuth_cluster_ASA = new double*[numberOfMobileStations];
	azimuth_ASA = new double***[numberOfMobileStations];
	
	// Generate azimuth angles of arrival
	for(int i = 0; i < numberOfMobileStations; i++){
		azimuth_ASA[i] = new double**[neighbourPositions.size()];
		int itIdx = 1;
		for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
			if(it->first == bsId){
				if(LOSCondition[i][0]){
					azimuth_ASA[i][0] = new double*[N_cluster_LOS];
					azimuth_cluster_ASA[i] = new double[N_cluster_LOS];
			
					double X_1 = uniform(-1.0,1.0);
					double Y_1 = normal(0.0, (pow(sigma_asA_LOS[i][0] / 7.0,2)));
			
					for(int j = 0; j < N_cluster_LOS; j++){
						azimuth_ASA[i][0][j] = new double[numOfRays_LOS];
				
						//Formula 7.47
						azimuth_cluster_ASA[i][j] = (sigma_asA_LOS[i][0] / (0.7 * C_AS(N_cluster_LOS, true, i))) * sqrt( -1.0 * log(clusterPowers[i][0][j] / *(std::max_element(clusterPowers[i][0], clusterPowers[i][0] + N_cluster_LOS))) );
						double X_N = uniform(-1.0,1.0);
						double Y_N = normal(0.0, (pow(sigma_asA_LOS[i][0] / 7.0,2)));
				
						// First Ray has geometric LOS direction?!
						azimuth_cluster_ASA[i][j] = azimuth_cluster_ASA[i][j] * (X_N - X_1) + (Y_N - Y_1) + (AoA_LOS_dir[i][0]*180.0/pi);  // since AOA_LOS_dir is in radians
				
						for(int k = 0; k < numOfRays_LOS; k++){
							double cluster_ASA = module->par("Cluster_ASA_LOS");
					
							// Final angle computation per ray (LOS)
							azimuth_ASA[i][0][j][k] = azimuth_cluster_ASA[i][j] + cluster_ASA * ray_offset[k];

							if(azimuth_ASA[i][0][j][k] < -360.0){
								azimuth_ASA[i][0][j][k]+=360.0;
							}
						}
					}
				}else{
					azimuth_ASA[i][0] = new double*[N_cluster_NLOS];
					azimuth_cluster_ASA[i] = new double[N_cluster_NLOS];
			
					for(int j = 0; j < N_cluster_NLOS; j++){
						azimuth_ASA[i][0][j] = new double[numOfRays_NLOS];
				
						// Formula 7.47
						azimuth_cluster_ASA[i][j] = (sigma_asA_NLOS[i][0] / (0.7 * C_AS(N_cluster_NLOS, false, i))) * sqrt( -1.0 * log(clusterPowers[i][0][j] / *(std::max_element(clusterPowers[i][0], clusterPowers[i][0] + N_cluster_NLOS))) );
						double X_N = uniform(-1.0,1.0);
						double Y_N = normal(0.0, (pow(sigma_asA_NLOS[i][0] / 7.0,2)));
				
						azimuth_cluster_ASA[i][j] = azimuth_cluster_ASA[i][j] * X_N + Y_N + (AoA_LOS_dir[i][0]*180.0/pi);   // since AOA_LOS_dir is in radians
				 
						for(int k = 0; k < numOfRays_NLOS; k++){
							double cluster_ASA = module->par("Cluster_ASA_NLOS");
					
							// Final angle computation per ray (NLOS)
							azimuth_ASA[i][0][j][k] = azimuth_cluster_ASA[i][j] + cluster_ASA * ray_offset[k];
						}
					}
				}
			}else{
				if(LOSCondition[i][itIdx]){
					azimuth_ASA[i][itIdx] = new double*[N_cluster_LOS];
					azimuth_cluster_ASA[i] = new double[N_cluster_LOS];
			
					double X_1 = uniform(-1.0,1.0);
					double Y_1 = normal(0.0, (pow(sigma_asA_LOS[i][itIdx] / 7.0,2)));
			
					for(int j = 0; j < N_cluster_LOS; j++){
						azimuth_ASA[i][itIdx][j] = new double[numOfRays_LOS];
				
						//Formula 7.47
						azimuth_cluster_ASA[i][j] = (sigma_asA_LOS[i][itIdx] / (0.7 * C_AS(N_cluster_LOS, true, i))) * sqrt( -1.0 * log(clusterPowers[i][itIdx][j] / *(std::max_element(clusterPowers[i][itIdx], clusterPowers[i][itIdx] + N_cluster_LOS))) );
						double X_N = uniform(-1.0,1.0);
						double Y_N = normal(0.0, (pow(sigma_asA_LOS[i][itIdx] / 7.0,2)));
				
						// First Ray has geometric LOS direction?!
						azimuth_cluster_ASA[i][j] = azimuth_cluster_ASA[i][j] * (X_N - X_1) + (Y_N - Y_1) + (AoA_LOS_dir[i][itIdx]*180.0/pi);   // since AOA_LOS_dir is in radians
  				
						for(int k = 0; k < numOfRays_LOS; k++){
							double cluster_ASA = module->par("Cluster_ASA_LOS");
					
							// Final angle computation per ray (LOS)
							azimuth_ASA[i][itIdx][j][k] = azimuth_cluster_ASA[i][j] + cluster_ASA * ray_offset[k];
						}
					}
					
				}else{
					azimuth_ASA[i][itIdx] = new double*[N_cluster_NLOS];
					azimuth_cluster_ASA[i] = new double[N_cluster_NLOS];
			
					for(int j = 0; j < N_cluster_NLOS; j++){
						azimuth_ASA[i][itIdx][j] = new double[numOfRays_NLOS];
				
						// Formula 7.47
						azimuth_cluster_ASA[i][j] = (sigma_asA_NLOS[i][itIdx] / (0.7 * C_AS(N_cluster_NLOS, false, i))) * sqrt( -1.0 * log(clusterPowers[i][itIdx][j] / *(std::max_element(clusterPowers[i][itIdx], clusterPowers[i][itIdx] + N_cluster_NLOS))) );
						double X_N = uniform(-1.0,1.0);
						double Y_N = normal(0.0, (pow(sigma_asA_NLOS[i][itIdx] / 7.0,2)));
				
						azimuth_cluster_ASA[i][j] = azimuth_cluster_ASA[i][j] * X_N + Y_N + (AoA_LOS_dir[i][itIdx]*180.0/pi);   // since AOA_LOS_dir is in radians
				
						for(int k = 0; k < numOfRays_NLOS; k++){
							double cluster_ASA = module->par("Cluster_ASA_NLOS");
					
							// Final angle computation per ray (NLOS)
							azimuth_ASA[i][itIdx][j][k] = azimuth_cluster_ASA[i][j] + cluster_ASA * ray_offset[k];
						}
					}
				}
				itIdx++;
			}
		}
	}
	
	double **azimuth_cluster_ASD;
	azimuth_cluster_ASD = new double*[numberOfMobileStations];
	azimuth_ASD = new double***[numberOfMobileStations];
	
	// Generate azimuth angles of departure in the same way 
	for(int i = 0; i < numberOfMobileStations; i++){
		azimuth_ASD[i] = new double**[neighbourPositions.size()];
		int itIdx = 1;
		for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
			if(it->first == bsId){
				if(LOSCondition[i][0]){
					azimuth_ASD[i][0] = new double*[N_cluster_LOS];
					azimuth_cluster_ASD[i] = new double[N_cluster_LOS];
			
					double X_1 = uniform(-1.0,1.0);
					double Y_1 = normal(0.0, (pow(sigma_asD_LOS[i][0] / 7.0,2)));
			
					for(int j = 0; j < N_cluster_LOS; j++){
						azimuth_ASD[i][0][j] = new double[numOfRays_LOS];
				
						//Formula 7.47
						azimuth_cluster_ASD[i][j] = (sigma_asD_LOS[i][0] / (0.7 * C_AS(N_cluster_LOS, true, i))) * sqrt( -1.0 * log(clusterPowers[i][0][j] / *(std::max_element(clusterPowers[i][0], clusterPowers[i][0] + N_cluster_LOS))) );
						double X_N = uniform(-1.0,1.0);
						double Y_N = normal(0.0, (pow(sigma_asD_LOS[i][0] / 7.0,2)));
				
						// First Ray has geometric LOS direction?!
						azimuth_cluster_ASD[i][j] = azimuth_cluster_ASD[i][j] * (X_N - X_1) + (Y_N - Y_1) + (AoD_LOS_dir[i][0]*180.0/pi);   // since AOD_LOS_dir is in radians
				
						for(int k = 0; k < numOfRays_LOS; k++){
							double cluster_ASD = module->par("Cluster_ASD_LOS");
					
							// Final angle computation per ray (LOS)
							azimuth_ASD[i][0][j][k] = azimuth_cluster_ASD[i][j] + cluster_ASD * ray_offset[k];

							if(azimuth_ASD[i][0][j][k] < -360.0){
								azimuth_ASD[i][0][j][k]+=360.0;
							}
						}
					}
				}else{
					azimuth_ASD[i][0] = new double*[N_cluster_NLOS];
					azimuth_cluster_ASD[i] = new double[N_cluster_NLOS];
			
					for(int j = 0; j < N_cluster_NLOS; j++){
						azimuth_ASD[i][0][j] = new double[numOfRays_NLOS];
				
						// Formula 7.47
						azimuth_cluster_ASD[i][j] = (sigma_asD_NLOS[i][0] / (0.7 * C_AS(N_cluster_NLOS, false, i))) * sqrt( -1.0 * log(clusterPowers[i][0][j] / *(std::max_element(clusterPowers[i][0], clusterPowers[i][0] + N_cluster_NLOS))) );
						double X_N = uniform(-1.0,1.0);
						double Y_N = normal(0.0, (pow(sigma_asD_NLOS[i][0] / 7.0,2)));
				
						azimuth_cluster_ASD[i][j] = azimuth_cluster_ASD[i][j] * X_N + Y_N + (AoD_LOS_dir[i][0]*180.0/pi);    // since AOD_LOS_dir is in radians
				
						for(int k = 0; k < numOfRays_NLOS; k++){
							double cluster_ASD = module->par("Cluster_ASD_NLOS");
					
							// Final angle computation per ray (NLOS)
							azimuth_ASD[i][0][j][k] = azimuth_cluster_ASD[i][j] + cluster_ASD * ray_offset[k];
						}
					}
				}
			}else{
				if(LOSCondition[i][itIdx]){
					azimuth_ASD[i][itIdx] = new double*[N_cluster_LOS];
					azimuth_cluster_ASD[i] = new double[N_cluster_LOS];
			
					double X_1 = uniform(-1.0,1.0);
					double Y_1 = normal(0.0, (pow(sigma_asD_LOS[i][itIdx] / 7.0,2)));
			
					for(int j = 0; j < N_cluster_LOS; j++){
						azimuth_ASD[i][itIdx][j] = new double[numOfRays_LOS];
				
						//Formula 7.47
						azimuth_cluster_ASD[i][j] = (sigma_asD_LOS[i][itIdx] / (0.7 * C_AS(N_cluster_LOS, true, i))) * sqrt( -1.0 * log(clusterPowers[i][itIdx][j] / *(std::max_element(clusterPowers[i][itIdx], clusterPowers[i][itIdx] + N_cluster_LOS))) );
						double X_N = uniform(-1.0,1.0);
						double Y_N = normal(0.0, (pow(sigma_asD_LOS[i][itIdx] / 7.0,2)));
				
						// First Ray has geometric LOS direction?!
						azimuth_cluster_ASD[i][j] = azimuth_cluster_ASD[i][j] * (X_N - X_1) + (Y_N - Y_1) + (AoD_LOS_dir[i][itIdx]*180.0/pi);   // since AOD_LOS_dir is in radians
				
						for(int k = 0; k < numOfRays_LOS; k++){
							double cluster_ASD = module->par("Cluster_ASD_LOS");
					
							// Final angle computation per ray (LOS)
							azimuth_ASD[i][itIdx][j][k] = azimuth_cluster_ASD[i][j] + cluster_ASD * ray_offset[k];
						}
					}
					
				}else{
					azimuth_ASD[i][itIdx] = new double*[N_cluster_NLOS];
					azimuth_cluster_ASD[i] = new double[N_cluster_NLOS];
			
					for(int j = 0; j < N_cluster_NLOS; j++){
						azimuth_ASD[i][itIdx][j] = new double[numOfRays_NLOS];
				
						// Formula 7.47
						azimuth_cluster_ASD[i][j] = (sigma_asD_NLOS[i][itIdx] / (0.7 * C_AS(N_cluster_NLOS, false, i))) * sqrt( -1.0 * log(clusterPowers[i][itIdx][j] / *(std::max_element(clusterPowers[i][itIdx], clusterPowers[i][itIdx] + N_cluster_NLOS))) );
						double X_N = uniform(-1.0,1.0);
						double Y_N = normal(0.0, (pow(sigma_asD_NLOS[i][itIdx] / 7.0,2)));
				
						azimuth_cluster_ASD[i][j] = azimuth_cluster_ASD[i][j] * X_N + Y_N + (AoD_LOS_dir[i][itIdx]*180.0/pi);    // since AOD_LOS_dir is in radians
				
						for(int k = 0; k < numOfRays_NLOS; k++){
							double cluster_ASD = module->par("Cluster_ASD_NLOS");
					
							// Final angle computation per ray (NLOS)
							azimuth_ASD[i][itIdx][j][k] = azimuth_cluster_ASD[i][j] + cluster_ASD * ray_offset[k];
						}
					}
				}
				itIdx++;
			}
		}
	}
	
	// Generate Elevation angles (for 2D channel model, set all elevtaion angles to 90 deg; otherwise use the code below for elevation angle generation according to METIS 1.2 or 1.4)
	elevation_ASA = new double***[numberOfMobileStations];
	elevation_ASD = new double***[numberOfMobileStations];
	for(int i = 0; i < numberOfMobileStations; i++){
		elevation_ASA[i] = new double**[neighbourPositions.size()];
		elevation_ASD[i] = new double**[neighbourPositions.size()];
		int itIdx = 1;
		for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
			if(it->first == bsId){
				elevation_ASA[i][0] = new double*[N_cluster_NLOS];
				elevation_ASD[i][0] = new double*[N_cluster_NLOS];
				for(int j = 0; j < N_cluster_NLOS; j++){
					elevation_ASA[i][0][j] = new double[numOfRays_NLOS];
					elevation_ASD[i][0][j] = new double[numOfRays_NLOS];
					for(int k = 0; k < numOfRays_NLOS; k++){
						elevation_ASA[i][0][j][k] = 90;
						elevation_ASD[i][0][j][k] = 90;
					}
				}
			}else{
				elevation_ASA[i][itIdx] = new double*[N_cluster_NLOS];
				elevation_ASD[i][itIdx] = new double*[N_cluster_NLOS];
				for(int j = 0; j < N_cluster_NLOS; j++){
					elevation_ASA[i][itIdx][j] = new double[numOfRays_NLOS];
					elevation_ASD[i][itIdx][j] = new double[numOfRays_NLOS];
					for(int k = 0; k < numOfRays_NLOS; k++){
						elevation_ASA[i][itIdx][j][k] = 90;
						elevation_ASD[i][itIdx][j][k] = 90;
					}
				}
				itIdx++;
			}
		}
	} 
/*
	double **elevation_cluster_ASA;
	elevation_cluster_ASA = new double*[numberOfMobileStations];
	elevation_ASA = new double***[numberOfMobileStations];
	
	// Generate elevation angles of arrival
	for(int i = 0; i < numberOfMobileStations; i++){
		elevation_ASA[i] = new double**[numOfInterferers];
		int itIdx = 1;
		for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
			if(it->first == bsId){
				if(LOSCondition[i][0]){
					elevation_ASA[i][0] = new double*[N_cluster_LOS];
					elevation_cluster_ASA[i] = new double[N_cluster_LOS];
			
					double X_1 = uniform(-1.0,1.0);
					double Y_1 = normal(0.0, (pow(sigma_zsA_LOS[i][0] / 7.0,2)));
			
					for(int j = 0; j < N_cluster_LOS; j++){
						elevation_ASA[i][0][j] = new double[numOfRays_LOS];
				
						//Formula 7.47
						elevation_cluster_ASA[i][j] = (sigma_zsA_LOS[i][0] / C_ZS(N_cluster_LOS, true)) * -1.0 * log(clusterPowers[i][0][j] / *(std::max_element(clusterPowers[i][0], clusterPowers[i][0] + N_cluster_LOS)));
						double X_N = uniform(-1.0,1.0);
						double Y_N = normal(0.0, (pow(sigma_zsA_LOS[i][0] / 7.0,2)));
				
						// First Ray has geometric LOS direction?!
						elevation_cluster_ASA[i][j] = elevation_cluster_ASA[i][j] * (X_N - X_1) + (Y_N - Y_1) + (ZoA_LOS_dir[i][0]*180.0/pi);  // since ZoA_LOS_dir is in radians
				
						for(int k = 0; k < numOfRays_LOS; k++){
							double cluster_ZSA = module->par("Cluster_ZSA_LOS");
					
							// Final angle computation per ray (LOS)
							elevation_ASA[i][0][j][k] = elevation_cluster_ASA[i][j] + cluster_ZSA * ray_offset[k];

							if(elevation_ASA[i][0][j][k] > 180.0){
								elevation_ASA[i][0][j][k] = 360.0 - elevation_ASA[i][0][j][k];
							}
						}
					}
				}else{
					elevation_ASA[i][0] = new double*[N_cluster_NLOS];
					elevation_cluster_ASA[i] = new double[N_cluster_NLOS];
			
					for(int j = 0; j < N_cluster_NLOS; j++){
						elevation_ASA[i][0][j] = new double[numOfRays_NLOS];
				
						// Formula 7.47
						elevation_cluster_ASA[i][j] = (sigma_zsA_NLOS[i][0] / C_ZS(N_cluster_NLOS, false)) * -1.0 * log(clusterPowers[i][0][j] / *(std::max_element(clusterPowers[i][0], clusterPowers[i][0] + N_cluster_NLOS)));
						double X_N = uniform(-1.0,1.0);
						double Y_N = normal(0.0, (pow(sigma_zsA_NLOS[i][0] / 7.0,2)));
				
						elevation_cluster_ASA[i][j] = elevation_cluster_ASA[i][j] * X_N + Y_N + (ZoA_LOS_dir[i][0]*180.0/pi);   // since ZoA_LOS_dir is in radians
				
						for(int k = 0; k < numOfRays_NLOS; k++){
							double cluster_ZSA = module->par("Cluster_ZSA_NLOS");
					
							// Final angle computation per ray (NLOS)
							elevation_ASA[i][0][j][k] = elevation_cluster_ASA[i][j] + cluster_ZSA * ray_offset[k];

							if(elevation_ASA[i][0][j][k] > 180.0){
								elevation_ASA[i][0][j][k] = 360.0 - elevation_ASA[i][0][j][k];
							}
						}
					}
				}
			}else{
				if(LOSCondition[i][itIdx]){
					elevation_ASA[i][itIdx] = new double*[N_cluster_LOS];
					elevation_cluster_ASA[i] = new double[N_cluster_LOS];
			
					double X_1 = uniform(-1.0,1.0);
					double Y_1 = normal(0.0, (pow(sigma_zsA_LOS[i][itIdx] / 7.0,2)));
			
					for(int j = 0; j < N_cluster_LOS; j++){
						elevation_ASA[i][itIdx][j] = new double[numOfRays_LOS];
				
						//Formula 7.47
						elevation_cluster_ASA[i][j] = (sigma_zsA_LOS[i][itIdx] / C_ZS(N_cluster_LOS, true)) * -1.0 * log(clusterPowers[i][itIdx][j] / *(std::max_element(clusterPowers[i][itIdx], clusterPowers[i][itIdx] + N_cluster_LOS)));
						double X_N = uniform(-1.0,1.0);
						double Y_N = normal(0.0, (pow(sigma_zsA_LOS[i][itIdx] / 7.0,2)));
				
						// First Ray has geometric LOS direction?!
						elevation_cluster_ASA[i][j] = elevation_cluster_ASA[i][j] * (X_N - X_1) + (Y_N - Y_1) + (ZoA_LOS_dir[i][itIdx]*180.0/pi);  // since ZoA_LOS_dir is in radians
				
						for(int k = 0; k < numOfRays_LOS; k++){
							double cluster_ZSA = module->par("Cluster_ZSA_LOS");
					
							// Final angle computation per ray (LOS)
							elevation_ASA[i][itIdx][j][k] = elevation_cluster_ASA[i][j] + cluster_ZSA * ray_offset[k];

							if(elevation_ASA[i][itIdx][j][k] > 180.0){
								elevation_ASA[i][itIdx][j][k] = 360.0 - elevation_ASA[i][itIdx][j][k];
							}
						}
					}
					
				}else{
					elevation_ASA[i][itIdx] = new double*[N_cluster_NLOS];
					elevation_cluster_ASA[i] = new double[N_cluster_NLOS];
			
					for(int j = 0; j < N_cluster_NLOS; j++){
						elevation_ASA[i][itIdx][j] = new double[numOfRays_NLOS];
				
						// Formula 7.47
						elevation_cluster_ASA[i][j] = (sigma_zsA_NLOS[i][itIdx] / C_ZS(N_cluster_NLOS, false)) * -1.0 * log(clusterPowers[i][itIdx][j] / *(std::max_element(clusterPowers[i][itIdx], clusterPowers[i][itIdx] + N_cluster_NLOS)));
						double X_N = uniform(-1.0,1.0);
						double Y_N = normal(0.0, (pow(sigma_zsA_NLOS[i][itIdx] / 7.0,2)));
				
						elevation_cluster_ASA[i][j] = elevation_cluster_ASA[i][j] * X_N + Y_N + (ZoA_LOS_dir[i][itIdx]*180.0/pi);   // since ZoA_LOS_dir is in radians
				
						for(int k = 0; k < numOfRays_NLOS; k++){
							double cluster_ZSA = module->par("Cluster_ZSA_NLOS");
					
							// Final angle computation per ray (NLOS)
							elevation_ASA[i][itIdx][j][k] = elevation_cluster_ASA[i][j] + cluster_ZSA * ray_offset[k];

							if(elevation_ASA[i][itIdx][j][k] > 180.0){
								elevation_ASA[i][itIdx][j][k] = 360.0 - elevation_ASA[i][itIdx][j][k];
							}
						}
					}
				}
				itIdx++;
			}
		}
	}
	
	double **elevation_cluster_ASD;
	elevation_cluster_ASD = new double*[numberOfMobileStations];
	elevation_ASD = new double***[numberOfMobileStations];
	
	// Generate azimuth angles of departure in the same way 
	for(int i = 0; i < numberOfMobileStations; i++){
		elevation_ASD[i] = new double**[numOfInterferers];
		int itIdx = 1;
		for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
			if(it->first == bsId){
				if(LOSCondition[i][0]){
					elevation_ASD[i][0] = new double*[N_cluster_LOS];
					elevation_cluster_ASD[i] = new double[N_cluster_LOS];
			
					double X_1 = uniform(-1.0,1.0);
					double Y_1 = normal(0.0, (pow(sigma_zsD_LOS[i][0] / 7.0,2)));
			
					for(int j = 0; j < N_cluster_LOS; j++){
						elevation_ASD[i][0][j] = new double[numOfRays_LOS];
				
						//Formula 7.47
						elevation_cluster_ASD[i][j] = (sigma_zsD_LOS[i][0] / C_ZS(N_cluster_LOS, true)) * -1.0 * log(clusterPowers[i][0][j] / *(std::max_element(clusterPowers[i][0], clusterPowers[i][0] + N_cluster_LOS)));
						double X_N = uniform(-1.0,1.0);
						double Y_N = normal(0.0, (pow(sigma_zsD_LOS[i][0] / 7.0,2)));
				
						// First Ray has geometric LOS direction?!
						elevation_cluster_ASD[i][j] = elevation_cluster_ASD[i][j] * (X_N - X_1) + (Y_N - Y_1) + (ZoD_LOS_dir[i][0]*180.0/pi);   // since ZoD_LOS_dir is in radians
				
						for(int k = 0; k < numOfRays_LOS; k++){
							dist2D = sqrt(pow((xPos - MSPos[i].x),2) + pow((yPos - MSPos[i].y),2));
							double mu_ZSD = mean_ZSD(dist2D, heightUE, true);
					
							// Final angle computation per ray (LOS)
							elevation_ASD[i][0][j][k] = elevation_cluster_ASD[i][j] + ((3.0/8.0) * ray_offset[k] * pow(10.0, mu_ZSD));

							if(elevation_ASD[i][0][j][k] > 180.0){
								elevation_ASD[i][0][j][k] = 360.0 - elevation_ASD[i][0][j][k];
							}
						}
					}
				}else{
					elevation_ASD[i][0] = new double*[N_cluster_NLOS];
					elevation_cluster_ASD[i] = new double[N_cluster_NLOS];
			
					for(int j = 0; j < N_cluster_NLOS; j++){
						elevation_ASD[i][0][j] = new double[numOfRays_NLOS];
				
						// Formula 7.47
						elevation_cluster_ASD[i][j] = (sigma_zsD_NLOS[i][0] / C_ZS(N_cluster_NLOS, false)) * -1.0 * log(clusterPowers[i][0][j] / *(std::max_element(clusterPowers[i][0], clusterPowers[i][0] + N_cluster_NLOS)));
						double X_N = uniform(-1.0,1.0);
						double Y_N = normal(0.0, (pow(sigma_zsD_NLOS[i][0] / 7.0,2)));
						dist2D = sqrt(pow((xPos - MSPos[i].x),2) + pow((yPos - MSPos[i].y),2));
						double muOffsetZoD = -1.0 * pow(10.0, (-1.0 * 0.62 * log10(std::max(10.0, dist2D))) + 2.035 - (0.07 * heightUE));
				
						elevation_cluster_ASD[i][j] = elevation_cluster_ASD[i][j] * X_N + Y_N + (ZoD_LOS_dir[i][0]*180.0/pi) + muOffsetZoD;  // since ZoD_LOS_dir is in radians
				
						for(int k = 0; k < numOfRays_NLOS; k++){
							double mu_ZSD = mean_ZSD(dist2D, heightUE, false);
					
							// Final angle computation per ray (NLOS)
							elevation_ASD[i][0][j][k] = elevation_cluster_ASD[i][j] + ((3.0/8.0) * ray_offset[k] * pow(10.0, mu_ZSD));

							if(elevation_ASD[i][0][j][k] > 180.0){
								elevation_ASD[i][0][j][k] = 360.0 - elevation_ASD[i][0][j][k];
							}
						}
					}
				}
			}else{
				if(LOSCondition[i][itIdx]){
					elevation_ASD[i][itIdx] = new double*[N_cluster_LOS];
					elevation_cluster_ASD[i] = new double[N_cluster_LOS];
			
					double X_1 = uniform(-1.0,1.0);
					double Y_1 = normal(0.0, (pow(sigma_zsD_LOS[i][itIdx] / 7.0,2)));
			
					for(int j = 0; j < N_cluster_LOS; j++){
						elevation_ASD[i][itIdx][j] = new double[numOfRays_LOS];
				
						//Formula 7.47
						elevation_cluster_ASD[i][j] = (sigma_zsD_LOS[i][itIdx] / C_ZS(N_cluster_LOS, true)) * -1.0 * log(clusterPowers[i][itIdx][j] / *(std::max_element(clusterPowers[i][itIdx], clusterPowers[i][itIdx] + N_cluster_LOS)));
						double X_N = uniform(-1.0,1.0);
						double Y_N = normal(0.0, (pow(sigma_zsD_LOS[i][itIdx] / 7.0,2)));
				
						// First Ray has geometric LOS direction?!
						elevation_cluster_ASD[i][j] = elevation_cluster_ASD[i][j] * (X_N - X_1) + (Y_N - Y_1) + (ZoD_LOS_dir[i][itIdx]*180.0/pi);   // since ZoD_LOS_dir is in radians
				
						for(int k = 0; k < numOfRays_LOS; k++){
							dist2D = sqrt(pow((it->second.x - MSPos[i].x),2) + pow((it->second.y - MSPos[i].y),2));
							double mu_ZSD = mean_ZSD(dist2D, heightUE, true);
					
							// Final angle computation per ray (LOS)
							elevation_ASD[i][itIdx][j][k] = elevation_cluster_ASD[i][j] + ((3.0/8.0) * ray_offset[k] * pow(10.0, mu_ZSD));

							if(elevation_ASD[i][itIdx][j][k] > 180.0){
								elevation_ASD[i][itIdx][j][k] = 360.0 - elevation_ASD[i][itIdx][j][k];
							}
						}
					}
					
				}else{
					elevation_ASD[i][itIdx] = new double*[N_cluster_NLOS];
					elevation_cluster_ASD[i] = new double[N_cluster_NLOS];
			
					for(int j = 0; j < N_cluster_NLOS; j++){
						elevation_ASD[i][itIdx][j] = new double[numOfRays_NLOS];
				
						// Formula 7.47
						elevation_cluster_ASD[i][j] = (sigma_zsD_NLOS[i][itIdx] / C_ZS(N_cluster_NLOS, false)) * -1.0 * log(clusterPowers[i][itIdx][j] / *(std::max_element(clusterPowers[i][itIdx], clusterPowers[i][itIdx] + N_cluster_NLOS)));
						double X_N = uniform(-1.0,1.0);
						double Y_N = normal(0.0, (pow(sigma_zsD_NLOS[i][itIdx] / 7.0,2)));
						dist2D = sqrt(pow((it->second.x - MSPos[i].x),2) + pow((it->second.y - MSPos[i].y),2));
						double muOffsetZoD = -1.0 * pow(10.0, (-1.0 * 0.62 * log10(std::max(10.0, dist2D))) + 2.035 - (0.07 * heightUE));
				
						elevation_cluster_ASD[i][j] = elevation_cluster_ASD[i][j] * X_N + Y_N + (ZoD_LOS_dir[i][itIdx]*180.0/pi) + muOffsetZoD;  // since ZoD_LOS_dir is in radians
				
						for(int k = 0; k < numOfRays_NLOS; k++){
							double mu_ZSD = mean_ZSD(dist2D, heightUE, false);
					
							// Final angle computation per ray (NLOS)
							elevation_ASD[i][itIdx][j][k] = elevation_cluster_ASD[i][j] + ((3.0/8.0) * ray_offset[k] * pow(10.0, mu_ZSD));

							if(elevation_ASD[i][itIdx][j][k] > 180.0){
								elevation_ASD[i][itIdx][j][k] = 360.0 - elevation_ASD[i][itIdx][j][k];
							}
						}
					}
				}
				itIdx++;
			}
		}
	}
*/
	// Generate random phases (7.3.17)
	randomPhase = new double****[numberOfMobileStations];
	randomPhase_LOS = new double*[numberOfMobileStations];					// Random phase for LOS component
	for(int i = 0; i < numberOfMobileStations; i++){
		int itIdx = 1;
		randomPhase[i] = new double***[neighbourPositions.size()];
		randomPhase_LOS[i] = new double[neighbourPositions.size()];
		for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
			if(it->first == bsId){
				if(LOSCondition[i][0]){
					randomPhase[i][0] = new double**[N_cluster_LOS];
					randomPhase_LOS[i][0] = uniform(-1.0*pi, pi);
			
					for(int j = 0; j < N_cluster_LOS; j++){
						randomPhase[i][0][j] = new double*[numOfRays_LOS];
				
						for(int k = 0; k < numOfRays_LOS; k++){
							randomPhase[i][0][j][k] = new double[4];
							
							for(int l = 0; l < 4; l++){
								randomPhase[i][0][j][k][l] = uniform(-1.0*pi, pi);				// for the random phases of NLOS component in equation 7-61 of METIS 1.2	
							}
						}
					}
				} else {
					randomPhase[i][0] = new double**[N_cluster_NLOS];
			
					for(int j = 0; j < N_cluster_NLOS; j++){
						randomPhase[i][0][j] = new double*[numOfRays_NLOS];
				
						for(int k = 0; k < numOfRays_NLOS; k++){
							randomPhase[i][0][j][k] = new double[4];
							
							for(int l = 0; l < 4; l++){
								randomPhase[i][0][j][k][l] = uniform(-1.0*pi, pi);				// for the random phases as in equation 7-59 of METIS 1.2 (one each for ThetaTheta, ThetaPhi, PhiTheta and PhiPhi)
							}
						}
					}
				}
			} else {
				if(LOSCondition[i][itIdx]){
					randomPhase[i][itIdx] = new double**[N_cluster_LOS];
					randomPhase_LOS[i][itIdx] = uniform(-1.0*pi, pi);
			
					for(int j = 0; j < N_cluster_LOS; j++){
						randomPhase[i][itIdx][j] = new double*[numOfRays_LOS];
				
						for(int k = 0; k < numOfRays_LOS; k++){
							randomPhase[i][itIdx][j][k] = new double[4];
							
							for(int l = 0; l < 4; l++){
								randomPhase[i][itIdx][j][k][l] = uniform(-1.0*pi, pi);			// for the random phases of NLOS component in equation 7-61 of METIS 1.2
							}
						}
					}
				}else{
					randomPhase[i][itIdx] = new double**[N_cluster_NLOS];
			
					for(int j = 0; j < N_cluster_NLOS; j++){
						randomPhase[i][itIdx][j] = new double*[numOfRays_NLOS];
				
						for(int k = 0; k < numOfRays_NLOS; k++){
							randomPhase[i][itIdx][j][k] = new double[4];
							
							for(int l = 0; l < 4; l++){
								randomPhase[i][itIdx][j][k][l] = uniform(-1.0*pi, pi);			// for the random phases as in equation 7-59 of METIS 1.2 (one each for ThetaTheta, ThetaPhi, PhiTheta and PhiPhi)
							}
						}
					}
				}
				itIdx++;
			}
		}
	}
	
	// Generate cross polarization values
	double XPR_Mean_LOS = module->par("XPR_Mean_LOS");
	double XPR_Std_LOS = module->par("XPR_Std_LOS");
	double XPR_Mean_NLOS = module->par("XPR_Mean_NLOS");
	double XPR_Std_NLOS = module->par("XPR_Std_NLOS");
	Xn_m = new double***[numberOfMobileStations];
	
	for(int i = 0; i < numberOfMobileStations; i++){
		Xn_m[i] = new double**[neighbourPositions.size()];
		int id_BS =1;
		for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
			if(it->first == bsId){
				if(LOSCondition[i][0]){
					Xn_m[i][0] = new double*[N_cluster_LOS];
					for(int j = 0; j < N_cluster_LOS; j++){
						Xn_m[i][0][j] = new double[numOfRays_LOS];
				
						for(int k = 0; k < numOfRays_LOS; k++){
							Xn_m[i][0][j][k] = normal(XPR_Mean_LOS, pow(XPR_Std_LOS,2) );
						}
					}
				} else {
					Xn_m[i][0] = new double*[N_cluster_NLOS];
					for(int j = 0; j < N_cluster_NLOS; j++){
						Xn_m[i][0][j] = new double[numOfRays_NLOS];
				
						for(int k = 0; k < numOfRays_NLOS; k++){
							Xn_m[i][0][j][k] = normal(XPR_Mean_NLOS, pow(XPR_Std_NLOS,2) );
						}
					}
				}
			} else {
				if(LOSCondition[i][id_BS]){
					Xn_m[i][id_BS] = new double*[N_cluster_LOS];
					for(int j = 0; j < N_cluster_LOS; j++){
						Xn_m[i][id_BS][j] = new double[numOfRays_LOS];
				
						for(int k = 0; k < numOfRays_LOS; k++){
							Xn_m[i][id_BS][j][k] = normal(XPR_Mean_LOS, pow(XPR_Std_LOS,2) );
						}
					}
				} else {
					Xn_m[i][id_BS] = new double*[N_cluster_NLOS];
					for(int j = 0; j < N_cluster_NLOS; j++){
						Xn_m[i][id_BS][j] = new double[numOfRays_NLOS];
				
						for(int k = 0; k < numOfRays_NLOS; k++){
							Xn_m[i][id_BS][j][k] = normal(XPR_Mean_NLOS, pow(XPR_Std_NLOS,2) );
						}
					}
				}
				id_BS++;
			}
		}
	}
	
	//output2 << "Init 2 METIS at BS " << bsId << " with rand: " << normal(0,1) << std::endl;
	
   	//---------------------------------------------------------------
   	/*
   	// Own Implementation:
	// The first and second exp expression in final H Computation are
	// equivalent to what is called "MIMO Matrix" within Cost2100. This
	// is identical, even if SISO is used, however it is not a matrix
	// then, but a scalar (yet complex) value. Lets call it AoA_AoD_Comp.
	// Formula: AoA_AoD_Comp = exp(j * k_0 * dot(AoA,r_rx)) * exp(j * k_0 * dot(AoD,r_tx))
	// where:
	// k_0 = wavenumber = (2 * pi) / lambda (lambda = speedOfLight / carrierFreq)
	// r_rx = position vector (x,y,z)_rx of receiver antenna (Matrix for MIMO)
	// R-tx = position vector (x,y,z)_tx of transmitter antenna (Matrix for MIMO)
	// AoA = spherical unit vector of Angle of Arrival
	// AoD = spherical unit vector of Angle of Departure

	// First Compute constant values once
	//double k_0 = (2 * pi) / (speedOfLight / freq_c);
	
	// TODO: Add for Loops/indices.
	double AoA[3];
	AoA[0] = sin(azimuth_ASA[][][]) * cos(elevation_ASA[][][]);
	AoA[1] = sin(azimuth_ASA[][][]) * sin(elevation_ASA[][][]);
	AoA[2] = cos(azimuth_ASA[][][]);
	
	double AoD[3];
	AoD[0] = sin(azimuth_ASD[][][]) * cos(elevation_ASD[][][]);
	AoD[1] = sin(azimuth_ASD[][][]) * sin(elevation_ASD[][][]);
	AoD[2] = cos(azimuth_ASD[][][]);
	
	complex<double> AoA_AoD_Comp = 	exp( complex(0, k_0 * (MSPos.x * AoA[0] + MSPos.y * AoA[1] + heightUE * AoA[2])) ) * 
									exp( complex(0, k_0 * (xPos * AoD[0] + yPos * AoD[1] + heightBS * AoD[2])) );
									
	// The Last exp expression is the Doppler component.
	// Formula: doppler = exp( j * k_0 * dot(AoA,V_rx) * t)
	// Where:
	// k_0 = wavenumber = (2 * pi) / lambda (lambda = speedOfLight / carrierFreq)
	// AoA = spherical unit vector of Angle of Arrival
	// V_rx = Velocity vector of receiver (MS) in global coordinate system
	// t = linear time vector (0.0,...,endTime) in timeSamples steps
	
	complex<double> doppler = exp( complex(0, k_0 * (MSVelDir[][0] * AoA[0] + MSVelDir[][1] * AoA[1] + MSVelDir[][2] * AoA[2])) * timeVector[][] );
	
	// The depolarization Matrix in different in WINNER II and METIS Model.
	// Not only in LOS case, but also in NLOS case.
	// Formula: 
	// [ exp( j * uniform(0,2*pi) )					...		exp( j * uniform(0,2*pi) ) / sqrt(Xn_m) ]
	// [ exp( j * uniform(0,2*pi) ) / sqrt(Xn_m)	...		exp( j * uniform(0,2*pi) ) 				]
	// where:
	// Xn_m = Cross Polarization ratio
	double depolMatrix[2][2];
	depolMatrix[0][0] = exp( complex<double>(0,uniform(0,2*pi)) );
	depolMatrix[1][1] = exp( complex<double>(0,uniform(0,2*pi)) );

	depolMatrix[1][0] = exp( complex<double>(0,uniform(0,2*pi)) / );
	depolMatrix[0][1] = exp( complex<double>(0,uniform(0,2*pi)) / );
	
	// The last part is the MSGain and BSGain. This is to be multiplied with
	// the depolarization Matrix once from the left and once from the right,
	// resulting in a scalar (complex) value.
	
	// Finally we add the cluster power.
	// Formula:
	// H_n = sqrt(P_n / M) * sumRays_n
	//double sq_P_over_M = sqrt(P_n / numOfRays_NLOS);
	*/
	//------------------------------------------------------------------------
	
	
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
	complex<double> *****raySum_LOS;
	int clusterIdx_LOS;
	
	raySum_LOS = new complex<double>****[numberOfMobileStations];
	
	for(int m = 0; m < numberOfMobileStations; m++){
		// +4 Clusters, because 2 are subdivided into 6, resulting in 4 more.
		raySum_LOS[m] = new complex<double>***[(N_cluster_LOS + 4)];
		for(int i = 0; i < (N_cluster_LOS + 4); i++){
			raySum_LOS[m][i] = new complex<double>**[timeSamples];
			for(int j = 0; j < timeSamples; j++){
				raySum_LOS[m][i][j] = new complex<double>*[NumRxAntenna];
				for(int k = 0; k < NumRxAntenna; k++){
					raySum_LOS[m][i][j][k] = new complex<double>[NumTxAntenna]();
				}
			}
		}
	}
	
	// For NLOS case:
	complex<double> *****raySum;
	int clusterIdx;
	raySum = new complex<double>****[numberOfMobileStations];

	for(int m = 0; m < numberOfMobileStations; m++){
		// +4 Clusters, because 2 are subdivided into 6, resulting in 4 more.
		raySum[m] = new complex<double>***[(N_cluster_NLOS + 4)];
			for(int i = 0; i < (N_cluster_NLOS + 4); i++){
			raySum[m][i] = new complex<double>**[timeSamples];
			for(int j = 0; j < timeSamples; j++){
				raySum[m][i][j] = new complex<double>*[NumRxAntenna];
				for(int k = 0; k < NumRxAntenna; k++){
					raySum[m][i][j][k] = new complex<double>[NumTxAntenna]();
				}
			}
		}
	}
	
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
						raySumInterferer_LOS[m][i][j][k] = new complex<double>*[NumRxAntenna];
						for(int l = 0; l < NumRxAntenna; l++){
							raySumInterferer_LOS[m][i][j][k][l] = new complex<double>[NumTxAntenna]();
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
						raySumInterferer[m][i][j][k] = new complex<double>*[NumRxAntenna];
						for(int l = 0; l < NumRxAntenna; l++){
							raySumInterferer[m][i][j][k][l] = new complex<double>[NumTxAntenna]();
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
					for(int u = 0; u < NumRxAntenna; u++){
						// Cycle through all Transmitter antennas (BS)
						for(int s = 0; s < NumTxAntenna; s++){
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
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2] )) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
														//std::cout << "raySum[i][clusterIdx][t][u][s]: " << raySum[i][clusterIdx][t][u][s] << std::endl;
													}
												} // End time axis							
											} // End cycle Rays
											clusterIdx++;
										} else if (l == 1){ // for sub-cluster 2 of first two clusters
											sq_P_over_M = sqrt(P_n[l] / size_SC_2);
											//std::cout << "P_n[l]: " << P_n[l] << std::endl;
											//std::cout << "sq_P_over_M: " << sq_P_over_M << std::endl;
											// Cycle through all NLOS Rays
											for(int m = 0; m < size_SC_2; m++){		
												AoA[0] = sin(elevation_ASA[i][0][n][SC_2[m]]*pi/180) * cos(azimuth_ASA[i][0][n][SC_2[m]]*pi/180);
												AoA[1] = sin(elevation_ASA[i][0][n][SC_2[m]]*pi/180) * sin(azimuth_ASA[i][0][n][SC_2[m]]*pi/180);
												AoA[2] = cos(elevation_ASA[i][0][n][SC_2[m]]*pi/180);
										
												AoD[0] = sin(elevation_ASD[i][0][n][SC_2[m]]*pi/180) * cos(azimuth_ASD[i][0][n][SC_2[m]]*pi/180);
												AoD[1] = sin(elevation_ASD[i][0][n][SC_2[m]]*pi/180) * sin(azimuth_ASD[i][0][n][SC_2[m]]*pi/180);
												AoD[2] = cos(elevation_ASD[i][0][n][SC_2[m]]*pi/180);
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
														//std::cout << "raySum[i][clusterIdx][t][u][s]: " << raySum[i][clusterIdx][t][u][s] << std::endl;
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
												AoA[0] = sin(elevation_ASA[i][0][n][SC_3[m]]*pi/180) * cos(azimuth_ASA[i][0][n][SC_3[m]]*pi/180);
												AoA[1] = sin(elevation_ASA[i][0][n][SC_3[m]]*pi/180) * sin(azimuth_ASA[i][0][n][SC_3[m]]*pi/180);
												AoA[2] = cos(elevation_ASA[i][0][n][SC_3[m]]*pi/180);
										
												AoD[0] = sin(elevation_ASD[i][0][n][SC_3[m]]*pi/180) * cos(azimuth_ASD[i][0][n][SC_3[m]]*pi/180);
												AoD[1] = sin(elevation_ASD[i][0][n][SC_3[m]]*pi/180) * sin(azimuth_ASD[i][0][n][SC_3[m]]*pi/180);
												AoD[2] = cos(elevation_ASD[i][0][n][SC_3[m]]*pi/180);
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
														//std::cout << "raySum[i][clusterIdx][t][u][s]: " << raySum[i][clusterIdx][t][u][s] << std::endl;
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
										
											complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2])) );
											complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
													//std::cout << "raySum[i][clusterIdx][t][u][s]: " << raySum[i][clusterIdx][t][u][s] << std::endl;
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
					for(int u = 0; u < NumRxAntenna; u++){
						// Cycle through all Transmitter antennas (BS)
						for(int s = 0; s < NumTxAntenna; s++){
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
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
														//std::cout << "raySum_LOS[i][clusterIdx][t][u][s]: " << raySum_LOS[i][clusterIdx][t][u][s] << std::endl;
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
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
														//std::cout << "raySum_LOS[i][clusterIdx_LOS][t][u][s]: " << raySum_LOS[i][clusterIdx_LOS][t][u][s] << std::endl;
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
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
														//std::cout << "raySum_LOS[i][clusterIdx_LOS][t][u][s]: " << raySum_LOS[i][clusterIdx_LOS][t][u][s] << std::endl;
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
										
											complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2])) );
											complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
													//std::cout << "raySum_LOS[i][clusterIdx_LOS][t][u][s]: " << raySum_LOS[i][clusterIdx_LOS][t][u][s] << std::endl;
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
										
									complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2])) );
									complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
					for(int u = 0; u < NumRxAntenna; u++){
						// Cycle through all Transmitter antennas (BS)
						for(int s = 0; s < NumTxAntenna; s++){
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
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
										
											complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2])) );
											complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
					for(int u = 0; u < NumRxAntenna; u++){
						// Cycle through all Transmitter antennas (BS)
						for(int s = 0; s < NumTxAntenna; s++){
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
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
										
												complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2])) );
												complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
										
											complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2])) );
											complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
										
									complex<double> exp_arrival = exp( complex<double>(0.0,k_0 * (AoA[0] * RxAntennaPosition[i][u][0] + AoA[1] * RxAntennaPosition[i][u][1] + AoA[2] * RxAntennaPosition[i][u][2])) );
									complex<double> exp_departure = exp( complex<double>(0.0,k_0 * (AoD[0] * TxAntennaPosition[0][s][0] + AoD[1] * TxAntennaPosition[0][s][1] + AoD[2] * TxAntennaPosition[0][s][2])) );
										
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
		TimeFrequency.open ("/home/sahari/GSCM_Channel_model/abstractLTEChannelModel/results/TimeFrequency_BS_" + std::to_string( (long long) bsId ) + "_" + std::to_string((MSPos[i].x * 100)) + ".txt");
		dist2D = sqrt(pow((xPos - MSPos[i].x),2) + pow((yPos - MSPos[i].y),2));
		dist3D = sqrt(pow(dist2D,2) + pow((heightBS - heightUE),2));
		complex<double> res = complex<double>(0.0,0.0);
		for(int t = 0; t < timeSamples; t++){
			for(int f = 0; f < downRBs; f++){
				double freq_ = freq_c + f*180000;
				res = complex<double>(0.0,0.0);
				if(LOSCondition[i][0]){
					pathloss = CalcPathloss(dist2D, dist3D, true);
					for(int u = 0; u < NumRxAntenna; u++){
						for(int s = 0; s < NumTxAntenna; s++){
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
					TimeFrequency << pathloss * tempRes << " ";
				}else{
					pathloss = CalcPathloss(dist2D, dist3D, false);
					for(int u = 0; u < NumRxAntenna; u++){
						for(int s = 0; s < NumTxAntenna; s++){
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
					TimeFrequency << pathloss * tempRes << " ";
				} //End if condition for LOS/NLOS
			} //End loop for RBs
			TimeFrequency << "\n\n";
		}
		TimeFrequency << pathloss << "\n";
		TimeFrequency.close();
	}
	
	if(numOfInterferers>0){
		//Only compute SINRneighbour values if there actually are any 
		// neighbours to influence transmissions
		SINRneighbour = new double***[numberOfMobileStations];
		for(int i = 0; i < numberOfMobileStations; i++){
			SINRneighbour[i] = new double**[numOfInterferers];
			for(int j = 0; j < numOfInterferers; j++){
				SINRneighbour[i][j] = new double*[timeSamples];
				for(int k = 0; k < timeSamples; k++){
					SINRneighbour[i][j][k] = new double[downRBs];
				}
			}
		}
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
								for(int u = 0; u < NumRxAntenna; u++){
									for(int s = 0; s < NumTxAntenna; s++){
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
								for(int u = 0; u < NumRxAntenna; u++){
									for(int s = 0; s < NumTxAntenna; s++){
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
	
	//output2 << "Init 3 METIS at BS " << bsId << " with rand: " << normal(0,1) << std::endl;

	
	// To write the channel coefficients to respective files
	
	std::cout << "START WRITING CHANNEL COEFFICIENTS for BS: " << bsId << std::endl;
	
	ofstream H_Time;
	//ofstream H_TimeInterferer;
	
	for(int i = 0; i < numberOfMobileStations; i++){
		H_Time.open ("/home/sahari/GSCM_Channel_model/abstractLTEChannelModel/results/H_Time_BS_" + std::to_string( (long long) bsId )  + "_" + std::to_string((MSPos[i].x * 100)) + ".txt");;
		complex<double> res = complex<double>(0.0,0.0);
		for(int t = 0; t < timeSamples; t++){
			for(int u = 0; u < NumRxAntenna; u++){
				for(int s = 0; s < NumTxAntenna; s++){
					res = complex<double>(0.0,0.0);
					if(LOSCondition[i][0]){
						for(int n = 0; n < N_cluster_LOS; n++){
							if(n < 2){
								res = res + raySum_LOS[i][n][t][u][s];
								res = res + raySum_LOS[i][n+1][t][u][s];
								res = res + raySum_LOS[i][n+2][t][u][s];
								res = res + raySum_LOS[i][n+3][t][u][s];
								res = res + raySum_LOS[i][n+4][t][u][s];
								res = res + raySum_LOS[i][n+5][t][u][s];
								n++;
							}else{
								res = res + raySum_LOS[i][n + 4][t][u][s];
							}
						}
					}else{
						for(int n = 0; n < N_cluster_NLOS; n++){
							if(n < 2){ 
								res = res + raySum[i][n][t][u][s];
								res = res + raySum[i][n+1][t][u][s];
								res = res + raySum[i][n+2][t][u][s];
								res = res + raySum[i][n+3][t][u][s];
								res = res + raySum[i][n+4][t][u][s];
								res = res + raySum[i][n+5][t][u][s];
								n++;
							}else{
								res = res + raySum[i][n + 4][t][u][s];
							}
						}
					}
					H_Time << res << " ";
				} //end loop TxAntenna
				H_Time << "\n";
			} //end loop RxAntenna
			//H_Time << "\n";
		}
		H_Time << "\n" << i+1 << "\n";
		H_Time.close();
	}
	
	for(int i = 0; i < numberOfMobileStations; i++){
		int idIdx = 0;
		//Channel Coefficients for interferers:
		for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
			if(it->first == bsId){
				// Skip own BS
				continue;
			}else{
				//H_TimeInterferer.open ("/home/sahari/GSCM_Channel_model/abstractLTEChannelModel/results/H_Time_Interferer_BS_" + std::to_string( (long long) bsId )  + "_" + std::to_string((MSPos[i].x * 100)) + ".txt");
				complex<double> res = complex<double>(0.0,0.0);
				for(int t = 0; t < timeSamples; t++){
					for(int u = 0; u < NumRxAntenna; u++){
						for(int s = 0; s < NumTxAntenna; s++){
							res = complex<double>(0.0,0.0);
							if(LOSCondition[i][idIdx+1]){
								for(int n = 0; n < N_cluster_LOS; n++){
									if(n < 2){
										res = res + raySumInterferer_LOS[i][idIdx][n][t][u][s];
										res = res + raySumInterferer_LOS[i][idIdx][n+1][t][u][s];
										res = res + raySumInterferer_LOS[i][idIdx][n+2][t][u][s];
										res = res + raySumInterferer_LOS[i][idIdx][n+3][t][u][s];
										res = res + raySumInterferer_LOS[i][idIdx][n+4][t][u][s];
										res = res + raySumInterferer_LOS[i][idIdx][n+5][t][u][s];
										n++;
									}else{
										res = res + raySumInterferer_LOS[i][idIdx][n + 4][t][u][s];
									}
								}
							}else{
								for(int n = 0; n < N_cluster_NLOS; n++){
									if(n < 2){ 
										res = res + raySumInterferer[i][idIdx][n][t][u][s];
										res = res + raySumInterferer[i][idIdx][n+1][t][u][s];
										res = res + raySumInterferer[i][idIdx][n+2][t][u][s];
										res = res + raySumInterferer[i][idIdx][n+3][t][u][s];
										res = res + raySumInterferer[i][idIdx][n+4][t][u][s];
										res = res + raySumInterferer[i][idIdx][n+5][t][u][s];
										n++;
									}else{
										res = res + raySumInterferer[i][idIdx][n + 4][t][u][s];
									}
								}
							}
							//H_TimeInterferer << res << " ";
						} //end loop TxAntenna
						//H_TimeInterferer << "\n";
					} // end loop RxAntenna
					//H_TimeInterferer << "\n\n";
				}
				//H_TimeInterferer << "\n" << i+1 << " " << idIdx+1 << "\n";
				idIdx++;
			}
		}
		//H_TimeInterferer.close();
	}
	
	std::cout << "FINISHED WRITING CHANNEL COEFFICIENTS for BS: " << bsId << std::endl;
	
	//------------------------------------------------------------------
	
	//output2 << "Init 3 METIS at BS " << bsId << " with rand: " << normal(0,1) << std::endl;

	// Free memory allocated on heap for local variables
	//
	delete neighbourIdMatching; 

	for(int m = 0; m < numberOfMobileStations; m++){
		for(int i = 0; i < (N_cluster_LOS + 4); i++){
			for(int j = 0; j < timeSamples; j++){
				for(int k = 0; k < NumRxAntenna; k++){
					delete[] raySum_LOS[m][i][j][k];
				}
				delete[] raySum_LOS[m][i][j];
			}
			delete[] raySum_LOS[m][i];
		}
		delete[] raySum_LOS[m];
		for(int i = 0; i < (N_cluster_NLOS + 4); i++){
			for(int j = 0; j < timeSamples; j++){
				for(int k = 0; k < NumRxAntenna; k++){
					delete[] raySum[m][i][j][k];
				}
				delete[] raySum[m][i][j];
			}
			delete[] raySum[m][i];
		}
		delete[] raySum[m];
	}
	delete[] raySum_LOS;
	delete[] raySum;

	if(numOfInterferers>0){
		for (int m = 0; m < numberOfMobileStations; m++){
			for(int i = 0; i < numOfInterferers; i++){
				for(int j = 0; j < (N_cluster_LOS + 4); j++){
					for(int k = 0; k < timeSamples; k++){
						for(int l = 0; l < NumRxAntenna; l++){
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

	for(int i = 0; i < numberOfMobileStations; i++){
		for(int s=0; s<neighbourPositions.size(); s++){
			for(int j=0;j<N_cluster_NLOS;j++){
				if(!LOSCondition[i][s]){
					for(int k = 0; k < numOfRays_NLOS; k++){
						delete[] randomPhase[i][0][j][k];
					}
					delete[] randomPhase[i][s][j];
				}
				delete[] elevation_ASA[i][s][j];
				delete[] elevation_ASD[i][s][j];
			}
			delete[] elevation_ASA[i][s];
			delete[] elevation_ASD[i][s];
			delete[] clusterDelays[i][s];
			if(LOSCondition[i][s]){
				for(int j=0;j<N_cluster_LOS;j++){
					for(int k = 0; k < numOfRays_LOS; k++){
						delete[] randomPhase[i][0][j][k];
					}
					delete[] randomPhase[i][s][j];
				}
				delete[] clusterDelays_LOS[i][s];
			}
			delete[] randomPhase[i][s];
			delete[] clusterPowers[i][s];
		}
		delete[] azimuth_cluster_ASA[i];
		delete[] azimuth_cluster_ASD[i];
		delete[] randomPhase[i];
		delete[] elevation_ASA[i];
		delete[] elevation_ASD[i];
		delete[] clusterDelays[i];
		delete[] clusterDelays_LOS[i];
	 	delete[] clusterPowers[i];
	}
	delete[] azimuth_cluster_ASA;
	delete[] azimuth_cluster_ASD;
	delete[] randomPhase;
	delete[] elevation_ASA;
	delete[] elevation_ASD;
	delete[] clusterDelays;
	delete[] clusterDelays_LOS;
	delete[] clusterPowers;

	return true;
}

/**
* Function that computes the scaling factor according to table 7.5 and formula 7.48
* @param numCluster The number of clusters per link
* @param LOS true iff LOS Condition, false otherwise
* @return Scaling factor
*/
double METISChannel::C_AS(int numCluster, bool LOS, int i){
	if(LOS){
		double K = 10.0 * log10(abs(sigma_kf_LOS[i][0]));
		double k_factor_LOS = 1.1035 - 0.028 * K - 0.002 * pow(K,2) + 0.0001 * pow(K,3);
		
		return C_AS(numCluster, false, i) * k_factor_LOS;
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
double METISChannel::C_ZS(int numCluster, bool LOS){
	int i = 0;
	if(LOS){
		double K = 10.0 * log10(abs(sigma_kf_LOS[i][0]));
		double k_factor_LOS = 1.3086 - 0.0339 * K - 0.0077 * pow(K,2) + 0.0002 * pow(K,3);
		
		return C_ZS(numCluster, false) * k_factor_LOS;
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
void METISChannel::generateAutoCorrelation_LOS(){
	int r = initModule->par("cellRadiusMETIS");
	std::cout << "Radius: " << r << std::endl;
	autoCorrelation_LOS = new double*[7];
	for(int i = 0; i < 7; i++){
		autoCorrelation_LOS[i] = new double[numberOfMobileStations];
		for(int j = 0; j < numberOfMobileStations; j++){
			autoCorrelation_LOS[i][j] = 0;
		}
	}
	
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
	double **filter; 
	double **sum; 
	
	grid = new double**[7];
	tmpX = new double**[7];
	tmpY = new double**[7];
	filter = new double*[7];
	sum = new double*[7];
	
	for(int i = 0; i < 7; i++){
		filter[i] = new double[11];
		sum[i] = new double[11];
	}
	
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
	for(int a = 0; a < numberOfMobileStations; a++){
		for(int i = 0; i < 7; i++){
			autoCorrelation_LOS[i][a] = tmpY[i][(int) (MSPos[a].x + middle - xPos)][(int) (MSPos[a].y + middle - yPos)];
			//autoCorrelation_LOS[i][a] = normal(0,1);
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
	for(int i=0; i<7; i++){
		delete[] filter[i];
		delete[] sum[i];
	}
	delete[] filter;
	delete[] sum;
}

/**
* Function that computes the exponential auto correlation for NLOS.
* @return True if succesful. False otherwise.
*/
void METISChannel::generateAutoCorrelation_NLOS(){
	int r = initModule->par("cellRadiusMETIS");
	std::cout << "Radius: " << r << std::endl;
	autoCorrelation_NLOS = new double*[6];
	for(int i = 0; i < 6; i++){
		autoCorrelation_NLOS[i] = new double[numberOfMobileStations];
		for(int j = 0; j < numberOfMobileStations; j++){
			autoCorrelation_NLOS[i][j] = 0;
		}
	}
	
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
	for(int a = 0; a < numberOfMobileStations; a++){
		for(int i = 0; i < 6; i++){
			autoCorrelation_NLOS[i][a] = tmpY[i][(int) (MSPos[a].x + middle - xPos)][(int) (MSPos[a].y + middle - yPos)];
			//autoCorrelation_NLOS[i][a] = normal(0,1);
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
	//Placeholder
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
	delete sigma_ds_LOS;					
	delete sigma_asD_LOS;					
	delete sigma_asA_LOS;
	delete sigma_zsD_LOS;
	delete sigma_zsA_LOS;					
	delete sigma_sf_LOS;					
	delete sigma_kf_LOS;					
	delete sigma_ds_NLOS;					
	delete sigma_asD_NLOS;					
	delete sigma_asA_NLOS;
	delete sigma_zsD_NLOS;
	delete sigma_zsA_NLOS;					
	delete sigma_sf_NLOS;						
	delete bs_antenna_bearing;				
	delete bs_antenna_downtilt;			
	delete bs_antenna_slant;				
	delete ms_antenna_bearing;				
	delete ms_antenna_downtilt;			
	delete ms_antenna_slant;				
}
