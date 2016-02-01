/*
 * ChannelAlternative.cc
 *
 *  Created on: Jul 22, 2014
 *      Author: Thomas Prinz
 * 
 */
 
#include "ChannelAlternative.h"
#include <stdio.h>

using namespace std;
using namespace itpp;

// Initialize your Channel through ini access via module pointer.
bool ChannelAlternative::init(cSimpleModule* module, Position** msPositions, std::map <int,Position> neighbourPositions){
	int nrCells = module->getParentModule()->getParentModule()->getParentModule()->par("numberOfCells");
	tenlogk = module->par("tenlogk");
	alpha = module->par("alpha");
	FADING_PATHS = module->par("FADING_PATHS");
	nrMSs = module->par("numberOfMobileStations");
	stationaryDoppler = module->par("stationaryDoppler");
	DELAY_RMS = module->par("delay_rms");
	RANDOM_SEED_ID = module->par("CHNrandom_seed_id");
	myBsID = module->par("bsId");
	int maxNumberOfNeighbours = module->par("maxNumberOfNeighbours");
	myID = module->getIndex();
	
	cModule *cell = module->getParentModule()->getParentModule();
    neighbourIdMatching = new NeighbourIdMatching(myBsID, maxNumberOfNeighbours, cell);
	
	for(int i = 0; i < nrMSs; i++){
		velocity.push_back(0);
	}
	
	// Initialize the arrays.
	bsXposition = new double [nrCells];
	bsYposition = new double [nrCells];
	msXposition = new double [nrMSs];
	msYposition = new double [nrMSs];
	shadowLoss  = new double [nrCells];
	
	int fromBsId = neighbourIdMatching->getDataStrId(myBsID);
	for(int i = 0; i < nrMSs; i++){
		msXposition[i] = msPositions[fromBsId][i].x;
		msYposition[i] = msPositions[fromBsId][i].y;
	}
	
	for(int i = 0; i < nrCells; i++){
		bsXposition[i] = neighbourPositions[i].x;
		bsYposition[i] = neighbourPositions[i].y;
		cout << "Cell " << i << " has Position: " << bsXposition[i] << " " << bsYposition[i] << endl;
	}
	
	if(myBsID == 3){
		ofstream output;
		string out = "pathloss.txt";
		output.open(out);
		for(int i = 0; i < nrMSs; i++){
			double dist = sqrt( pow(bsXposition[3] - msXposition[i],2) + pow(bsYposition[3] - msYposition[i],2) );
			double dist1 = sqrt( pow(bsXposition[0] - msXposition[i],2) + pow(bsYposition[0] - msYposition[i],2) );
			double dist2 = sqrt( pow(bsXposition[1] - msXposition[i],2) + pow(bsYposition[1] - msYposition[i],2) );
			double dist3 = sqrt( pow(bsXposition[2] - msXposition[i],2) + pow(bsYposition[2] - msYposition[i],2) );
			double dist4 = sqrt( pow(bsXposition[4] - msXposition[i],2) + pow(bsYposition[4] - msYposition[i],2) );
			double dist5 = sqrt( pow(bsXposition[5] - msXposition[i],2) + pow(bsYposition[5] - msYposition[i],2) );
			double dist6 = sqrt( pow(bsXposition[6] - msXposition[i],2) + pow(bsYposition[6] - msYposition[i],2) );
			//output << dist << " " << dist1 << " " << dist2 << " " << dist3 << " " << dist4 << " " << dist5 << " " << dist6 << " " << endl;
			output << calcPathloss(dist) << " " << calcPathloss(dist1) << " " << calcPathloss(dist2) << " " << calcPathloss(dist3) << " " << calcPathloss(dist4) << " " << calcPathloss(dist5) << " " << calcPathloss(dist6) << " " << endl;
		}
		output.close();
	}

	//double u, x;

	// Init Arrrays
	angle_of_arrival = new double **[nrMSs];
	delay = new double **[nrMSs];

	for (int pos_cell_ID = 0; pos_cell_ID < nrMSs; pos_cell_ID++)
	{
		angle_of_arrival[pos_cell_ID]  = new double *[nrCells];
		delay[pos_cell_ID]  = new double *[nrCells];
		for (int signal_cell_ID = 0; signal_cell_ID < nrCells; signal_cell_ID++)
		{
			angle_of_arrival[pos_cell_ID][signal_cell_ID] = new double[FADING_PATHS];
			delay[pos_cell_ID][signal_cell_ID] = new double[FADING_PATHS];
			for(int pathID=0; pathID<FADING_PATHS; pathID++)
			{
				angle_of_arrival[pos_cell_ID][signal_cell_ID][pathID] =0;
				delay[pos_cell_ID][signal_cell_ID][pathID] = 0;
			}
		}
	}

	for (int pos_cell_ID = 0; pos_cell_ID < nrMSs; pos_cell_ID++)
	{
		for (int signal_cell_ID = 0; signal_cell_ID < nrCells; signal_cell_ID++)
		{
			for (int i = 0; i < FADING_PATHS; ++i)
			{
				//for(int msID=0; msID<nrMSs; msID++) //msID is a unique value
				//{
				//angle of arrival on path i, for channel rbID used for doppler_shift calculation
				//might be also subband independent, i.e. per s
				//angle_of_arrival[pos_cell_ID][signal_cell_ID][i] = cos(uniform(0, M_PI, myID+RANDOM_SEED_ID*nrMSs));
				angle_of_arrival[pos_cell_ID][signal_cell_ID][i] = cos( M_PI / (double)FADING_PATHS * (double)i);
				//std::cout << "Angle: " << angle_of_arrival[pos_cell_ID][signal_cell_ID][i] << std::endl;
				//delay on path i
				//might be also subband independent, i.e. per s
				delay[pos_cell_ID][signal_cell_ID][i] = (double)exponential(DELAY_RMS, myBsID + RANDOM_SEED_ID);
				//std::cout << "Delay: " << delay[pos_cell_ID][signal_cell_ID][i] << std::endl;

				//u=(i+1)/(double)(FADING_PATHS+1);
				//x=-1*DELAY_RMS*log(1-u);

				//delay[pos_cell_ID][signal_cell_ID][i] = (double)exponential(DELAY_RMS, 0);

				// fprintf(resultFile1, "%d:%d:%d:%lf \n", pos_cell_ID, signal_cell_ID, i, delay[pos_cell_ID][signal_cell_ID][i]);
				//}
			}
		}
	}

	if(myBsID == 3){
		ofstream output;
		ofstream output2;
		output2.open ("Eval_PL.txt");
		for(int i = 0; i < nrMSs; i++){
			string out = "validation_" + to_string(static_cast<long long> (i)) + (string) ".txt";
			output.open (out);
			output << "Basestation: " << myBsID << endl;
			output << "Basestation Position: " << bsXposition[myBsID] << " " << bsYposition[myBsID] << endl;
			output << "Mobilestation: " << i << endl;
			output << "Mobilestation Position: " << msXposition[i] << " " << msYposition[i] << endl;
			output << "Distance: " << sqrt( pow(bsXposition[myBsID] - msXposition[i],2) + pow(bsYposition[myBsID] - msYposition[i],2) ) << endl;
			output << "Velocity: " << velocity[i]<< endl;
			output.close();
			
			double dist = sqrt( pow(bsXposition[3] - msXposition[i],2) + pow(bsYposition[3] - msYposition[i],2) );
			double dist1 = sqrt( pow(bsXposition[0] - msXposition[i],2) + pow(bsYposition[0] - msYposition[i],2) );
			double dist2 = sqrt( pow(bsXposition[1] - msXposition[i],2) + pow(bsYposition[1] - msYposition[i],2) );
			double dist3 = sqrt( pow(bsXposition[2] - msXposition[i],2) + pow(bsYposition[2] - msYposition[i],2) );
			double dist4 = sqrt( pow(bsXposition[4] - msXposition[i],2) + pow(bsYposition[4] - msYposition[i],2) );
			double dist5 = sqrt( pow(bsXposition[5] - msXposition[i],2) + pow(bsYposition[5] - msYposition[i],2) );
			double dist6 = sqrt( pow(bsXposition[6] - msXposition[i],2) + pow(bsYposition[6] - msYposition[i],2) );			
			output2 << (i+1) << " " << msXposition[i] << " " << msYposition[i] << " " << pow(10,(calcPathloss(dist)/10)) << " " << pow(10,(calcPathloss(dist1)/10)) << " " << pow(10,(calcPathloss(dist2)/10)) << " " << pow(10,(calcPathloss(dist3)/10)) << " " << pow(10,(calcPathloss(dist4)/10)) << " " << pow(10,(calcPathloss(dist5)/10)) << " " << pow(10,(calcPathloss(dist6)/10)) << " " << endl;
		}
		output2.close();
	}
	
	return true;
}

// It may be necessary for the Channel to receive Message from other LPs.
void ChannelAlternative::handleMessage(cMessage* msg){
	// not needed here.
}



// Computes the pathloss for a given distance using an arbitrary model.
double ChannelAlternative::calcPathloss(double dist){
	//cout << "Tenlogk: " << TENLOGK << " Alpha: " << ALPHA << endl;

	//double pathloss(tenlogk - 10 * alpha * log10(dist));
	//cout << "Pathloss: " << pathloss << endl;

	// This makes sure that the pathloss is never greater than 1 (Can be a problem with very close MS)
	// Note that the MS should not be that close to the BS anyway....
	//return (pathloss > 1)? 1 : pathloss;
	
	return -1.0 * (22*log10(dist) + 28 + 20*log10(2.5));
}

// Computes the Termal Noise.
double ChannelAlternative::getTermalNoise(double temp, double bandwidth){
	// 112dbm according to Paper.
	//return -142.0;
	return toDB(5.88e-15);
}

// Calculates the current SINR for given interferers.
double ChannelAlternative::calcSINR(int RB, vector<double> &power, vector<Position> &pos, vector<int> &bsId_, bool up, int msId){
	// Calc Signal Power
	double signal = calculateLoss(msId, 2500000000 + RB*180000, velocity[msId], myBsID);
	// Calc Interferer and Noise Power
	double interferer = 0;
	//cout << "Inteferer Size: " << power.size() << endl;
	for(uint i = 0; i < power.size() - 2; i++){
		interferer = interferer + toLinear(calculateLoss(msId, 2500000000 + RB*180000, velocity[msId], bsId_[i + 2]));
	}
	interferer = toDB(interferer + toLinear(getTermalNoise(300,180000)));
	//cout << "interferer Power: " << interferer << endl;
	//interferer = interferer + getTermalNoise(300, 180000);
	// Calc Ratio / SINR
	//if(myBsID == 3)
	//cout << "SINR: " << signal - interferer << endl;
	return (signal - interferer);
}

// Calculates the current SINR for given interferers.
vec ChannelAlternative::calcSINR(vector<double> &power, vector<Position> &pos, vector<int> &bsId_, bool up, int msId){
	vec sinr(25);
	if(myBsID == 3){
		//cout << "Basestation: " << myBsID << endl;
		//cout << "Basestation Position: " << bsXposition[myBsID] << " " << bsYposition[myBsID] << endl;
		//cout << "Mobilestation: " << msId << endl;
		//cout << "Mobilestation Position: " << msXposition[msId] << " " << msYposition[msId] << endl;
		//cout << "Distance: " << sqrt( pow(bsXposition[myBsID] - msXposition[msId],2) + pow(bsYposition[myBsID] - msYposition[msId],2) ) << endl;
		//cout << "SINR: " << sinr << endl;
		for(int i = 0; i < 25; i++){
			sinr.set(i,calcSINR(i, power, pos, bsId_, up, msId));
		}
		ofstream output;
		string out = "validation_" + to_string(static_cast<long long> (msId)) + (string) ".txt";
		output.open (out, fstream::app);
		if(simTime() > 2.0){
			output << simTime() << ": " << sinr << endl << endl;
		}
		output.close();
	}else{
		sinr = ones(25);
	}
	return sinr;
}

// Updates the Channel if necessary for moving MS
void ChannelAlternative::updateChannel(Position** msPos){
	int fromBsId = neighbourIdMatching->getDataStrId(myBsID);
	for(int i = 0; i < nrMSs; i++){
		msXposition[i] = msPos[fromBsId][i].x;
		msYposition[i] = msPos[fromBsId][i].y;
	}
}

//Destructor
ChannelAlternative::~ChannelAlternative(){
	delete[] delay;
	delete[] angle_of_arrival;
}

double ChannelAlternative::calculateLoss(int msId, double frequency, double mobile_speed, int intbsID){
	double dist = sqrt( pow(bsXposition[intbsID] - msXposition[msId],2) + pow(bsYposition[intbsID] - msYposition[msId],2) );
	double fading = calculateFading(msId, frequency, mobile_speed, intbsID);
	//cout << "Fading: " << fading << " Pathloss: " << calcPathloss(dist) << endl;
	//return (fading + calcPathloss(dist));
	return toDB(exponential(pow(10, calcPathloss(dist)/10), RANDOM_SEED_ID));
}

// ============================================================================================
/*
 * Jakes-like  method, Frequency in HERTZ!
 * With OFDM subbands (numbered from low to high f):
 * frequency = CENTER_FREQUENCY - SUBBANDS * FREQUENCY_SPACING / 2 + sb * FREQUENCY_SPACING + FREQUENCY_SPACING / 2
 */
double ChannelAlternative::calculateFading(int msId, double frequency, double mobile_speed, int intbsID){
	double phi_d = 0;
	double phi_i = 0;
	double phi = 0;
	//  double phi_sum = 0;

	double re_h = 0;
	double im_h = 0;
	
	double LIGHTSPEED = 299792458.0;

	double doppler_shift = mobile_speed * (frequency) / LIGHTSPEED;
	
	//cout << "Calc Fading for BS: " << myBsID << " for signal source: " << intbsID << " with Shift: " << doppler_shift << endl;

	//Although the MS might not move, there are other moving objects (reflectors) in the environment
	//cause a minimal doppler frequency
	if(doppler_shift < stationaryDoppler){
		doppler_shift = stationaryDoppler;
	}
	//cout << " Maximum Doppler Frequency is: "<< doppler_shift << " at speed " << mobile_speed << endl;

	for (int i = 0; i < FADING_PATHS; i++) {
		//some math for complex numbers:
		//z = a + ib        cartesian form
		//z = p * e ^ i(phi)    polar form
		//a = p * cos(phi)
		//b = p * sin(phi)
		//z1 * z2 = p1 * p2 * e ^ i(phi1 + phi2)

		//msAssocList[msID] gives the BS that MS is associated to
		phi_d = angle_of_arrival[msId][intbsID][i] * doppler_shift;    // phase shift due to doppler => t-selectivity
		phi_i = delay[msId][intbsID][i] * frequency;                   // phase shift due to delay spread => f-selectivity

		phi = 2.00 * M_PI * (phi_d * simTime().dbl() - phi_i);    // calculate resulting phase due to t-and f-selective fading
		//phi = 2.00 * M_PI * phi_d * simTime() - phi_i;

		//one ring model/Clarke's model plus f-selectivity according to Cavers:
		//due to isotropic antenna gain pattern on all paths only a^2 can be received on all paths
		//since we are interested in attenuation a:=1, attenuation per path is then:
		double attenuation = (1.00 / sqrt((double)FADING_PATHS));
		//cout<<"phi_d :"<<phi_d<<"phi_i :"<<phi_i<<endl;
		//convert to cartesian form and aggregate {Re, Im} over all fading paths
		re_h = re_h + attenuation * cos(phi);
		im_h = im_h - attenuation * sin(phi);
	}

	//output: |H_f|^2= absolute channel impulse response due to fading in dB
	//note that this may be >0dB due to constructive interference
	//cout << "fadingloss :" << 10 * log10(re_h * re_h + im_h * im_h) << endl;
	return 10 * log10(re_h * re_h + im_h * im_h);
	//return 10 * log10(exponential(1, myID));

}

// ===============================================================================
// Convert the double value val to its dB value.
double ChannelAlternative::toDB(double val)
{
	// Compute and return the approriate dB value.
	return (10*log10(val));
}

// ===============================================================================
// Convert the double value val to its Linear value.
double ChannelAlternative::toLinear(double val)
{
	// Compute and return the approriate Linear value.
	return pow(10,val/10);
}
