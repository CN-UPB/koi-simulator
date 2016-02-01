/*
 * Cost2100Channel.cc
 *
 *  Created on: Jul 15, 2014
 *      Author: Thomas Prinz
 * 
 * Encapsulated Implementation of the Cost 2100 Channel Model.
 */

#include <algorithm> 
#include "Cost2100Channel.h"

#if defined(__INTEL_COMPILER)
    #define thread_local __thread
#endif

//using namespace std;

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

/*
 * Distance between two points in 3D Space
 */
inline double dist(vec const &pos1, vec const &pos2){
	return sqrt(pow(pos1(0) - pos2(0),2) + pow(pos1(1) - pos2(1),2) + pow(pos1(2) - pos2(2),2));
}

bool Cost2100Channel::init(cSimpleModule* module, Position** msPositions, std::map <int,Position> neighbourPositions){
	
	//cout << "Start Initialization of Channel: " << module->getParentModule()->getParentModule()->getIndex() << endl;
	
	//! BEGIN Read in all necessary Parameters
	
	double power_factor_mean = module->par("factorLOS_mean");
	double power_factor_std = module->par("factorLOS_std");
	maxNumberOfNeighbours = module->par("maxNumberOfNeighbours");
	power_factor = power_factor_mean * pow(10, (randn() * power_factor_std / 10));
	LOS_cutoff_distance = module->par("LOS_cutoff_distance");
	VR_radius_los = module->par("VR_radius_los");
	TR_radius_los = module->par("TR_radius_los");
	vrRadius = module->par("vrRadius");
	trRadius = module->par("trRadius");
	tx_num = module->par("tx_num");
	rx_num = module->par("rx_num");
	startfrequency = module->par("startfrequency");
	endfrequency = module->par("endfrequency");
	carrierFreq = (startfrequency + endfrequency) / 2;
	subcarrier_bandwidth = module->par("subcarrier_bandwidth");
	excess_delay = module->par("excessDelay");
	cluster_power = module->par("clusterPower");
	bsId = module->par("bsId");
	tti = module->par("tti");
	double averageFarCluster = module->par( "averageFarCluster" );
    double cellRadius = module->par("cellRadius");
    numberOfMobileStations = module->par("numberOfMobileStations");
	int numberBs = module->getParentModule()->getParentModule()->getParentModule()->par("numberOfCells");
	lambda = 299792458.0 / ((startfrequency + endfrequency) / 2.0);
	initModule = module;
	neighbourPos = neighbourPositions;
	
	// Each Mobilestation need 2 Slots in this vector.
	// senderPosition.reserve(numberOfMobileStations*20);
	// targetPosition.reserve(numberOfMobileStations*20);
	// interfererPositions.reserve(numberOfMobileStations*20);
	// interfererPower.reserve(numberOfMobileStations*20);
	// senderPower.reserve(numberOfMobileStations*20);
	
	// We have already our own init Information. (therefore = 1)
	init_counter = 1;
	
	//find the neighbours and store the pair (bsId, position in data structures) in a map
	NeighbourIdMatching *neighbourIdMatching;
    cModule *cell = module->getParentModule()->getParentModule();
    neighbourIdMatching = new NeighbourIdMatching(bsId, maxNumberOfNeighbours, cell);
    numberOfNeighbours = neighbourIdMatching->numberOfNeighbours();
    
    //! END Read in all necessary Parameters
	
	//-----------------------------------
	
	//! BEGIN Generate VR's for Basestation
	
	//cout << "Start Generation of VRs for BS: " << bsId << endl;
        
	// The avg number of VR necessary to achieve on avg 'poissonCluster' active clusters
	int VR_avg = averageFarCluster * (PI * pow(cellRadius,2) / (PI * pow((vrRadius - trRadius),2)));
	// Final randomized number of VR's for this Basestation
	int VR_rand = binomial(VR_avg * numberBs,1.0/numberBs);
	//std::cout << "Number of VR's: " << VR_rand << " for BS " << bsId << std::endl;
	
	numberOfVR = VR_rand;
	
	double xPosition = module->par("xPos");
	double yPosition = module->par("yPos");
	double zPosition = module->par("zPos");
	
	neighbourPositions[bsId].x = xPosition;
	neighbourPositions[bsId].y = yPosition;
	
	// Uniformly distribute VR's within Cell
	for(int i = 0; i < VR_rand;i++){
		double x_pos, y_pos;
		do{
			x_pos = uniform(-1.0 * cellRadius,cellRadius);
			y_pos = uniform(-1.0 * cellRadius,cellRadius);
		}while(sqrt(pow(x_pos,2) + pow(y_pos,2)) > cellRadius);
		
		VisibilityRegion vr;
		vr.x = xPosition + x_pos;
		vr.y = yPosition + y_pos;
		vr.clusterId = -1;
		vr.last = false;
		vr.origin = bsId;
		VR[-1].push_back(vr);
	}
	
	// Generate Line of Sight VR for this BS
	double x, y;
	do{
		x = randu() * (LOS_cutoff_distance - VR_radius_los);
		y = randu() * (LOS_cutoff_distance - VR_radius_los);
	}while(sqrt(pow(x, 2) + pow(y, 2)) > LOS_cutoff_distance);
	
	LOS_VR.x = x + xPosition;
	LOS_VR.y = y + yPosition;
	
	//cout << "Finished VRs for BS: " << bsId << endl;
	
	//! END Generate VR's for Basestation
	
	//-----------------------------------
	
	//! BEGIN Generate local cluster for Basestation
	{
		//cout << "Start Generation of Local BS Cluster for BS: " << bsId << endl;
			
		// Cluster is positioned at the BS position.
		vec cluster_pos;
		cluster_pos.ins(0,xPosition);
		cluster_pos.ins(1,yPosition);
		cluster_pos.ins(2,zPosition);
		
		// The spread of clusters (delay, angles) and the shadow fading are
		// correlated "depending on the scenario". (Variable Cross Correlation)
		
		vector<double> correlation; 
		const char *crossCorrelation = module->par("crossCorrelation");
		mat crossCor = chol(mat(crossCorrelation));					// itpp converts string to matrix and does cholesky decomposition
		crossCor = randn(1, 6) * crossCor;
		
		// Shadow fading 's' follows a log normal distribution which means
		// for mean mu, standard derivation sigma and normal distributed
		// random variable z that s = e^(mu + sigma * z) (page 131)
		double sf_mu = module->par("shadowFading_mean");
		double sf_sigma = module->par("shadowFading_std");
		double shadow_fading = pow(M_E, (sf_mu + sf_sigma * crossCor.get(0,0)));
		
		// The cluster spread in delay (tau) as well as elevation (theta)
		// and azimuth (phi). They follow a log normal distribution with base 10
		// which results, for mean mu, standard derivation sigma and a normal
		// variable z, in x = mu * 10^(0.1 * z * sigma) (page 132)
		double delay_mu = module->par("delay_mean");
		double delay_sigma = module->par("delay_std");
		
		double tau_C = delay_mu * pow(10, (0.1 * delay_sigma * crossCor.get(0,1))); 	// delay spread
		double d_tau = (tau_C * speedOfLightVac) / 2.0;												// spatial delay spread
		
		// Link Delay is 0, as this is no Twin Cluster
		double tau_link_delay = 0;
		
		// The spatial cluster spread is, in this case (local cluster) determined
		// by the random delay. a = b = delay * c0 / 2 = d_tau (page 120)
		vec spread;
		spread.ins(0,d_tau);
		spread.ins(1,d_tau);
		spread.ins(2,0);		// This is set to d_tau / 20 in the other implementation. It should actually be 0 in my opinion. TODO: verify
		spread.ins(3,d_tau);
		spread.ins(4,d_tau);
		spread.ins(5,0);		// This is set to d_tau / 20 in the other implementation. It should actually be 0 in my opinion. TODO: verify
		
		// This is no Twin Cluster, so the position is the same for both
		localcluster_bs = Cluster(LOCAL_AT_BS,cluster_pos,cluster_pos,spread,shadow_fading,tau_link_delay,zeros(4), -1);
		
		//cout << "Finished Local BS Cluster for BS: " << bsId << endl;
	}
	//! END Generate local cluster for Basestation
	
	//-----------------------------------
	
	//! BEGIN Generate local cluster for all Mobile stations within Cell
	{
		//cout << "Startet Initialization of Local MS clusters for BS: " << bsId << endl;
		
		int fromBsId = neighbourIdMatching->getDataStrId(bsId);
		
		//std::cout << "MS Pos: " << msPositions[fromBsId][i].x << " " << msPositions[fromBsId][i].y << std::endl;
		
		// Cluster is positioned at the MS position.
		vec cluster_pos;
		cluster_pos.ins(0,msPositions[fromBsId][module->getIndex()].x);
		cluster_pos.ins(1,msPositions[fromBsId][module->getIndex()].y);
		cluster_pos.ins(2,0);
		
		// The spread of clusters (delay, angles) and the shadow fading are
		// correlated "depending on the scenario". (Variable Cross Correlation)
		
		vector<double> correlation; 
		const char *crossCorrelation = module->par("crossCorrelation");
		mat crossCor = chol(mat(crossCorrelation));					// itpp converts string to matrix and does cholesky decomposition
		crossCor = randn(1, 6) * crossCor;
		
		// Shadow fading 's' follows a log normal distribution which means
		// for mean mu, standard derivation sigma and normal distributed
		// random variable z that s = e^(mu + sigma * z) (page 131).
		// Other implementation takes 10 as a base here...
		double sf_mu = module->par("shadowFading_mean");
		double sf_sigma = module->par("shadowFading_std");
		double shadow_fading = pow(M_E, (sf_mu + sf_sigma * crossCor.get(0,0)));
		
		// The cluster spread in delay (tau) as well as elevation (theta)
		// and azimuth (phi). They follow a log normal distribution with base 10
		// which results, for mean mu, standard derivation sigma and a normal
		// variable z, in x = mu * 10^(0.1 * z * sigma) (page 132)
		double delay_mu = module->par("delay_mean");
		double delay_sigma = module->par("delay_std");
		
		double tau_C = delay_mu * pow(10, (0.1 * delay_sigma * crossCor.get(0,1))); 	// delay spread
		double d_tau = (tau_C * speedOfLightVac) / 2.0;									// spatial delay spread
		
		// Link Delay is 0, as this is no Twin Cluster
		double tau_link_delay = 0;
		
		// The spatial cluster spread is, in this case (local cluster) determined
		// by the random delay. a = b = delay * c0 / 2 = d_tau (page 120)
		vec spread;
		spread.ins(0,d_tau);
		spread.ins(1,d_tau);
		spread.ins(2,0);		// This is set to d_tau / 20 in the other implementation. It should actually be 0 in my opinion.
		spread.ins(3,d_tau);
		spread.ins(4,d_tau);
		spread.ins(5,0);		// This is set to d_tau / 20 in the other implementation. It should actually be 0 in my opinion.
		
		// This is no Twin Cluster, so the position is the same for both
		localclusters_ms.push_back(Cluster(LOCAL_AT_MS,cluster_pos,cluster_pos,spread,shadow_fading,tau_link_delay,zeros(4), -1));
		
		//cout << "Finished Initialization of Local MS clusters for BS: " << bsId << endl;
	}
	//! END Generate local cluster for all Mobile stations within Cell
	
	//-----------------------------------
	
	//! BEGIN Assign Basestations to the Visibility Regions.
	
	//cout << "Startet assignment of the BS to the VR. (local for BS " << bsId << ")" << endl;
        cout << "NumVR:" << numberOfVR << endl;	
	// Essentially building Commom Cluster Assignment table here.
	for(int i = 0; i < numberOfVR; i++){
		for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
			
			// TODO: 3D Position exchange
			vec neighbourPos;
			neighbourPos.ins(0,it->second.x);
			neighbourPos.ins(1,it->second.y);
			neighbourPos.ins(2,0);
			
			vec own_pos;
			own_pos.ins(0,xPosition);
			own_pos.ins(1,yPosition);
			own_pos.ins(2,zPosition);
			
			double dist_bs_neighbour = dist(own_pos, neighbourPos);
			
			//std::cout << "Distance between Neighbour and BS: " << dist_bs_neighbour << " ( Cellradius: " << cellradius << " )" << std::endl;
			
			// TODO: get probability for different ranges from input parameter
			// 0.01% off from actual probability because VR cannot be assigned to no BS.
			if((uniform(0,1) < 14.0/17.0) && (dist_bs_neighbour < 2 * cellRadius)){
				VR.at(-1)[i].bsIds.push_back(it->first);
			}
		}
	}
	
	//cout << "Finished assignment of the BS to the VR. (local for BS " << bsId << ")" << endl;
	
	//! END Assign Basestations to the Visibility Regions.
	
	//-----------------------------------
	
	//! BEGIN Choose BS that takes care of Cluster generation for the VR
	
	//cout << "Startet assigning the VR to a BS that generates the Remote Cluster for BS: " << bsId << endl;
	
	// Choose a Hop Distance for each VR. Then randomly choose a Neighbour
	// with that particular distance to take care of the cluster generation
	ivec hopDistance;
	double cRadius = module->par("cellRadius");
	
	// TODO: get from input parameter (currently 66%/34% 0/1 hop)
	for(int i = 0; i < numberOfVR; i++){
		double random_hop_distance = uniform(0,1);
		if(random_hop_distance < 0.66){
			hopDistance.ins(i,0);
		}else{
			hopDistance.ins(i,1);
		}
	}
	
	vec bs_pos(3);
	bs_pos.set(0,xPosition);
	bs_pos.set(1,yPosition);
	bs_pos.set(2,zPosition);		// Currently 0
	
	// Get 1 hop / 2 hop Neighbourhood
	std::vector<int> oneHop;
	std::vector<int> twoHop;
	for(std::map<int, Position>::iterator it = neighbourPositions.begin(); it != neighbourPositions.end(); it++){
		vec neighbourPos;
		neighbourPos.ins(0,it->second.x);
		neighbourPos.ins(1,it->second.y);
		neighbourPos.ins(2,0);
		double dist_temp = dist(bs_pos,neighbourPos);
		
		if(dist_temp == 0){
			continue;
		} else if(dist_temp < 2 * cRadius) {
			oneHop.push_back(it->first);
		} else {
			twoHop.push_back(it->first);
		}
	}
	
	for(int i = 0; i < hopDistance.size(); i++){
		VR.at(-1).at(i).assigner = -1;
		switch(hopDistance(i)){
			case 0:
				VR.at(-1).at(i).assigner = bsId;
				//ForeignVR.push_back(VR.at(i));
				break;
			case 1:
				VR.at(-1).at(i).assigner = oneHop.at(intuniform(0,oneHop.size() - 1));
				break;
			case 2:
				VR.at(-1).at(i).assigner = twoHop.at(intuniform(0,twoHop.size() - 1));
				break;
			default:
				cout << "reached default case!" << endl;
		}
	}
	
	//cout << "Finished assigning the VR to a BS that generates the Remote Cluster for BS: " << bsId << endl;
	
	//! END Choose BS that takes care of Cluster generation for the VR
	
	//-----------------------------------
	
	//! BEGIN Send VR's to the other Basestations
	
	// Send all VR's around, so that the assigner can assign it to some 
	// cluster. Send an empty VR Message with last=true to indicate that
	// this Basestations has sent everything.
	for(int i = 0; i < numberOfVR; i++){
		VisibilityRegionMessage *vrMessage = new VisibilityRegionMessage("CHANNEL_INFO",1);
		vrMessage->setVr(VR.at(-1).at(i));
		// We still have the module pointer and can therefore send stuff around
		// Messages called "CHANNEL_INFO" will be automatically forwarded to the Channel
		// and are processed there.
		module->send(vrMessage, "toPhy");
	}
	VisibilityRegion vrTemp;
	vrTemp.last = true;
	vrTemp.assigner = -1;
	
	VisibilityRegionMessage *vrMessage = new VisibilityRegionMessage("CHANNEL_INFO",1);
	vrMessage->setVr(vrTemp);
	module->send(vrMessage, "toPhy");
	
	//! End Send VR's to the other Basestations
	
	return true;
}

void Cost2100Channel::handleMessage(cMessage* msg){
	//cout << "Received Message in encapsulated Channel class." << endl;
	if(msg->isName("CHANNEL_INFO")){
		// We receive the Foreign VRs. (Message kind type 2)
		if(msg->getKind() == 2){
			VisibilityRegion vr = ((VisibilityRegionMessage *) msg)->getVr();
			
			//cout << "Assigner is: " << vr.assigner << " at BS: " << bsId << " with coords: " << vr.x << " " << vr.y << endl;
			
			if(vr.last == true){
				// Count the number of Neighbours who have transmitted all necessary data to you
				init_counter++;
			}
			
			// Received all information to start generate clusters
			if(init_counter == numberOfNeighbours && vr.last == true){
				//std::cout << "NEW CHANNEL: Received all Info! Received " << ForeignVR.size() << " at BS " << bsId << " different VR's" << std::endl;
				initModule->scheduleAt(simTime() + 100 * tti, new cMessage("CHANNEL_INFO", 4));
			}
			
			// Only Accept the VR as Foreign VR if this BS is the Assginer BS
			// TODO: only send the message to the right BS in the first place
			if(vr.assigner == bsId){
				ForeignVR.push_back(vr);
			}
		}else if(msg->getKind() == 4){
			//cout << "Start Generating remote Clusters for BS: " << bsId << endl;
			
			//! BEGIN Generate Remote Clusters for the received Visibility Regions
			
			// Average number of Visibility Regions per Cluster
			double cc = initModule->par("commonCluster");
			
			// TODO: Cluster Ids must be unique within entire simulation.
			// Works as far as less than 10000 Cluster per BS is guaranteed
			int clusterId = bsId*10000;
			int numCluster = 0;
			int numVr = ForeignVR.size();
			
			vector<VisibilityRegion> temporyVR;
			std::map<int, std::vector<VisibilityRegion> > commonClusterTable;
			
			// As long as there are VR without a cluster, continue assigning.
			while(numVr > 0){
				double vr_group = poisson(cc);
				while(vr_group > ForeignVR.size() || vr_group == 0){
					vr_group = poisson(cc);
				}
				
				for(int i = 0; i < vr_group; i++){
					int randVr = uniform(0,ForeignVR.size() - 1);
					commonClusterTable[clusterId].push_back(ForeignVR.at(randVr));
					ForeignVR.at(randVr).clusterId = clusterId;
					temporyVR.push_back(ForeignVR.at(randVr));
					ForeignVR.erase(ForeignVR.begin() + randVr);
					numVr--;
				}
				clusterId++;
				numCluster++;
			}
			
			//std::cout << "Number Cluster: " << numCluster << std::endl;
			
			// ForeignVR is now empty. Set it back to the set of VR that are now assigned
			ForeignVR = temporyVR;
			
			// Delay Parameters
			double delay_mu = initModule->par("delay_mean");
			double delay_sigma = initModule->par("delay_std");
			
			// Get the Azimuth and Elevation Spread for Both sides (mean / std) in degree
			double mu_theta_BS = initModule->par("EoD_mean");		// mean Elevation Spread at BS side
			double mu_phi_BS = initModule->par("AoD_mean");			// mean Azimuth Spread at BS side
			double mu_theta_MS = initModule->par("EoA_mean");		// mean Elevation Spread at MS side
			double mu_phi_MS = initModule->par("AoA_mean");			// mean Azimuth Spread at MS side
			
			double sigma_theta_BS = initModule->par("EoD_std");		// std Elevation Spread at BS side
			double sigma_phi_BS = initModule->par("AoD_std");		// std Azimuth Spread at BS side
			double sigma_theta_MS = initModule->par("EoA_std");		// std Elevation Spread at MS side
			double sigma_phi_MS = initModule->par("AoA_std");		// std Azimuth Spread at MS side
			
			// Shadow Fading Parameters.
			double sf_mu = initModule->par("shadowFading_mean");
			double sf_sigma = initModule->par("shadowFading_std");
			
			std::map<int, std::vector<VisibilityRegion> >::iterator it = commonClusterTable.begin();
					
			for(int i = 0;i < numCluster; i++){
				// Compute the reference VR and reference BS
				std::vector<VisibilityRegion> VRTable = it->second;
				int numberBs = initModule->getParentModule()->getParentModule()->getParentModule()->par("numberOfCells");
				int BSCount[numberBs];
				
				for(int j = 0; j < numberBs; j++){
					BSCount[j] = 0;
				}
				
				// For this particular cluster, add up all BS of all VR it is connected with.
				// The BS which appears the most will be the reference BS.
				for(uint j = 0; j < VRTable.size(); j++){
					for(uint k = 0; k <VRTable.at(j).bsIds.size(); k++){
						BSCount[VRTable.at(j).bsIds.at(k)]++;
					}
				}
				
				int max = 0;
				int maxId = 0;
				for(int j = 0; j < numberBs; j++){
					if(BSCount[j] > max){
						max = BSCount[j];
						maxId = j;
					}
				}
				//std::cout << "Max: " << max << " at " << maxId << " (For Cluster " << (i + 10000*bsId) << ")" << std::endl;
				if(max == 0){
					maxId = bsId;
				}
				
				uint maxVR = 0;
				int VR_idx = 0;
				
				for(uint j = 0; j < VRTable.size(); j++){
					bool look = false;
					for(auto it = VRTable.at(j).bsIds.begin(); it != VRTable.at(j).bsIds.end(); it++){
						if (maxId == *it){
							look = true;
						}
					}
					if(!look){
						continue;
					}
					if(VRTable.at(j).bsIds.size() > maxVR){
						maxVR = VRTable.at(j).bsIds.size();
						VR_idx = j;
					}
				}

				vector<double> correlation; 
				const char *crossCorrelation = initModule->par("crossCorrelation");
				// itpp converts string to matrix and does cholesky decomposition
				mat crossCor = chol(mat(crossCorrelation));
				crossCor = randn(1, 6) * crossCor;
				
				// Shadow fading 's' follows a log normal distribution which means
				// for mean mu, standard derivation sigma and normal distributed
				// random variable z that s = e^(mu + sigma * z) (page 131).
				// Other implementation takes 10 as a base here...
				double shadow_fading = pow(M_E, (sf_mu + sf_sigma * crossCor.get(0,0)));
				
				// The cluster spread in delay (tau) as well as elevation (theta)
				// and azimuth (phi). They follow a log normal distribution with base 10
				// which results, for mean mu, standard derivation sigma and a normal
				// variable z, in x = mu * 10^(0.1 * z * sigma) (page 132)
				double c0 = 299792458; // Speed of light
				double tau_C = delay_mu * pow(10, (0.1 * delay_sigma * crossCor.get(0,1))); 	// delay spread
				double d_tau = (tau_C * c0) / 2.0;												// spatial delay spread
				
				double theta_C_BS = mu_theta_BS * pow(10, (0.1 * sigma_theta_BS * crossCor.get(0,2))); 		// elevation spread to BS
				double phi_C_BS = mu_phi_BS * pow(10, (0.1 * sigma_phi_BS * crossCor.get(0,3))); 			// azimuth spread to BS
				double theta_C_MS = mu_theta_MS * pow(10, (0.1 * sigma_theta_MS * crossCor.get(0,4))); 		// elevation spread to MS
				double phi_C_MS = mu_phi_MS * pow(10, (0.1 * sigma_phi_MS * crossCor.get(0,5))); 			// azimuth spread to MS
				vec VR_position = zeros(3);
				vec BS_position = zeros(3);
				
				//std::cout << "BS chosen as ref BS: " << maxId << std::endl; 
				
				// Get Reference BS Positions:
				BS_position.set(0,neighbourPos.at(maxId).x);
				BS_position.set(1,neighbourPos.at(maxId).y);
				BS_position.set(2,0);

				// Get Reference VR Positions:
				VR_position.set(0,VRTable.at(VR_idx).x);
				VR_position.set(1,VRTable.at(VR_idx).y);
				VR_position.set(2,0);
				
				// TODO: 3D
				vec BS_to_VR = (VR_position - BS_position);
				BS_to_VR.set(2, 0);
				
				// Get Spherical Coordinates
				vec BS_to_VR_Sph = Cart_to_Sph(BS_to_VR);

				// The vector is rotated with a Gaussian distributed angle relative to the
				// imaginary line between the BS and the VR (Page 121)
				// The distance follows a nonnegative distribution (e.g. exponential for macro cell)
				// bounded by a minimum distance r_min. (page 121)
				// The other implementation chose a uniform distribution for distance instead.
				vec BS_to_VR_Sph_Final;
				double azimuth_base = initModule->par("cluster_azimuth_base");
				double azimuth_range = initModule->par("cluster_azimuth_range");
				double elevation_base = initModule->par("cluster_elevation_base");
				double elevation_range = initModule->par("cluster_elevation_range");
				BS_to_VR_Sph_Final.ins(0, BS_to_VR_Sph.get(0) + azimuth_base + randn() * azimuth_range);
				BS_to_VR_Sph_Final.ins(1, BS_to_VR_Sph.get(1) + elevation_base + randn() * elevation_range);
				
				double distance_base = initModule->par("cluster_distance_base");
				double distance_range = initModule->par("cluster_distance_range");
				double dist_C_BS = exponential(distance_range);
				if(dist_C_BS < distance_base){dist_C_BS = distance_base;}
				BS_to_VR_Sph_Final.ins(2, dist_C_BS);
				
				// Transform spherical Coordinates back to Cartesian
				vec cluster_pos_BS = BS_position + Sph_to_Cart(BS_to_VR_Sph_Final);
				vec cluster_pos_MS;
				
				// Determine if Twin Cluster or Single Cluster
				double K_sel = initModule->par("singleClusterRatio");
				int cluster_type = (uniform(0,1) > K_sel) ? 1 : 0;
				CLUSTER_TYPE c_type;
				
				double tau_link_delay;
				// Twin Cluster
				if(cluster_type == 1){
					// Compute the second representation for the Twin Cluster, the same way as the first.
					// The distance has to match with the angular spread.
					// Therefore: dist_C_BS * tan ( angle_C_BS ) = dist_C_VR * tan( angle_C_VR ) (Page 121)
					vec VR_cluster_vec = zeros(3);
					VR_cluster_vec.set(0, BS_to_VR_Sph.get(0) + azimuth_base + randn() * azimuth_range);
					VR_cluster_vec.set(1, BS_to_VR_Sph.get(1) + elevation_base + randn() * elevation_range);
					VR_cluster_vec.set(2, dist_C_BS * tan(phi_C_BS / 180 * pi / 2) / tan(phi_C_MS / 180 * pi / 2));
					cluster_pos_MS = VR_position + Sph_to_Cart(VR_cluster_vec);
					
					// Compute the additional delay caused by multi bounce
					// The delay is log normal distributed with base 10
					double twin_delay_mean = initModule->par("linkDelay_mean");
					double twin_delay_std = initModule->par("linkDelay_std");
					tau_link_delay = (dist(cluster_pos_BS, cluster_pos_MS) / c0) + twin_delay_mean * pow(10,(0.1 * twin_delay_std * randn()));
					
					c_type = TWIN;
				// Single Cluster
				}else{
					tau_link_delay = 0;					// Single Cluster has no link delay
					cluster_pos_MS = cluster_pos_BS;	// Position is the same for Single Cluster
					c_type = SINGLE;
				}
				
				vec spread;
				// Compute the clusters actual spatial spread
				spread.ins(0, d_tau);
				// Here dist_C_BS = Adjacent and tan ( alpha ) = Opposite / Adjacent (Trigonometry)
				// Therefore Adjacent eliminates itself and the length of the Opposite remains.
				// which is in this case the elevation/azimuth spread of the ellipsoid that forms the cluster
				spread.ins(1, dist_C_BS * tan(phi_C_BS / 180 * pi / 2));
				spread.ins(2, dist_C_BS * tan(theta_C_BS / 180 * pi / 2));
				spread.ins(3, d_tau);
				spread.ins(4, dist_C_BS * tan(phi_C_BS / 180 * pi / 2));
				spread.ins(5, dist_C_BS * tan(theta_C_BS / 180 * pi / 2));
				
				// Compute the actual angle from BS/MS towards the cluster:
				vec BS_angle = Cart_to_Sph(BS_position - cluster_pos_BS).get(0, 1);
				vec MS_angle = Cart_to_Sph(VR_position - cluster_pos_MS).get(0, 1);
				vec angle_final = zeros(4);
				angle_final.set(0,BS_angle(0));
				angle_final.set(1,BS_angle(1));
				angle_final.set(2,MS_angle(0));
				angle_final.set(3,MS_angle(1));
				
				// Save Cluster
				//remoteCluster.push_back(Cluster(c_type, cluster_pos_BS, cluster_pos_MS, spread, shadow_fading, tau_link_delay, angle_final, i + bsId*10000));
				remoteCluster[i + bsId*10000] = Cluster(c_type, cluster_pos_BS, cluster_pos_MS, spread, shadow_fading, tau_link_delay, angle_final, i + bsId*10000);
				it++;
			}
			
			//cout << "Finished Generating remote Clusters for BS: " << bsId << endl;
			
			//! END Generate Remote Clusters for the received Visibility Regions
			
			//-----------------------------------
			
			//! BEGIN Send generated Clusters around
			
			// Send VR's back to Origin with ClusterId and send the necessary Cluster directly with it
			cout << "Foreign VR size: " << ForeignVR.size() << endl;
			for(uint i = 0; i < ForeignVR.size(); i++){
				VisibilityRegion tempVR = ForeignVR.at(i);
				
				// Send back VR to its Origin
				VisibilityRegionMessage *vrMessage = new VisibilityRegionMessage("CHANNEL_INFO",5);
				vrMessage->setVr(tempVR);
				initModule->sendDelayed(vrMessage, 5 * tti , "toPhy");
			}
			
			std::cout << "Remote Cluster Size: " << remoteCluster.size() << std::endl;
			for(auto it = remoteCluster.begin(); it != remoteCluster.end(); ++it){
				Cluster tempCluster = it->second;
				
				ClusterMessage *clusterMessage = new ClusterMessage("CHANNEL_INFO",7);
				clusterMessage->setC(tempCluster);
				
				initModule->sendDelayed(clusterMessage, 10 * tti , "toPhy");
			}
			/*
			for(uint i = 0; i < remoteCluster.size(); i++){
				Cluster tempCluster = remoteCluster.at(i);
				
				ClusterMessage *clusterMessage = new ClusterMessage("CHANNEL_INFO",7);
				clusterMessage->setC(tempCluster);
				
				initModule->sendDelayed(clusterMessage, 10 * tti , "toPhy");
			}*/
			
			// After all have been sent, clear all information,
			// the necessary Info is received via the messages just sent.
			ForeignVR.clear();
			VR.clear();
			remoteCluster.clear();

			//std::cout << "sent VR's back and Clusters with it from BS" << bsId << std::endl;
			
			//! END Send generated Clusters around
		}else if(msg->getKind() == 6){
			VisibilityRegion vr = ((VisibilityRegionMessage *) msg)->getVr();
			//if(vr.origin == bsId){
				VR[vr.origin].push_back(vr);
			//}
			//cout << "Received " << VR.size() << " VRs out of " << numberOfVR << endl;
		}else if(msg->getKind() == 8){
			Cluster c = ((ClusterMessage *) msg)->getC();
			remoteCluster[c.get_cluster_ID()] = c;
			// Add cluster if it belongs to one of the VR of this Basestation
			/*
			for(uint i = 0; i < VR.size(); i++){
				if(VR.at(i).clusterId == c.get_cluster_ID()){
					remoteCluster.push_back(c);
					break;
				}
			}
			* */
		}
		delete msg;
	}else{
		// cannot process msg.
		delete msg;
	}
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
double Cost2100Channel::calcPathloss(double dist) {
	double pathloss;
	
	// For Testing the Cost 273 parameters defined on page 130 of the Cost2100 book are used.
	// The cell radius here is defined as 1000m which fits our scenario.
	// TODO: get from external
    static double h_BS = 50;			// Basestation height
    static double h_rooftop = 15;		// Rooftop height
    static double freq_c = 2e9;		// center frequency
    static double phi_road = 45;		// Road orientation
    static double w_road = 25;			// Road width
    static double w_street = 50;		// Building Seperation
    static double h_MS = 1.5;			// Mobile station height
    
	static double delta_base = h_BS - h_rooftop;
	
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

	// The squareroot of the amount of power remains.
	return pathloss;
}

// Johnson Nyquist Noise
double Cost2100Channel::getTermalNoise(double temp, double bandwidth){
	return (temp * bandwidth * 1.3806488e-23);
}

/*
 * Computes attenuation according to Formula 3.40
 * tau_0 = Line of Sight Delay
 * tau_b = Cut off Delay
 * tau = Actual delay
 */
double Cost2100Channel::attenuation(double decay_factor, double tau_0, double tau_b, double tau){
	return std::max( exp(-1.0*decay_factor*(tau - tau_0)) , exp(-1.0*decay_factor*(tau_b - tau_0)) );
}

/*
 * Computes the Final Matrix over all active clusters.
 */
vec Cost2100Channel::H(Cluster &local_bs, vector<Cluster> &local_ms, unordered_map <int, Cluster> &cluster,std::map<int,vector<VisibilityRegion>> &vr, int num_mpc, int num_dmc, int bsId, Position &ms_Pos, Position &bs_Pos, VisibilityRegion &LOS_VR, int RB){
	vector<cmat> H;
	
	vec bsPosition(3);
	bsPosition.set(0, bs_Pos.x);
	bsPosition.set(1, bs_Pos.y);
	bsPosition.set(2,0);
	
	vec msPosition(3);
	msPosition.set(0, ms_Pos.x);
	msPosition.set(1, ms_Pos.y);
	msPosition.set(2,0);
	
	//! BEGIN Generate Signal Matrix
	//std::cout << "Computing Matrix for BS(" << bsPosition(0) << "," << bsPosition(1) << ") and MS(" <<  msPosition(0) << "," << msPosition(1) << ")" << std::endl;
	
	vector<VisibilityRegion> signalVRs;
	vector<Cluster> signalCluster;
	vector<double> VRGain;
	signalVRs.reserve(100);
	signalCluster.reserve(100);
	VRGain.reserve(100);
	
	//std::cout << "vr.size() " << vr[bsId].size() << "at BS: " << bsId << std::endl;
	
	for(uint i = 0; i < vr[bsId].size(); i++){
		vec VRPos(3);
		VRPos.set(0,vr.at(bsId).at(i).x);
		VRPos.set(1,vr.at(bsId).at(i).y);
		VRPos.set(2,0);
		
		if(dist(VRPos,msPosition) < vrRadius && (std::find(vr.at(bsId).at(i).bsIds.begin(), vr.at(bsId).at(i).bsIds.end(), bsId) != vr.at(bsId).at(i).bsIds.end())){
			double vrgain = calc_VRGain(lambda,trRadius,vrRadius,VRPos,msPosition);
			signalVRs.push_back(vr.at(bsId).at(i));
			VRGain.push_back(vrgain);
		}
	}
	
	//cout << "Found VRs: " << signalVRs.size();
	
	/*
	for(uint i = 0; i < signalVRs.size(); i++){
		for(uint j = 0; j < cluster.size(); j++){
			if(signalVRs.at(i).clusterId == cluster.at(j).get_cluster_ID()){
				signalCluster.push_back(cluster.at(j));
			}
		}
	}
	*/
	for(uint i = 0; i < signalVRs.size(); i++){
		signalCluster.push_back(cluster[signalVRs.at(i).clusterId]);
	}
	
	double pathloss = calcPathloss(dist(msPosition,bsPosition));
	//std::cout << "Pathloss (Cost 231 Walfish Ikegami Model): " << pathloss << std::endl;
	
	// Antennae spacing 1/2 wavelength in z direction
	mat antennae_pos = zeros(2,3);
	antennae_pos.set(1,2,lambda * 0.5);
	
	// Contains the double directional (Angles for BS and MS) impusle response (DDIR)
	mat channel_matrix;
	
	for(uint i = 0; i < signalCluster.size(); i++){
		mat channel = compute_channel_matrix(pathloss, signalCluster.at(i), num_mpc, num_dmc, bsPosition, msPosition, carrierFreq, VRGain.at(i));
		for(int j = 0; j < num_mpc + num_dmc; j++){
			channel_matrix.append_row(channel.get_row(j));
		}
	}
	
	// Compute the local MS Cluster
	/*
	double vr_gain_local_ms = calc_VRGain(lambda,trRadius,vrRadius,msPosition,msPosition);
	mat channel = compute_channel_matrix(pathloss, local_bs, num_mpc, num_dmc, bsPosition, msPosition, carrierFreq, vr_gain_local_ms);
	for(int j = 0; j < num_mpc + num_dmc; j++){
		channel_matrix.append_row(channel.get_row(j));
	}
	* */
	
	// Compute the local BS Cluster
	/*
	double vr_gain_local_bs = calc_VRGain(lambda,trRadius,vrRadius,bsPosition,msPosition);
	channel = compute_channel_matrix(pathloss, local_ms, num_mpc, num_dmc, bsPosition, msPosition, carrierFreq, vr_gain_local_bs);
	for(int j = 0; j < num_mpc + num_dmc; j++){
		channel_matrix.append_row(channel.get_row(j));
	}
	* */
	
	// Compute Line of Sight Component
	double power_los, power_factor_los;
	double d_ms_bs = dist(msPosition,bsPosition);
	vec LOS_VR_pos(3);
	LOS_VR_pos.set(0,LOS_VR.x);
	LOS_VR_pos.set(1,LOS_VR.y);
	LOS_VR_pos.set(2,0);
	
	// Power is Zero, if MS lies outside of the Line of Sight area (distance)
	if (d_ms_bs > LOS_cutoff_distance) {
		power_los = 0;
		power_factor_los = 0;
	} else {
		double d_ms_vr_los = dist(msPosition, LOS_VR_pos);
		// Power is Zero, if MS lies outside of the Line of Sight VR
		if (d_ms_vr_los > VR_radius_los) {
			power_los = 0;
			power_factor_los = 0;
		} else {
			// The Line of Sight Power equals the sum off all cluster powers 
			// multiplied by a factor that is log normal distributed (Section 3.6.1.3)
			double power_other = 0;
			power_factor_los = power_factor;
			power_other = pow(sum(channel_matrix.get_col(5)), 2) + pow(sum(channel_matrix.get_col(6)), 2);
			power_los = abs(power_factor_los) * power_other;
		}
	}
	
	vec los_bs_sph = Cart_to_Sph(msPosition - bsPosition);
	vec los_ms_sph = Cart_to_Sph(bsPosition - msPosition);

	double d_bs_ms_xy = dist(bsPosition, msPosition);	// distance BS to MS, in X-Y plane
	double tau_0_xy = d_bs_ms_xy / speedOfLightVac;								// delay of the LOS
	complex<double> phase_los(0,  -2 * PI * carrierFreq * tau_0_xy);
	complex<double> channel_amp_los = sqrt(power_los) * exp(phase_los);
	double channel_amp_los_real = channel_amp_los.real();
	double channel_amp_los_imag = channel_amp_los.imag();
	vec channel_los = concat(los_bs_sph.get(0, 1), los_ms_sph.get(0, 1), to_vec(tau_0_xy), to_vec(channel_amp_los_real), to_vec(channel_amp_los_imag));
	channel_matrix.append_row(channel_los);
	
	//cout << "Channel Matrix: " << channel_matrix << endl;
	
	// Generate own final Matrix for each ressourceblock frequency
	// In principle this Matrix is 7-dimensional: BS X MS X Cluster X MPC X RX_Antenna X TX_Antenna X Frequency
	// Since BS, MS, are fixed here and all MPC/Cluster are added up, we end up with a 3-dimensional Matrix
	// According to Formula 4.2
	if(RB == -1){
		for(double f = startfrequency + subcarrier_bandwidth*0.5; f < endfrequency;f = f + subcarrier_bandwidth){
			cmat H_current(tx_num,rx_num);
			for(int i = 0; i < tx_num; i++){
				for(int j = 0; j < rx_num; j++){
					H_current.set(i,j,complex<double>(0,0));
				}
			}
			for (int m = 0; m < channel_matrix.rows(); m++ ) {
				complex<double> delay_phase(0, -2.0 * pi * f * channel_matrix.get(m, 4));
				complex<double> delay_response = exp(delay_phase);
				complex<double> channel_amplitude(channel_matrix.get(m, 5), channel_matrix(m, 6));
				
				cmat MIMO_mat(2,2);
				calc_MIMO_Channel_Matrix(2, 2, antennae_pos,antennae_pos, channel_matrix.get(m, 2), channel_matrix.get(m, 3), channel_matrix.get(m, 0), channel_matrix.get(m, 1), lambda, MIMO_mat);
				
				// Final formula 4.2
				H_current = H_current + (MIMO_mat * delay_response * channel_amplitude);
			}
			//std::cout << std::endl << "Composite Channel Matrix for frequency " << f << ": " << std::endl << H_current << std::endl;
			H.push_back(H_current);
		}
	}else{
		for(double f = startfrequency + subcarrier_bandwidth*(RB+0.5); f < startfrequency + subcarrier_bandwidth*(RB+1) ;f = f + subcarrier_bandwidth){
			cmat H_current(tx_num,rx_num);
			for(int i = 0; i < tx_num; i++){
				for(int j = 0; j < rx_num; j++){
					H_current.set(i,j,complex<double>(0,0));
				}
			}
			for (int m = 0; m < channel_matrix.rows(); m++ ) {
				complex<double> delay_phase(0, -2.0 * pi * f * channel_matrix.get(m, 4));
				complex<double> delay_response = exp(delay_phase);
				complex<double> channel_amplitude(channel_matrix.get(m, 5), channel_matrix(m, 6));
				
				cmat MIMO_mat(2,2);
				calc_MIMO_Channel_Matrix(2, 2, antennae_pos,antennae_pos, channel_matrix.get(m, 2), channel_matrix.get(m, 3), channel_matrix.get(m, 0), channel_matrix.get(m, 1), lambda, MIMO_mat);
				
				// Final formula 4.2
				H_current = H_current + (MIMO_mat * delay_response * channel_amplitude);
			}
			//std::cout << std::endl << "Composite Channel Matrix for frequency " << f << ": " << std::endl << H_current << std::endl;
			H.push_back(H_current);
		}
	}

	//! END Generate Signal Matrix
	vec result(H.size());
	// SISO Case Hardcoded for now.
	for(uint i = 0; i < H.size(); i++){
		result.set(i,pow((abs(H.at(i).get(0,0))),2));
		//result.set(i,abs(H.at(i).get(0,0)));
	}
	//cout << result << endl << endl;
	//cout << "result size: " << result.size() << endl;
	return result;
}

/*
 * Computes a Matrix, that contains angles, delays and amplitude for all MPC's of a given cluster.
 * Formula 3.41 / 3.42 from the book.
 */
mat Cost2100Channel::compute_channel_matrix(double pathloss, Cluster &c, int num_mpc, int num_dmc, vec const &bs_pos, vec const &ms_pos, double center_freq, double VR_gain){
	mat mpc_bs_pos = c.get_cluster_MPC().get_MPC_pos(BS);
	mat mpc_ms_pos = c.get_cluster_MPC().get_MPC_pos(MS);
	mat channel_matrix;
	//cout << c.get_cluster_pos(BS) << c.get_cluster_pos(MS) << bs_pos << ms_pos << endl;
	double cluster_delay = (dist(c.get_cluster_pos(BS), bs_pos) + dist(c.get_cluster_pos(MS), ms_pos)) / speedOfLightVac + c.get_cluster_link_delay();
	double LOS_delay = dist(ms_pos, bs_pos) / speedOfLightVac;
	
	double cluster_attenuation = attenuation(cluster_power, LOS_delay, excess_delay, cluster_delay);
	double cluster_shadow_fading = c.get_cluster_shadowing_fading();
	
	// Specular MPC
	for(int i = 0; i < num_mpc; i++){
		// Exact Delay for this multipath component.
		double mpc_delay = (dist(mpc_bs_pos.get_row(i), bs_pos) + dist(mpc_ms_pos.get_row(i), ms_pos)) / speedOfLightVac + c.get_cluster_link_delay();
		
		complex<double> mpc_phase = exp(complex<double>(0,-2.0 * M_PI * mpc_delay * center_freq));
		
		// MPC complex amplitude according to formula 3.41
		complex<double> mpc_amplitude = pathloss * c.get_cluster_MPC().get_MPC_amplitude()(i) * sqrt(cluster_attenuation*cluster_shadow_fading) * mpc_phase;
		
		vec channel = concat(Cart_to_Sph(mpc_bs_pos.get_row(i) - bs_pos).get(0, 1), Cart_to_Sph(mpc_ms_pos.get_row(i) - ms_pos).get(0, 1), to_vec(mpc_delay), to_vec(mpc_amplitude.real()), to_vec(mpc_amplitude.imag()));
		
		// Contains MPC angle to BS, MPC angle to MS, MPC delay, MPC complex amplitude forall MPC of the given cluster.
		channel_matrix.append_row(channel);
	}
	
	mpc_bs_pos = c.get_cluster_DMC().get_MPC_pos(BS);
	mpc_ms_pos = c.get_cluster_DMC().get_MPC_pos(MS);
	
	/*
	// Dense/Diffuse Multipaths
	// Missing within the matlab implementation and therefore commented out for now.
	for(int i = 0; i < num_dmc; i++){
		// Exact Delay for this multipath component.
		//double dmc_delay = (dist(mpc_bs_pos.get_row(i), bs_pos) + dist(mpc_ms_pos.get_row(i), ms_pos)) / speedOfLightVac + c.get_cluster_link_delay() + c.get_cluster_DMC().get_DMC_delay_add()(i);
		double dmc_delay = (dist(mpc_bs_pos.get_row(i), bs_pos) + dist(mpc_ms_pos.get_row(i), ms_pos)) / speedOfLightVac + c.get_cluster_link_delay();
		
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
 * Computes the 2x2 Polarization Matrix according to Formula 3.43 and 3.44
 * The Phase is uniformly distributed between 0 and 2 PI.
 * Currently Polarization is not used.
 * TODO: get parameter from ini files
 */
cmat Cost2100Channel::polarization_matrix(){
	cmat pol_mat(2,2);
	
	double mu_xpdv = 0;
	double sigma_xpdv = 0;
	double mu_xpdh = 0;
	double sigma_xpdh = 0;
	double mu_cpr = 0;
	double sigma_cpr = 0;
	
	double XPDV = exp(mu_xpdv + randn() * sigma_xpdv);
	double XPDH = exp(mu_xpdh + randn() * sigma_xpdh);
	double CPR = exp(mu_cpr + randn() * sigma_cpr);
	
	vec rand_phase = zeros(4);
	rand_phase.set(0,uniform(0,2 * PI));
	rand_phase.set(1,uniform(0,2 * PI));
	rand_phase.set(2,uniform(0,2 * PI));
	rand_phase.set(3,uniform(0,2 * PI));
	
	pol_mat.set(0,0, exp(complex<double>(0,rand_phase(0))) );
	pol_mat.set(0,1, (1.0 / sqrt(XPDH * CPR)) * exp(complex<double>(0,rand_phase(2))) );
	pol_mat.set(1,0, (1.0 / sqrt(XPDV)) * exp(complex<double>(0,rand_phase(3))) );
	pol_mat.set(1,1, (1.0 / sqrt(CPR)) * exp(complex<double>(0,rand_phase(4))) );
	
	// Normalize Matrix usuing Frobenius/Euklidian norm (Formula 3.44)
	return (pol_mat / norm(pol_mat,"fro"));
}

/*
 * Computes the VR Gain according to formula 3.35
 */
double Cost2100Channel::calc_VRGain(double lambda, double TR_Radius, double VR_Radius, vec VR_Center, vec MS_Pos){
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
void Cost2100Channel::calc_MIMO_Channel_Matrix(	int num_rx_antennae, int num_tx_antennae,
								mat const &rx_positions, mat const &tx_positions,
								double azimuthOfArrival, double elevationOfArrival,
								double azimuthOfDeparture, double elevationOfDeparture,	
								double lambda, cmat &result){
	
	//std::cout << std::endl << "Azimuth of Arrival: " << azimuthOfArrival << std::endl;
	//std::cout << "Elevation of Arrival: " << elevationOfArrival << std::endl;
	//std::cout << "Azimuth of Departure: " << azimuthOfDeparture << std::endl;
	//std::cout << "Elevation of Departure: " << elevationOfDeparture << std::endl << std::endl;
									
	// Compute unit direction of wave vector for arrival
	vec waveVectorDirection(3);
	waveVectorDirection.set(0,sin(azimuthOfArrival) * cos(elevationOfArrival));
	waveVectorDirection.set(1,sin(azimuthOfArrival) * sin(elevationOfArrival));
	waveVectorDirection.set(2,cos(azimuthOfArrival));
	
	double temp = ((2 * M_PI) / lambda);
	
	// Compute the wave vector assuming a plane wave
	vec k_arrival =  temp * waveVectorDirection;
	//std::cout << "k_arrival: " << k_arrival << std::endl;
	
	// Compute unit direction of wave vector for departure
	vec waveVectorDirection2(3);
	waveVectorDirection2.set(0,sin(azimuthOfDeparture) * cos(elevationOfDeparture));
	waveVectorDirection2.set(1,sin(azimuthOfDeparture) * sin(elevationOfDeparture));
	waveVectorDirection2.set(2,cos(azimuthOfDeparture));
	
	// Compute TX Steering vector
	vec k_departure =  temp * waveVectorDirection2;
	//std::cout << "k_departure: " << k_departure << std::endl;
	
	// Compute Receiver (MS) Steering vector
	cvec rx_steering(num_rx_antennae);
	for(int i = 0; i < rx_positions.rows();i++){
		std::complex<double> entry(0,dot(k_arrival,rx_positions.get_row(i)));
		rx_steering.set(i,exp(-1.0 * entry));
	}
	
	// Compute Transmitter (BS) Steering vector
	cvec tx_steering(num_tx_antennae);
	for(int i = 0; i < tx_positions.rows();i++){
		std::complex<double> entry(0,dot(k_departure,tx_positions.get_row(i)));
		tx_steering.set(i,exp(-1.0 * entry));
	}
	//std::cout << "RX Steering: " << rx_steering << std::endl;
	//std::cout << "TX Steering: " << tx_steering << std::endl;
	result = ((cmat) rx_steering) * tx_steering.transpose();
}

/*
 * Computes the actual SINR for a given Ressource Block
 * according to the Cost 2100 Channel Model
 */
double Cost2100Channel::calcSINR(int RB, vector<double> &power, vector<Position> &pos, vector<int> &bsId_, bool up, int msId){
    vec SenderGain = power[0] * H(localcluster_bs, localclusters_ms, remoteCluster, VR, 5, 0, bsId, pos[0], pos[1], LOS_VR, RB);
     
    //cout << "Test, BSID: " << bsId << endl;

    //calculate the sum of the inverse of all interferer pathloss at the target
    int numberOfInterferer = pos.size() - 2;
    //cout << "num Interferer: " << numberOfInterferer << endl;
    vec sumOfInterfererGain = zeros(SenderGain.size());
    //std::cout << "Number of interferers: " << numberOfInterferer << std::endl;
    //std::cout << "Number of interferer Powers: " << power.size() << std::endl;
    //std::cout << "Number of interferer Positions: " << pos.size() << std::endl;
    for(int i = 0; i < numberOfInterferer; i++)  {
        sumOfInterfererGain += power[i + 2] * H(localcluster_bs, localclusters_ms, remoteCluster, VR, 5, 0, bsId_[i + 2], pos[i + 2], pos[0], LOS_VR, RB);
    }
    // Assume Room Temperature and one LTE RB Bandwidth for adding Noise
    sumOfInterfererGain = sumOfInterfererGain + getTermalNoise(300, 180000);
    
    //std::cout << "Number of Interferers: " << numberOfInterferer << std::endl;
    //std::cout << "Signal Gain: " << SenderGain << std::endl;
    //std::cout << "Interferer Gain: " << sumOfInterfererGain << std::endl;

	// Elementwise Division (Divide Signal by Interference per ressourceblock)
    vec sinr = 10 * log10(itpp::elem_div(SenderGain,sumOfInterfererGain));
    // Return SINR in DB
    return sinr(0);
    
    int i = uniform(10000,20000);
    double result = 0;
    for(int j = 0; j < i; j++){
        result += uniform(1,2);
    }
    return result;
}

/*
 * Computes the actual SINR for a all Ressource Blocks
 * according to the Cost 2100 Channel Model
 */
vec Cost2100Channel::calcSINR(vector<double> &power, vector<Position> &pos, vector<int> &bsId_, bool up, int msId){
    
    vec SenderGain = power[0] * H(localcluster_bs, localclusters_ms, remoteCluster, VR, 5, 0, bsId, pos[0], pos[1], LOS_VR);
    
    //cout << "Test, BSID: " << bsId << endl;

    //calculate the sum of the inverse of all interferer pathloss at the target
    int numberOfInterferer = pos.size() - 2;
    //cout << "num Interferer: " << numberOfInterferer << endl;
    vec sumOfInterfererGain = zeros(SenderGain.size());
    //std::cout << "Number of interferers: " << numberOfInterferer << std::endl;
    //std::cout << "Number of interferer Powers: " << power.size() << std::endl;
    //std::cout << "Number of interferer Positions: " << pos.size() << std::endl;
    for(int i = 0; i < numberOfInterferer; i++)  {
        sumOfInterfererGain += power[i + 2] * H(localcluster_bs, localclusters_ms, remoteCluster, VR, 5, 0, bsId_[i + 2], pos[i + 2], pos[0], LOS_VR);
    }
    // Assume Room Temperature and one LTE RB Bandwidth for adding Noise
    sumOfInterfererGain = sumOfInterfererGain + getTermalNoise(300, 180000);
    
    //std::cout << "Number of Interferers: " << numberOfInterferer << std::endl;
    //std::cout << "Signal Gain: " << SenderGain << std::endl;
    //std::cout << "Interferer Gain: " << sumOfInterfererGain << std::endl;

	// Elementwise Division (Divide Signal by Interference per ressourceblock)
    vec sinr = 10 * log10(itpp::elem_div(SenderGain,sumOfInterfererGain));
    // Return SINR in DB
    //cout << "length: " << sinr.size() << endl;
    return sinr;
    
    int i = uniform(10000,20000);
    vec result = zeros(100);
    for(int j = 0; j < i; j++){
        result += uniform(1,2) * randn(100);
    }
    return result;
}

void Cost2100Channel::updateChannel(Position** msPos){
	//Position[numberOfMobileStations] = msPos[bsId] // datastr instead of BS ID???
	//TODO:
}

Cost2100Channel::~Cost2100Channel(){
}
