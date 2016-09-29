// Communication channel for incoming packets from MobileStations
//
// The BsChannel is responsible for handling packets send to the BaseStation.
// It calculates the SINR value for each incoming packet and decides whether 
// it was received successfully or not. If the packet was received, it is then 
// forwarded to the BsPhy for further handling. Otherwise, the packet is 
// dropped.
//
// Computing the SINR values requires the BsChannel to know which other
// transmissions occured in the current TTI. To that end, all BsChannel instances 
// receive TransInfo messages from neighbouring cells, carrying information 
// about possible interference.
//
// <h3>The METIS channel</h3>
// To compute SINR values for packets, the BsChannel uses the physical 
// layer model developed by the METIS project. All documentation can be found 
// at the <a href="https://www.metis2020.com/">project homepage</a>.
//
// To keep the documentation relevant to the simulation model, those parameters 
// used only by the METIS implementation are not documented here. Please refer 
// to the <a href="https://www.metis2020.com/wp-content/uploads/METIS_D1.4_v3.pdf">
// METIS Deliverable 1.4</a> for an in depth presentation and definition of the 
// model.
//
// To make certain that there is only one instance of the METISChannel class 
// per LteCell, only the very first BsChannel instance with the index 0 
// instantiates the channel. Then, it sends PointerExchange messages 
// to all other BsChannel instances as well as all local MobileStations,
// containing a pointer to that instance. 
simple BsChannel  {
	parameters:
		// The number of neighbouring LteCells
		int maxNumberOfNeighbours;
		// The number of MobileStations in the local cell
		int numberOfMobileStations;
		// The number of resource blocks for uplink transmissions
		int upResourceBlocks;
		// The number of resource blocks for downlink transmissions
		int downResourceBlocks;
		// The id of the BaseStation. This id is unique.
		int bsId;
		// X coordinate of the BaseStation's position
		double xPos;
		// Y coordinate of the BaseStation's position
		double yPos;
		// Z coordinate of the BaseStation's position
		double zPos = default(0);
		bool useSimpleChannelCalc;
		int simpleChannelCalcNops;
		int FADING_PATHS;
		int ga_seed;
		int random_seed_id;
		int CHNrandom_seed_id;
		bool calculatePathLoss;
		bool calculateShadowing;
		bool calculateFading;
		double stationaryDoppler;
		double delay_rms;
		double tenlogk;
		double alpha;
		double packetLoss;
		// The transmission time interval
		double tti @unit(s);
		// A small offset used control message order in cases where ordering 
		// cannot be guaranteed by generation order.
		double epsilon @unit(s);
		// Offset from the simulation start reserved for initialization messages.
		// No traffic is generated/send before it.
		double initOffset @unit(s);
		double averageFarCluster;
		double vrRadius;
		double trRadius;
		// The radius of the LteCell
		double cellRadius;
		// The radius of the LteCell
		double cellRadiusMETIS;
		double commonCluster;
		double cluster_azimuth_base;
		double cluster_azimuth_range;
		double cluster_azimuth_distr;
		double cluster_elevation_base;
		double cluster_elevation_range;
		double cluster_elevation_distr;
		double cluster_distance_base;
		double cluster_distance_range;
		double cluster_distance_distr;
		double linkDelay_mean;
		double linkDelay_std;
		double delay_mean;
		double delay_std;
		double singleClusterRatio;
		double clusterPower;
		double excessDelay;
		int activeLocalClusterBS;
		int activeLocalClusterMS;
		double averageLocalCluster;
		double NumberOfMPC;
		double AoD_mean;
		double AoD_std;
		double EoD_mean;
		double EoD_std;
		double AoA_mean;
		double AoA_std;
		double EoA_mean;
		double EoA_std;
		double shadowFading_mean;
		double shadowFading_std;
		double shadowFadingAutoCorrelation;
		double factorLOS_std;
		double factorLOS_mean;
		double LOS_cutoff_distance;
		double VR_radius_los;
		double TR_radius_los;
		int tx_num;
		int rx_num;
		double startfrequency;
		double endfrequency;
		double subcarrier_bandwidth;
		string crossCorrelation;
		// Table of EESM values for CQI computation
		string eesm_beta;
		// Bit Level Error Rate table for packet loss computations
		string bler_table;
		double cross_a;
		double cross_b;
		double cross_c;
		double cross_d;
		double cross_e;
		double cross_f;
		double cross_m;
		double cross_n;
		double cross_q;
		double cross_r;
		double cross_s;
		double cross_t;
		double cross_u;
		double cross_v;
		double cross_w;
		double cross_a_LOS;
		double cross_b_LOS;
		double cross_c_LOS;
		double cross_d_LOS;
		double cross_e_LOS;
		double cross_f_LOS;
		double cross_g_LOS;
		double cross_h_LOS;
		double cross_k_LOS;
		double cross_l_LOS;
		double cross_m_LOS;
		double cross_n_LOS;
		double cross_o_LOS;
		double cross_p_LOS;
		double cross_q_LOS;
		double cross_r_LOS;
		double cross_s_LOS;
		double cross_t_LOS;
		double cross_u_LOS;
		double cross_v_LOS;
		double cross_w_LOS;
		double Decorr_LOS_DS;
		double Decorr_LOS_ASD;
		double Decorr_LOS_ASA;
		double Decorr_LOS_ZSD;
		double Decorr_LOS_ZSA;
		double Decorr_LOS_SF;
		double Decorr_LOS_K;
		double Decorr_NLOS_DS;
		double Decorr_NLOS_ASD;
		double Decorr_NLOS_ASA;
		double Decorr_NLOS_ZSD;
		double Decorr_NLOS_ZSA;
		double Decorr_NLOS_SF;
		double DS_mu_NLOS;
		double DS_eps_NLOS;
		double AoD_mu_NLOS;
		double AoD_eps_NLOS;
		double AoA_mu_NLOS;
		double AoA_eps_NLOS;
		double ZoA_mu_NLOS;
		double ZoA_eps_NLOS;
		double SF_sigma_NLOS;
		double DS_mu_LOS;
		double DS_eps_LOS;
		double AoD_mu_LOS;
		double AoD_eps_LOS;
		double AoA_mu_LOS;
		double AoA_eps_LOS;
		double ZoA_mu_LOS;
		double ZoA_eps_LOS;
		double SF_sigma_LOS;
		double K_mu;
		double K_sigma;
		int NumberOfClusters_LOS;
		int NumberOfClusters_NLOS;
		double DelayScaling_LOS;
		double DelayScaling_NLOS;
		double PerClusterShadowing_LOS;
		double PerClusterShadowing_NLOS;
		int NumberOfRays_NLOS;
		int NumberOfRays_LOS;
		double Cluster_ASD_NLOS;
		double Cluster_ASA_NLOS;
		double Cluster_ZSA_NLOS;
		double Cluster_ASD_LOS;
		double Cluster_ASA_LOS;
		double Cluster_ZSA_LOS;
		double XPR_Mean_NLOS;
		double XPR_Std_NLOS;
		double XPR_Mean_LOS;
		double XPR_Std_LOS;
		double positionResendInterval;
		double CarrierFrequency;
		// Number of antennas on BaseStation
		int NumBsAntenna;
		// Number of antennas on MobileStation
		int NumMsAntenna;
		double BsHeight;
		double OutdoorHeightUE;
		double Velocity;
		// The 'mean' parameter for the exponential channel model
		double expMean = default(1.0);
		// The channel model to use. Possible values:
		// 0: METIS channel model (default)
		// 1: Exponential random channel model
		int channelModel = default(0);
		@display("i=block/rxtx");
		
	gates:
		// For simulation traffic KoiData packets from local 
		// MobileStations.
		input fromMs;
		// For information messages from BsMac.
		input fromBsMac;
		// For forwarding successfully received packets to the 
		// BsMac.
		output toPhy;
}
