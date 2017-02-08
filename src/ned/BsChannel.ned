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
		// Should interference be considered in SINR calculations?
		bool considerInterference = default(true); 
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
		// The transmission time interval
		double tti @unit(s);
		// A small offset used control message order in cases where ordering 
		// cannot be guaranteed by generation order.
		double epsilon @unit(s);
		// Offset from the simulation start reserved for initialization messages.
		// No traffic is generated/send before it.
		double initOffset @unit(s);
		// Channel bandwidth
		double chnBandwidth = default(180000Hz) @unit(Hz);
		// Cross correlation parameters for METIS
		double cross_a = default(0.0);
		double cross_b = default(0.0);
		double cross_c = default(0.0);
		double cross_d = default(0.0);
		double cross_e = default(0.0);
		double cross_f = default(0.0);
		double cross_m = default(0.0);
		double cross_n = default(0.0);
		double cross_q = default(0.0);
		double cross_r = default(0.0);
		double cross_s = default(0.0);
		double cross_t = default(0.0);
		double cross_u = default(0.0);
		double cross_v = default(0.0);
		double cross_w = default(0.0);
		double cross_a_LOS = default(0.0);
		double cross_b_LOS = default(0.0);
		double cross_c_LOS = default(0.0);
		double cross_d_LOS = default(0.0);
		double cross_e_LOS = default(0.0);
		double cross_f_LOS = default(0.0);
		double cross_g_LOS = default(0.0);
		double cross_h_LOS = default(0.0);
		double cross_k_LOS = default(0.0);
		double cross_l_LOS = default(0.0);
		double cross_m_LOS = default(0.0);
		double cross_n_LOS = default(0.0);
		double cross_o_LOS = default(0.0);
		double cross_p_LOS = default(0.0);
		double cross_q_LOS = default(0.0);
		double cross_r_LOS = default(0.0);
		double cross_s_LOS = default(0.0);
		double cross_t_LOS = default(0.0);
		double cross_u_LOS = default(0.0);
		double cross_v_LOS = default(0.0);
		double cross_w_LOS = default(0.0);
		double Decorr_LOS_DS = default(0.0);
		double Decorr_LOS_ASD = default(0.0);
		double Decorr_LOS_ASA = default(0.0);
		double Decorr_LOS_ZSD = default(0.0);
		double Decorr_LOS_ZSA = default(0.0);
		double Decorr_LOS_SF = default(0.0);
		double Decorr_LOS_K = default(0.0);
		double Decorr_NLOS_DS = default(0.0);
		double Decorr_NLOS_ASD = default(0.0);
		double Decorr_NLOS_ASA = default(0.0);
		double Decorr_NLOS_ZSD = default(0.0);
		double Decorr_NLOS_ZSA = default(0.0);
		double Decorr_NLOS_SF = default(0.0);
		double DS_mu_NLOS = default(0.0);
		double DS_eps_NLOS = default(0.0);
		double AoD_mu_NLOS = default(0.0);
		double AoD_eps_NLOS = default(0.0);
		double AoA_mu_NLOS = default(0.0);
		double AoA_eps_NLOS = default(0.0);
		double ZoA_mu_NLOS = default(0.0);
		double ZoA_eps_NLOS = default(0.0);
		double SF_sigma_NLOS = default(0.0);
		double DS_mu_LOS = default(0.0);
		double DS_eps_LOS = default(0.0);
		double AoD_mu_LOS = default(0.0);
		double AoD_eps_LOS = default(0.0);
		double AoA_mu_LOS = default(0.0);
		double AoA_eps_LOS = default(0.0);
		double ZoA_mu_LOS = default(0.0);
		double ZoA_eps_LOS = default(0.0);
		double SF_sigma_LOS = default(0.0);
		double K_mu = default(0.0);
		double K_sigma = default(0.0);
		int NumberOfClusters_LOS = default(0);
		int NumberOfClusters_NLOS = default(0);
		double DelayScaling_LOS = default(0.0);
		double DelayScaling_NLOS = default(0.0);
		double PerClusterShadowing_LOS = default(0.0);
		double PerClusterShadowing_NLOS = default(0.0);
		int NumberOfRays_NLOS = default(0);
		int NumberOfRays_LOS = default(0);
		double Cluster_ASD_NLOS = default(0.0);
		double Cluster_ASA_NLOS = default(0.0);
		double Cluster_ZSA_NLOS = default(0.0);
		double Cluster_ASD_LOS = default(0.0);
		double Cluster_ASA_LOS = default(0.0);
		double Cluster_ZSA_LOS = default(0.0);
		double XPR_Mean_NLOS = default(0.0);
		double XPR_Std_NLOS = default(0.0);
		double XPR_Mean_LOS = default(0.0);
		double XPR_Std_LOS = default(0.0);
		double CarrierFrequency = default(5.2e9Hz) @unit(Hz);
		// Number of antennas on BaseStation
		int NumBsAntenna = default(8);
		// Number of antennas on MobileStation
		int NumMsAntenna = default(2);
		double BsHeight = default(0.0);
		double OutdoorHeightUE = default(0.0);
		// K mean for Ricean fading models
		double kMean = default(1.0);
		// K standard deviation for Ricean fading models
		double kSigma = default(1.0);
		// The 'mean' parameter for fast fading channel models
		double expMean = default(1.0);
		// Pathloss exponent
		double plExp = default(3.0);
		// Reference distance for pathloss models
		double d0 = default(1m) @unit(m);
		// Pathloss at the reference distance for pathloss models
		double pl0 = default(1);
		// The channel model to use. Possible values:
		// 0: METIS channel model (default)
		// 1: Exponential random channel model
		// 2: Factory channel model, Tanghe et al.
		int channelModel = default(0);
		// Table with min SINR for Modulation and Coding Schemes
		string MCSTable;
		// Bandwidth available to each resource block
		double bandwidthPerRB @unit(Hz);
		// Gain of BS antennas in dB
		double bsGain = default(1.0dB) @unit(dB);
		// Gain of MS antennas in dB
		double msGain = default(1.0dB) @unit(dB);
		// Transmission power in Watts
		double transmissionPower = default(1.0W) @unit(W);
		// Standard deviation of the normal distribution for Shadowing
		double shadowSigma = default(1.0);
		// Do Debug output?
		bool debug = default(false);
		// Should D2D values be computed?
		bool d2dActive = default(true);
		// Velocity of the mobile stations in m/s
		double msVel;
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
