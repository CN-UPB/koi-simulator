simple MsChannel  {
	parameters:
		int maxNumberOfNeighbours;
		int downResourceBlocks;
		int upResourceBlocks;
		int currentChannel;
		bool useSimpleChannelCalc;
		int bsId;
		int msId;
		double epsilon @unit(s);
		double tti @unit(s);
		int simpleChannelCalcNops;
		double packetLoss;
		@display("i=block/rxtx");
		double averageFarCluster;
		double vrRadius;
		double trRadius;
		double cellRadius;
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
		string eesm_beta;
		string bler_table;
		
	gates:
		input fromBs;
		output toPhy;
		input fromPhy; //just for own position exchange

		// For receival of D2D messages from local MS
		input fromMs[];
}
