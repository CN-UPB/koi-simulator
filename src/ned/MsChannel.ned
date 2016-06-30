simple MsChannel  {
	parameters:
		int maxNumberOfNeighbours;
		int downResourceBlocks;
		int upResourceBlocks;
		int bsId;
		int msId;
		double epsilon @unit(s);
		double tti @unit(s);
		double packetLoss;
		@display("i=block/rxtx");
		string eesm_beta;
		string bler_table;
	gates:
		input fromBs;
		output toPhy;
		input fromPhy; //just for own position exchange

		// For receival of D2D messages from local MS
		input fromMs[];
}
