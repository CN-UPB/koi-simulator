simple BsMac  {
    parameters:
    	int maxNumberOfNeighbours;
    	int numberOfMobileStations;
    	int bsId;
        int currentChannel;
	int downResourceBlocks;
        int SlidingWindowSize;
        double SmoothingFactor;
        bool useSinr;
        string AverageType;
        string eesm_beta;
        string bler_table;
    	double xPos;
    	double yPos;
    	bool uniform;
	double offset;
	int N_s;
	int N_c;
	double runtime @unit(s);
    	
    	double epsilon @unit(s);

    	//TTI - Transmission Time Interval
    	double tti @unit(s);
    	double initOffset @unit(s);
	double transmissionPower = default(1.0);
    
        @display("i=block/mac");
        
    gates:
    	//packets for/from the application
        input fromApp[];
        output toApp[];
        
        //for MSs position exchange with the other cells
        input fromCell[];
        output toCell[];
        
        //for schedule exchange with the MSs
        input fromMsMac[];
        output toMsMac[];
        
        //input/output to MS through DataChannel
        input fromPhy;
        output toPhy;
        
        //position exchange with the BsChannels
        output toBsChannel[];

	// Access to Streaming scheduler
	input fromScheduler;
	output toScheduler;
}
