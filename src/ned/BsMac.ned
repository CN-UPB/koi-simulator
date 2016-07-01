simple BsMac  {
    parameters:
    	int maxNumberOfNeighbours;
    	int numberOfMobileStations;
    	int bsId;
        int currentChannel;
	int downResourceBlocks;
    	double xPos;
    	double yPos;
    	double epsilon @unit(s);
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
