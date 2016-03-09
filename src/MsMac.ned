simple MsMac  {
    parameters:
    	int bsId;
    	int msId;
    	double initBsXPos;
    	double initBsYPos;
	int initQuadrant;
        int currentChannel;
    	double initOffset @unit(s);
    	double tti @unit(s);
    	int positionResendInterval;
    	double initPosAlpha;
	double initPosBeta;
	double initPosGamma;
	double epsilon @unit(s);
	double radius;
	int PosSeed;
	int downResourceBlocks;
	int packetLength = default(1500b) @unit(b);
    	
        @display("i=block/mac");
        
    gates:
        input fromApp;
        output toApp;
        
        //packet exchange with BS through DataChannel
        input fromPhy;
        output toPhy;
        
        //direct schedule exchange with BS
        input fromBsMac;
        output toBsMac;

	// Access to Streaming scheduler
	input fromScheduler;
	output toScheduler;
}
