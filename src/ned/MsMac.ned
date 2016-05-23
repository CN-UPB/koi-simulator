simple MsMac  {
    parameters:
    	int bsId;
    	int msId;
    	double initBsXPos;
    	double initBsYPos;
	double initMsXPos = default(0.0);
	double initMsYPos = default(0.0);
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
	double transmissionPower = default(1.0);
	// The positioning algorithm for Mobile Stations inside their cell
	// 0 - Uniformely random anywhere in the cell
	// 1 - Use the positions given as OMNet parameters (initMsXPos,initMsXPos)
	// 2 - Use initial placement by sector
	// 3 - Use linear placement simulating a straight road
	int positioning = default(0);
    	
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
