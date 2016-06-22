module MobileStation
{
    parameters:
	int bsId;
	int msId;
	int currentChannel;
	double initBsXPos;
	double initBsYPos;
	int initQuadrant;
	double initPosAlpha;
	double initPosBeta;
	double initPosGamma;
	double radius;

	@display("bgb=200,300;i=device/cellphone");

    gates:
    	//packet exchange with the BS through DataChannel
        input fromBs;
        output toBs;

	//direct schedule exchange with BS
        input fromBsMac;
        output toBsMac;

	// Access to the Stream scheduling module
	input fromScheduler;
	output toScheduler;

    	// Links to MS in same cell for D2D transmissions
    	input fromMs[];
    	output toMs[];

    submodules:
        app: TrafficGen  {
        	parameters:
        		bsId = bsId;
        		msId = msId;
        		@display("p=50,50");
       	}
        mac: MsMac  {
        	parameters:
        		bsId = bsId;
        		msId = msId;
                currentChannel = currentChannel;
                radius = radius;
        		initBsXPos = initBsXPos;
    			initBsYPos = initBsYPos;
    			initQuadrant = initQuadrant;
    			initPosAlpha = initPosAlpha;
	    		initPosBeta = initPosBeta;
	    		initPosGamma = initPosGamma;
        		@display("p=50,120");
        }
        phy: MsPhy  {
        	parameters:
        		@display("p=50,190");
		gates:
			toMs[sizeof(fromMs)];
        }
        dataChn: MsChannel  {
        	parameters:
			msId = msId;
        		bsId = bsId;
                	currentChannel = currentChannel;
        		@display("p=50,260");
		gates:
			fromMs[sizeof(fromMs)];
        }
    connections:
        app.toMac --> mac.fromApp;
        mac.toApp --> app.fromMac;
        mac.toPhy --> phy.fromMac;
        mac.toBsMac --> toBsMac;
        fromBsMac --> mac.fromBsMac;
        phy.toMac --> mac.fromPhy;
        dataChn.toPhy --> phy.fromChannel;
        phy.toChannel --> toBs;
        phy.toMsChannel --> dataChn.fromPhy;
        fromBs --> dataChn.fromBs;
	fromScheduler --> mac.fromScheduler;
	mac.toScheduler --> toScheduler;

	// For D2D message traffic
	for i=0..sizeof(fromMs)-1  {
		phy.toMs[i] --> toMs[i];
		fromMs[i] --> dataChn.fromMs[i];
	}
}
