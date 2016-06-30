// A mobile sensor/actor station
// 
// This module generates all the traffic in the model and 
// communicates with other MobileStations via the local cell's BaseStation 
// or directly via D2D. Any MobileStation is part of exactly one LteCell.
//
// <h3>Submodules</h3>
// A MobileStation's submodules fullfill the following tasks:<br>
// 
// <p><strong>app:</strong> The <em>TrafficGen</em> module is responsible for 
// generating and consuming traffic. To that end, it communicates with the 
// <em>mac</em> to send and receive packets.
//				
// <p><strong>mac:</strong> The <em>MsMac</em> is the core module of the MobileStation.
// It holds all central information, like current position. It also holds 
// all packets awaiting transmission and serves as a communication hub,
// forwarding traffic from other compound modules to the intended MobileStation 
// submodule. It uses the <em>phy</em> module to send simulated traffic to other 
// modules, either a BaseStation or a neighbouring MobileStation. It also has 
// a direct connection to the local BaseStation for information exchange
// without having to go through the <em>phy</em> module. To determine 
// which packet to send, it is connected to the local cell's Scheduler.
//
// <p><strong>phy:</strong> The <em>BsPhy</em> module is merely a relay, 
// forwarding messages to other MobileStation or BaseStation modules and 
// forwarding messages received from those via the <em>dataChn</em> to the <em>mac</em>.
//
// <p><strong>dataChn</strong> The <em>dataChn</em> module serves as the 
// receiving end of any simulated traffic packets to this mobile station. 
// Whenever a packet is received, the module uses the <tt>METISChannel</tt> 
// to compute SINR values and determine whether the it was received successfully.
// In that case, they are then forwarded to the <em>phy</em> module.
module MobileStation
{
    parameters:
    	// The ID of the MobileStation's current cell.
	int bsId;
	// The ID of this MobileStation. MobileStation IDs are unique only 
	// within any given cell. Thus, to uniquely identify a MobileStation,
	// it's <tt>bsId</tt> and it's <tt>msId</tt> are needed.
	int msId;
	int currentChannel;
	// The X position of the local BaseStation.
	double initBsXPos;
	// The Y position of the local BaseStation.
	double initBsYPos;
	// The quadrant for initial positioning by sector.
	int initQuadrant;
	// The alpha angle for initial positioning by sector.
	double initPosAlpha;
	// The beta angle for initial positioning by sector.
	double initPosBeta;
	// The gamma angle for initial positioning by sector.
	double initPosGamma;
	// The radius of the local cell
	double radius;

	@display("bgb=200,300;i=device/cellphone");

    gates:
    	// Incoming DataPacketBundle from the local BaseStation. This is the 
	// gate to be used for simulated traffic, as opposed to purely 
	// implementation specific control/information messages.
        input fromBs;
    	// Outgoing DataPacketBundle to the local BaseStation. This is the 
	// gate to be used for simulated traffic, as opposed to purely 
	// implementation specific control/information messages.
        output toBs;

	// Exchange for control messages with the local BaseStation. This 
	// gate is not intended for simulated traffic.
        input fromBsMac;
	// Exchange for control messages with the local BaseStation. This 
	// gate is not intended for simulated traffic.
        output toBsMac;

	// Access to the Stream scheduling module.
	input fromScheduler;
	// Access to the Stream scheduling module.
	output toScheduler;

    	// Links to MS in same cell for D2D transmissions. This gate is 
	// intended for simulated traffic between two mobile stations.
    	input fromMs[];
    	// Links to MS in same cell for D2D transmissions. This gate is 
	// intended for simulated traffic between two mobile stations.
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
