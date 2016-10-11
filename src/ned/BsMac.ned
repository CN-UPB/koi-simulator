// Medium Access Controller and central module for the BaseStation
//
// The BsMac submodule is the core of the BaseStation. It handles decisions 
// about simulation traffic transmissions and routes messages from and to 
// neighbouring cells.
//
// At the beginning of each simulation, the BsMac waits for all local 
// MobileStations to send their position. As soon as all messages have arrived,
// the BsMac forwards a combined message with all MobileStation positions to 
// all neighbouring cells. In addition, local MobileStation positions are 
// also forwarded to all other MobileStations.
// When the MobileStation positions from other cells arrive, they are forwarded 
// to the BsChannel modules, as they are needed for SINR computations.
// 
// Also at the start of each simulation, the BsMac forwards the BaseStation 
// position to all local MobileStations and all neighbouring BaseStations.
//
// In the area of cell local communication, the BsMac is mainly concerned with 
// forwarding non-D2D traffic packets from the source to the recipient. When
// a traffic packet arrives from the BsPhy, the following steps are taken:
// -# Unpack the DataPacketBundle and sort the KoiData packets into queues by 
// their stream id.
// -# At the start of the next <em>tti</em>, send a StreamTransReq to the local 
// cell's Scheduler module for each packet queue/stream id with packets.
// -# When the StreamTransSched messages arrive from the Scheduler,
// take the corresponding messages out of their queues and wrap them in 
// DataPacketBundle messages.
// -# For each DataPacketBundle, set the resource block to the one prescribed 
// by the StreamTransSched. 
// -# For each DataPacketBundle, also generate a TransInfo message and send 
// it to all neighbouring BaseStations. 
// -# Finally, forward the DataPacketBundles to the BsPhy module for 
// transmission.
simple BsMac  {
	parameters:
		// The number of neighbouring LteCells
		int maxNumberOfNeighbours;
		// The number of MobileStations in the local cell
		int numberOfMobileStations;
		// The id of the BaseStation. This id is unique.
		int bsId;
		int currentChannel;
		// The number of resource blocks for downlink transmissions
		int downResourceBlocks;
		// X coordinate of the BaseStation's position
		double xPos;
		// Y coordinate of the BaseStation's position
		double yPos;
		// A small offset used control message order in cases where ordering 
		// cannot be guaranteed by generation order.
		double epsilon @unit(s);
		// The transmission time interval
		double tti @unit(s);
		// Offset from the simulation start reserved for initialization messages.
		// No traffic is generated/send before it.
		double initOffset @unit(s);
		// The power setting used for transmissions.
		double transmissionPower = default(1.0);
		@display("i=block/mac");
        
    gates:
	// For information exchange with neighbouring cells. 
	// <strong>Note:</strong> Gate indices do <em>NOT</em> correspond to 
	// BaseStation ids!
        input fromCell[];
	// For information exchange with neighbouring cells. 
	// <strong>Note:</strong> Gate indices do <em>NOT</em> correspond to 
	// BaseStation ids!
        output toCell[];
        
    	// For direct information exchange with MsMacs. The gate index 
	// corresponds to the MobileStation id.
        input fromMsMac[];
    	// For direct information exchange with MsMacs. The gate index 
	// corresponds to the MobileStation id.
        output toMsMac[];
        
        // For transmission/receival of simulation traffic packets
        input fromPhy;
        // For transmission/receival of simulation traffic packets
        output toPhy;
        
        // For information messages to BsChannel
        output toBsChannel[];

    	// Communication with Scheduler for transmission schedules
	input fromScheduler;
    	// Communication with Scheduler for transmission schedules
	output toScheduler;
}
