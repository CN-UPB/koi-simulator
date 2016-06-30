// The transmission layer of the MobileStation.
//
// This MobileStation module serves as the transmission layer, sending and 
// receiving packets. It does not do any computations, but currently only 
// forwards messages. Received messages are always forwarded to the MsMac.
// Packets intended for transmission are send either to the BsChannel of the 
// local BaseStation for uplink communication or to the MsChannel of the 
// destination MobileStation for a D2D stream.
simple MsPhy  {
	parameters:
		// Number of MobileStations in the current cell
		int numberOfMobileStations;
		@display("i=block/ifcard");
    	
    gates:
	// For the receival of simulation traffic packets from local 
	// BaseStation and MobileStations.
    	input fromChannel;
	// For uplink transmission of simulation traffic packets to the 
	// local BaseStation, connected to BsChannel.
    	output toChannel;
	// For information messages to the MsChannel module
    	output toMsChannel;
	// Communication with the MsMac module
    	input fromMac;
	// Communication with the MsMac module
    	output toMac;
	// Connection to channels of MS in same cell, for transmission of 
	// simulation traffic packets.
	output toMs[numberOfMobileStations];
}
