simple MsPhy  {
	parameters:
		int numberOfMobileStations;
    	@display("i=block/ifcard");
    	
    gates:
    	input fromChannel;
    	output toChannel; //for data packets
    	output toMsChannel; //for own position informations
    	
    	input fromMac;
    	output toMac;
	// Connection to channels of MS in same cell
	output toMs[numberOfMobileStations];
}
