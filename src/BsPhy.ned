simple BsPhy  {
	parameters:
		int numberOfMobileStations;
    	@display("i=block/ifcard");
    	
    gates:
    	input fromChannel[];
    	output toChannel[];
    	
    	input fromMac;
    	output toMac;
}
