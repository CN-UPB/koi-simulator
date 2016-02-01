simple MsPhy  {
	parameters:
    	@display("i=block/ifcard");
    	
    gates:
    	input fromChannel;
    	output toChannel; //for data packets
    	output toMsChannel; //for own position informations
    	
    	input fromMac;
    	output toMac;
}
