// The transmission layer of the BaseStation
//
// This BaseStation submodule forwards messages and is the main 
// receiving/transmitting module for the BaseStation. It receives simulation 
// traffic packets successfully transmitted from local MobileStations from the 
// BsChannel instance and forwards packets destined for the local MobileStations 
// to their respective MsChannel modules.
simple BsPhy  {
	parameters:
		// Number of MobileStations in the current cell
		int numberOfMobileStations;
		@display("i=block/ifcard");
    	
    gates:
	// For simulation traffic KoiData packets successfully received by 
	// the BsChannel. The gate id corresponds to MobileStation id of the 
	// originating MobileStation.
    	input fromChannel[];
	// For simulation traffic KoiData packets destined for local 
	// MobleStations. The gate id corresponds to the destination 
	// MobileStation id.
    	output toChannel[];
	// For messages from the BsMac.
    	input fromMac;
	// For messages to the BsMac.
    	output toMac;
}
