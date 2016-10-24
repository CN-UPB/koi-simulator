simple EDFRBSched like RBScheduler{
	parameters:
		int rbNumber; // This is the ID of the RB this scheduler handles
		// Number of subcarriers for OFDM
		int numSubcarriers = default(12);
		@class("EDFRBSched");
	gates:
		inout scheduler;
}
