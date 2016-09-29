simple EDFRBSched like RBScheduler{
	parameters:
		int rbNumber; // This is the ID of the RB this scheduler handles
		@class("EDFRBSched");
	gates:
		inout scheduler;
}
