simple RoundRobinRBSched like RBScheduler{
	parameters:
		int rbNumber; // This is the ID of the RB this scheduler handles
		@class("RBScheduler");
	gates:
		input fromScheduler;
		output toScheduler;
}
