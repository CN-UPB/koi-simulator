// The RBScheduler is responsible for scheduling all packets handled by 
// a single Ressource block, as assigned by the StreamScheduler.

moduleinterface RBScheduler {
	parameters:
		int rbNumber; // This is the ID of the RB this scheduler handles
		@display("i=block/subqueue");
	gates:
		inout scheduler;
}
