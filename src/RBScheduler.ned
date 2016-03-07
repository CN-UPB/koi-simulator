// The RBScheduler is responsible for scheduling all packets handled by 
// a single Ressource block, as assigned by the StreamScheduler.

simple RBScheduler {
	parameters:
		int rbNumber; // This is the ID of the RB this scheduler handles
	gates:
		input fromScheduler;
		output toScheduler;
}
