module Scheduler{
	parameters:
		double initOffset @unit(s);
		int resourceBlocks;
		int numberOfMobileStations;
		double streamSchedPeriod @unit(s);
		double tti @unit(s);
		string streamSchedType = default("RoundRobinStreamSched");
		string rbSchedType = default("RoundRobinRBSched");
	gates:
		input fromBs;
		input fromMs[numberOfMobileStations];
		output toBs;
		output toMs[numberOfMobileStations];
	submodules:
		streamSched: <streamSchedType> like StreamScheduler{
			initOffset = initOffset;
			resourceBlocks = resourceBlocks;
			numberOfMobileStations = numberOfMobileStations;
			streamSchedPeriod = streamSchedPeriod;
			tti = tti;
		}
		rbSchedulers[resourceBlocks]: <rbSchedType> like RBScheduler{
			rbNumber = index;
		}
	connections:
		for i=0..resourceBlocks-1 {
			rbSchedulers[i].toScheduler --> streamSched.fromRB[i];
			streamSched.toRB[i] --> rbSchedulers[i].fromScheduler;
		}

		fromBs --> streamSched.fromBs;
		streamSched.toBs --> toBs;

		for i=0..numberOfMobileStations-1 {
			fromMs[i] --> streamSched.fromMs[i];
			streamSched.toMs[i] --> toMs[i];
		}
}
