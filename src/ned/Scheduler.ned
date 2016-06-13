module Scheduler{
	parameters:
		double initOffset @unit(s);
		int upResourceBlocks;
		int downResourceBlocks;
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
			numberOfMobileStations = numberOfMobileStations;
			streamSchedPeriod = streamSchedPeriod;
			tti = tti;
		}
		rbSchedulersUp[upResourceBlocks]: <rbSchedType> like RBScheduler{
			rbNumber = index;
		}
		rbSchedulersDown[downResourceBlocks]: <rbSchedType> like RBScheduler{
			rbNumber = index;
		}
	connections:
		for i=0..upResourceBlocks-1 {
			rbSchedulersUp[i].scheduler <--> streamSched.upRB[i];
		}

		for i=0..downResourceBlocks-1 {
			rbSchedulersDown[i].scheduler <--> streamSched.downRB[i];
		}

		fromBs --> streamSched.fromBs;
		streamSched.toBs --> toBs;

		for i=0..numberOfMobileStations-1 {
			fromMs[i] --> streamSched.fromMs[i];
			streamSched.toMs[i] --> toMs[i];
		}
}
