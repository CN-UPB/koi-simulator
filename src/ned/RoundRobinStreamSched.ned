simple RoundRobinStreamSched like StreamScheduler{
	parameters:
		double initOffset @unit(s);
		int upResourceBlocks;
		int downResourceBlocks;
		int numberOfMobileStations;
		double streamSchedPeriod @unit(s);
		double tti @unit(s);
		@class("StreamScheduler");
	gates:
		input fromBs;
		input fromMs[numberOfMobileStations];
		output toBs;
		output toMs[numberOfMobileStations];
		inout upRB[upResourceBlocks];
		inout downRB[downResourceBlocks];
}
