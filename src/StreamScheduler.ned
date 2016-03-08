simple StreamScheduler {
	parameters:
		double initOffset @unit(s);
		int resourceBlocks;
		int numberOfMobileStations;
		double streamSchedPeriod @unit(s);
		double tti @unit(s);
	gates:
		input fromBs;
		input fromMs[numberOfMobileStations];
		input fromRB[resourceBlocks];
		output toBs;
		output toMs[numberOfMobileStations];
		output toRB[resourceBlocks];
}
