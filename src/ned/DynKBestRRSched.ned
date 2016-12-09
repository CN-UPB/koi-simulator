simple DynKBestRRSched like StreamScheduler{
	parameters:
		// Offset from the simulation start reserved for initialization messages.
		// No traffic is generated/send before it.
		double initOffset @unit(s);
		// The number of resource blocks for uplink transmissions
		int upResourceBlocks;
		// The number of resource blocks for downlink transmissions
		int downResourceBlocks;
		// The number of MobileStations in the local cell
		int numberOfMobileStations;
		// The period after which assignments of streams to resource blocks 
		// are recomputed.
		double streamSchedPeriod @unit(s);
		// Length of single transmission time interval
		double tti @unit(s);
		// A small offset used control message order in cases where ordering 
		// cannot be guaranteed by generation order.
		double epsilon @unit(s);
		// Number of best RB assigned to senders each TTI for UP transmissions
		int upK;
		// Number of best RB assigned to senders each TTI for DOWN transmissions
		int downK;
		int bsId;
		// Do Debug output?
		bool debug = default(false);
		// Length of packets in bits
		int packetLength @unit(b);
		// Should UP resource blocks be scheduled statically?
		bool upStatic = default(false);
		// Should DOWN resource blocks be scheduled statically?
		bool downStatic = default(false);
		// Length of the static schedule, in TTI
		int staticSchedLength = default(1);
		@class("DynKBestRRSched");
	gates:
		// For schedule exchange with the local BaseStation.
		input fromBs;
		// For Scheduling communication with the local MobileStations.
		input fromMs[numberOfMobileStations];
		// For schedule exchange with the local BaseStation.
		output toBs;
		// For Scheduling communication with the local MobileStations.
		output toMs[numberOfMobileStations];
		// For communication with uplink resource schedulers.
		inout upRB[upResourceBlocks];
		// For communication with downlink resource schedulers.
		inout downRB[downResourceBlocks];
}
