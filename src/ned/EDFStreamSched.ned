// Earliest Deadline First scheduling for the StreamScheduler
//
// This stream scheduler implements Earliest Deadline First scheduling for 
// stream to resource block assignments. It tries to utilize the smallest 
// possible number of resource blocks while ensuring that all streams can 
// meet their respective deadlines. 
//
// For details, see S. Auroux, D. Parruca and H. Karl, "Joint real-time 
// scheduling and interference coordination for wireless factory automation."
simple EDFStreamSched like StreamScheduler{
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
		// The inter-scheduling period is computed dynamically 
		// in this scheduler
		double streamSchedPeriod = -1s @unit(s);
		// Length of single transmission time interval
		double tti @unit(s);
		// A small offset used control message order in cases where ordering 
		// cannot be guaranteed by generation order.
		double epsilon @unit(s);
		int bsId;
		// Do Debug output?
		bool debug = default(false);
		// Should UP resource blocks be scheduled statically?
		bool upStatic = default(false);
		// Should DOWN resource blocks be scheduled statically?
		bool downStatic = default(false);
		// Length of the static schedule, in TTI
		int staticSchedLength = default(1);
		// UP Resource blocks assigned to the local cell
		string assignedUpRB = default("");
		// DOWN Resource blocks assigned to the local cell
		string assignedDownRB = default("");
		@class("EDFStreamScheduler");
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
