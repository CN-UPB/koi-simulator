// Scheduler for simulation traffic packets
//
// The Scheduler module is responsible for scheduling packet transmission for 
// a single LteCell. It has to ensure that no two transmissions occur on the 
// same resource block in the same TTI in the local cell. 
// To ensure this, it receives StreamInfo messages from all MobileStations in 
// it's cell, carrying information about the period between packets of that 
// stream and the communication partners as well as the transmission mode. 
//
// <h3>Submodules</h3>
// To produce schedules, the Scheduler module uses two submodules, the 
// StreamScheduler, which assigns resource blocks on a by stream basis,
// and several RBScheduler modules, which decide which packet should be 
// transmitted in each TTI for their resource block. The RBScheduler array 
// is divided into two groupes, one for each transmission direction, uplink
// and downlink.
//
// To allow easy replacement of scheduling algorithms, both submodules 
// have the <tt>moduleinterface</tt> type. Which implementation of those 
// interfaces is to be used for a given simulation is set via the 
// <tt>streamSchedType</tt> and <tt>rbSchedType</tt> parameters.
//
// Currently, there is only one implementation for each interface. The 
// RoundRobinStreamSched interface for the StreamScheduler and the 
// RoundRobinRBSched for the RBScheduler.
//
// A schedule for a single TTI works as follows:
// -# All MobileStations and BaseStations send their StreamTransReq messages, 
// for each of the streams which currently have packets awaiting transmission
// at the beginning of the TTI.
// -# After <tt>tti+epsilon</tt> time, all StreamTransReq messages are 
// forwarded to the RBScheduler currently assigned to handling their stream.
// -# Every RBScheduler decides on which packet from which StreamTransReq for 
// their streams should be send on it's resource block and notes that packet 
// in it's StreamTransSched message. That message is then forwarded to the 
// originating MobileStation.
module Scheduler
{
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
	// The module implementation to use for the StreamScheduler.
        string streamSchedType = default("RoundRobinStreamSched");
	// The module implementation to use for the RBScheduler.
        string rbSchedType = default("RoundRobinRBSched");
				int bsId;
        @display("i=block/classifier;bgb=366,173");
    gates:
	// For Scheduling communication with the local BaseStation.
        input fromBs;
	// For Scheduling communication with the local MobileStations.
        input fromMs[numberOfMobileStations];
	// For Scheduling communication with the local BaseStation.
        output toBs;
	// For Scheduling communication with the local MobileStations.
        output toMs[numberOfMobileStations];
    submodules:
        streamSched: <streamSchedType> like StreamScheduler {
            initOffset = initOffset;
            numberOfMobileStations = numberOfMobileStations;
            tti = tti;
						bsId = bsId;
            @display("p=69,79");
        }
        rbSchedulersUp[upResourceBlocks]: <rbSchedType> like RBScheduler {
            rbNumber = index;
            @display("p=220,125");
        }
        rbSchedulersDown[downResourceBlocks]: <rbSchedType> like RBScheduler {
            rbNumber = index;
            @display("p=220,27");
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
