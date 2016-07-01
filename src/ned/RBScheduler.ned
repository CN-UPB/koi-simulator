// Interface for resource block scheduler module
//
// The RBScheduler is responsible for scheduling all packets handled by 
// a single Ressource block, as assigned by the StreamScheduler. 
// In the StreamScheduler, each RBScheduler is assigned a number of streams 
// to handle. Each TTI, the RBScheduler receives a TransReqList for each 
// stream it has been assigned to handle. From these lists, the RBScheduler 
// chooses exactle one packet from exactly one stream, which might then be 
// transmitted on the resource block handled by this scheduler.
//
// After a packet has been choosen, the RBScheduler notes it in a 
// StreamTransSched message, together with it's resource block id, and forwards 
// the schedule to the StreamScheduler, which then sends it to the originator 
// of the corresponding StreamTransReq. 
//
// This is a module interface, intended to be extended by implementations for 
// different packet scheduling algorithms.
moduleinterface RBScheduler {
	parameters:
		// The id of the resource block this scheduler is responsible for.
		int rbNumber;
		@display("i=block/subqueue");
	gates:
		// For communication with the StreamScheduler.
		inout scheduler;
}
