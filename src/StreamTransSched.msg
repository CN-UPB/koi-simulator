// This message is returned from the Ressource Block schedulers to the 
// MS and BS, answering StreamTransReqs. It tells the MS/BS which packets for 
// a given stream should be send out in the current TTI.

cplusplus {{
#include "MessageTypes.h"
}}

enum MessageType;

message StreamTransSched  {
	int src; // The MS from which the packet was originally send
	int dest; // The MS to which the Packet should be send
	int rb; // The resource block this stream should use
	int packetIndex; // The index of the packet to be send in the packet queue
	bool bs; // Is this schedule intended for the BS or a MS?
	kind = MessageType::streamSched;
}
