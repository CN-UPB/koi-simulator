// This message is returned from the Ressource Block schedulers to the 
// MS and BS, answering StreamTransReqs. It tells the MS/BS which packets for 
// a given stream should be send out in the current TTI.

cplusplus {{
#include "MessageTypes.h"
}}

enum MessageType;
enum MessageDirection;

message StreamTransSched  {
	int src; // The MS from which the packet was originally send
	int dest; // The MS to which the Packet should be send
	int rb; // The resource block this stream should use
	int packetIndex; // The index of the packet to be send in the packet queue
	int messageDirection @enum(MessageDirection); // The transmission direction of this schedule
	unsigned long streamId; // The unique ID identifying the stream this schedule concerns
	kind = MessageType::streamSched;
}
