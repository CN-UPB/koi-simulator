// This message is returned from the Ressource Block schedulers to the 
// MS and BS, answering StreamTransReqs. It tells the MS/BS which packets for 
// a given stream should be send out in the current TTI.

cplusplus {{
#include "MessageTypes.h"
#include "KoiData_m.h"
#include <list>
using PacketIter = std::list<KoiData*>::const_iterator;
}}

class noncobject PacketIter;
enum MessageType;
enum MessageDirection;

message StreamTransSched  {
	int src; // The MS from which the packet was originally send
	int dest; // The MS to which the Packet should be send
	int rb; // The resource block this stream should use
	PacketIter packetPos; // Iterator pointing to the packet's position in the stream queue
	int messageDirection @enum(MessageDirection); // The transmission direction of this schedule
	unsigned long streamId; // The unique ID identifying the stream this schedule concerns
	kind = MessageType::streamSched;
}
