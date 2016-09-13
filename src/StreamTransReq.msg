// This message is send to the Scheduler by Mobile stations and Base stations 
// to request the use of a Ressource Block for sending in the current TTI.

cplusplus {{
#include "MessageTypes.h"
#include "KoiData_m.h"
#include <list>
using QueuePtr = std::list<KoiData*>*;
}}

class noncobject QueuePtr;
enum MessageType;
enum MessageDirection;

message StreamTransReq  {
	
	int src; // The MS from which the packet was originally send
	int dest; // The MS to which the Packet should be send
	simtime_t period; // The (average) time between two packets in this stream
	QueuePtr packets; // The current packet queue for this stream   
	int messageDirection @enum(MessageDirection); // The direction of the message, up, down or D2D
	unsigned long streamId; // The unique ID identifying the stream this schedule concerns
        int requestOrigin; // The id of the MS where the request originated, or -1 if it originated at the BS
	kind = MessageType::streamTransReq;
}
