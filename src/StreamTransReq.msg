// This message is send to the Scheduler by Mobile stations and Base stations 
// to request the use of a Ressource Block for sending in the current TTI.

cplusplus {{
#include <omnetpp.h>
#include "MessageTypes.h"
using QueuePtr = const cQueue*;
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
	kind = MessageType::streamTransReq;
}
