// This message is send to the Scheduler by Mobile stations and Base stations 
// to request the use of a Ressource Block for sending in the current TTI.

cplusplus {{
#include <omnetpp.h>
#include "MessageTypes.h"
using QueuePtr = const cQueue*;
}}

class noncobject QueuePtr;
enum MessageType;

message StreamTransReq  {
	
	int src; // The MS from which the packet was originally send
	int dest; // The MS to which the Packet should be send
	simtime_t period; // The (average) time between two packets in this stream
	QueuePtr packets; // The current packet queue for this stream   
	bool bs; // Is this request originating from the Base Station?
	kind = MessageType::streamTransReq;
}
