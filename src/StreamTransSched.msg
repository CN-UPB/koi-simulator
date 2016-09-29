// This message is returned from the Ressource Block schedulers to the 
// MS and BS, answering StreamTransReqs. It tells the MS/BS that all packets 
// scheduled to be send out in the current TTI have been marked in the stream's 
// packet queue.

cplusplus {{
#include "MessageTypes.h"
}}

enum MessageType;
enum MessageDirection;

message StreamTransSched  {
	int src; // The MS from which the packet was originally send
	int messageDirection @enum(MessageDirection); // The transmission direction of this schedule
	kind = MessageType::streamSched;
}
