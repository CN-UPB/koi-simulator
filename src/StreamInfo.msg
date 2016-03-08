// This message is send to the scheduler by all mobile stations to register 
// their communication streams.

cplusplus {{
#include "MessageTypes.h"
}}

enum MessageType;

message StreamInfo {
	int src;
	int dest;
	double interarrival; // The average time between two messages in this stream
}
