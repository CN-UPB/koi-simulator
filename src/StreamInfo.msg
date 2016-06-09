// This message is send to the scheduler by all mobile stations to register 
// their communication streams.

cplusplus {{
#include "MessageTypes.h"
}}

enum MessageType;

message StreamInfo {
	kind = MessageType::streamInfo;
	unsigned long streamId; // The unique ID of the stream
	int src;
	int dest;
	double interarrival; // The average time between two messages in this stream
}
