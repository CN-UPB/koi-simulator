// This message type is used as self messages for traffic generators

cplusplus {{
#include "MessageTypes.h"
}}

enum TrafficType;
enum MessageType;

message Traffic {
	kind = MessageType::traffic;
	int partner; // The ID of the mobile station for which to generate traffic
	int trafficType @enum(TrafficType); // The type of traffic, periodic or event based
	unsigned long streamId; // The unique ID identifying the stream this message generates traffic for
}
