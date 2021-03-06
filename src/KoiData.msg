// The message type used as simulated network traffic

cplusplus {{
#include "MessageTypes.h"
}}

enum TrafficType;
enum MessageDirection;
enum MessageType;

packet KoiData {
	int bsId; // The id of the base station servicing the sending MS
	bool critical = false;
	int cqi;
	bool d2d; // Does this package belong to a D2D stream
	simtime_t deadline; // The latest allowed arrival time for the message
	simtime_t interarrival; // The average time between two messages
	int resourceBlock;
	int src; // The ID of the mobile station which created the message
	int dest; // The ID of the target mobile station
	unsigned long streamId; // Unique ID of the stream this packet belongs to
	int trafficType @enum(TrafficType); // The type of traffic, periodic or event based
	double transPower; // The transmission power used for this packet
	int messageDirection @enum(MessageDirection);
	bool scheduled = false; // This packet is scheduled to be send in the current TTI
	simtime_t latestQueueEntry; // Most recent entry time into a queue
	simtime_t totalQueueDelay; // Total timew spend in Queues
	int mcs; // The MCS used in the most recent transmission of this packet
	kind = MessageType::koidata;
}
