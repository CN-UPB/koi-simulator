// This message is used by the StreamingScheduler to send all requests 
// for specific streams to the RBScheduler which handles them.

cplusplus {{
#include <vector>
#include "StreamTransReq_m.h"
#include "MessageTypes.h"
using ReqVec = std::vector<StreamTransReq*>;
}}

class noncobject ReqVec;
enum MessageType;
enum MessageDirection;

message TransReqList {
	ReqVec requests; // All requests for a given stream for the current TTI
	int messageDirection @enum(MessageDirection); // The transmission direction of the requests in this message
	kind = MessageType::transReqList;
}
