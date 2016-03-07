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

message TransReqList {
	ReqVec requests; // All requests for a given stream for the current TTI
	kind = MessageType::transReqList;
}
