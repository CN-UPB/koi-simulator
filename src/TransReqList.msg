// This message is used by the StreamingScheduler to send all requests 
// for specific streams to the RBScheduler which handles them.

cplusplus {{
#include <memory>
#include <vector>
#include <unordered_map>
#include "StreamTransReq_m.h"
#include "SINR_m.h"
#include "MessageTypes.h"
using ReqVec = std::vector<StreamTransReq*>;
using EstMap = std::shared_ptr<std::unordered_map<int,SINR*>>;
}}

class noncobject ReqVec;
class noncobject EstMap;
enum MessageType;
enum MessageDirection;

message TransReqList {
	ReqVec requests; // All requests for a given stream for the current TTI
	EstMap estimates; // All estimates for stations with streams in this request list
	int messageDirection @enum(MessageDirection); // The transmission direction of the requests in this message
	kind = MessageType::transReqList;
}
