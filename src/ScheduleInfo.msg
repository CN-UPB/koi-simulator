// This message is send to all BsMac/MsMac to indicate the sort order 
// required for their packet queues as well as a flag indicating whether 
// schedules for up/down direction are static or dynamic.

cplusplus {{
#include "MessageTypes.h"
#include "KoiData_m.h"
#include <functional>
using SortFn = std::function<bool(const KoiData*,const KoiData*)>;
}}

class noncobject SortFn;
enum MessageType;

message ScheduleInfo  {
	SortFn sortfn; // Comparator for the sorting of packets
	bool upStatic; // Is the schedule for UP transmission generated statically?
	bool downStatic; // Is the schedulue for DOWN transmissions generated statically?
	kind = MessageType::scheduleInfo;
}
