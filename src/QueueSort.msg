// This message is send to all BsMac/MsMac to indicate the sort order 
// required for their packet queues.

cplusplus {{
#include "MessageTypes.h"
#include "KoiData_m.h"
#include <functional>
using SortFn = std::function<bool(const KoiData*,const KoiData*)>;
}}

class noncobject SortFn;
enum MessageType;

message QueueSort  {
	SortFn sortfn; // Comparator for the sorting of packets
	kind = MessageType::sortOrder;
}
