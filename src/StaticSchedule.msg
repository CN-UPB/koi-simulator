// This message holds a static schedule for a cell's local stations.

cplusplus {{
#include "MessageTypes.h"
#include <utility>
#include <vector>
using TTISchedule = std::pair<int,std::vector<int>>;
using ScheduleList = std::vector<TTISchedule>;
}}

class noncobject ScheduleList;
enum MessageType;

message StaticSchedule  {
	
	int origin; // The id of the MS for which the schedule is intended, or -1 for the BS
	int scheduleLength; // The number of TTIs after which this schedule is repeated
	ScheduleList schedule; // The Schedule elements, showing which RB to use in which TTI
	kind = MessageType::staticSchedule;
}
