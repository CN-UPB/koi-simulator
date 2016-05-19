// This message is used by Mobile stations to inform neighbours that they are
// sending in the current TTI. This information is necessary to compute SINR
// interference values.

cplusplus {{
#include "MessageTypes.h"
#include "TransInfoBs_m.h"
}}

enum MessageType;
class TransInfoBs;

message TransInfoMs extends TransInfoBs {
	int msId; // The sending mobile station
	kind = MessageType::transInfoMs;
}
