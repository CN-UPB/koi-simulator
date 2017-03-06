// This message is used by Base and Mobile Stations to send information 
// about SINR values they saw in the previous round.

cplusplus {{
#include "MessageTypes.h"
}}

enum MessageType;
message SINR  {
	int bsId; // The cell ID of the sending mobile station 
	int msId; // The ID of the sending mobile station. -1 if the message originated at a base station
	double up[]; // Expected SINR values for each UP resource block
	double down[]; // Expected SINR values for each DOWN resource block
	int rUp[]; // Expected rate in bit/tti values for each UP resource block
	int rDown[]; // Expected bit/tti values for each DOWN resource block
	int mcsUp[]; // Best possible MCS for each UP resource block
	int mcsDown[]; // Best possible MCS for each DOWN resource block
	kind = MessageType::sinrEst;
}
