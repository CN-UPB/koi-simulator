// This message is used by Base stations to inform neighbours that they are
// sending in the current TTI. This information is necessary to compute SINR
// interference values.

cplusplus {{
#include "MessageTypes.h"
}}

enum MessageType;
enum MessageDirection;

message TransInfo {
	int bsId; // The sending base station
	double power; // The transmission power of the sending station
	int rb; // The ressource block used for the transmission
	kind = MessageType::transInfo;
	int msId; // The sending mobile station
	int messageDirection @enum(MessageDirection); // The transmission direction of this message
}
